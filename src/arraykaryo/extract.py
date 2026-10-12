"""arraykaryo extract: one pass over the export into the cache (cache.py), resumable per chromosome, chromosomes in parallel.

Two routes, timed on the mock at cohort scale (docs/METHODS.md records the numbers; `--route` picks one, `auto` the faster
available): (a) `bcftools query` streaming one chromosome as text into a vectorised numpy parser (the tabs of a block of
lines are found at once; a genotype is read from two bytes; the floats are parsed in bulk, digit columns against powers of
ten); (b) cyvcf2, whose per-record numpy arrays (gt_types, format('LRR'), format('BAF')) are filled in C. pysam's
per-sample record access is far slower and is kept for `inspect` only. Records are written in blocks (BLOCK records at a
time, transposed into the sample-major arrays). Every record of a chromosome is kept, in file order: multiallelic and indel
records and intensity-only probes are flagged in the probe table, not dropped (their LRR is usable)."""
import hashlib
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pysam

from . import __version__
from .cache import finish_arrays, open_arrays
from .genome import genome_for, normalize_contig
from .inspect import REQUIRED, index_counts
from .util import log, write_json, write_tsv

BLOCK = 8192                 # records per block written into the arrays (one positioned write per sample per block)
CHUNK_BYTES = 1 << 24        # bytes of text read from bcftools at a time (route a)
INFO_FIELDS = ("INTENSITY_ONLY", "GC", "GenTrain_Score", "Cluster_Sep", "ALLELE_A", "ALLELE_B")
POW10 = 10.0 ** np.arange(-40, 41).astype(np.float64)


def header_info(path, samples=None):
    """What the extraction needs from the header: samples (the subset in export order), the contig names and lengths, the INFO
    fields of interest present, the md5 of the header text, and the assembly."""
    vf = pysam.VariantFile(path)
    names = list(vf.header.samples)
    if samples is not None:
        want = set(samples)
        missing = sorted(want - set(names))
        if missing:
            raise SystemExit("samples not in the export: " + ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))
        names = [s for s in names if s in want]
    contigs = {c.name: c.length for c in vf.header.contigs.values()}
    formats = list(vf.header.formats)
    missing = [f for f in REQUIRED if f not in formats]
    if missing:
        raise SystemExit("the export lacks FORMAT/%s; arraykaryo needs GT, LRR and BAF (arraykaryo inspect says what it carries)" % ", FORMAT/".join(missing))
    infos = [k for k in INFO_FIELDS if k in vf.header.info]
    md5 = hashlib.md5(str(vf.header).encode()).hexdigest()
    vf.close()
    return dict(samples=names, contigs=contigs, infos=infos, header_md5=md5)


# ----------------------------------------------------------------------------------------------------------- route (a)

def parse_floats(a, starts, ends):
    """Decimal tokens (as bcftools prints floats: -0.0480331, 1.5e-05, nan, .) at a[starts:ends] -> float64 array. Vectorised:
    the tokens are laid out as a (n, width) byte matrix; each digit column contributes digit x 10^(its place relative to the
    point); an exponent part scales the mantissa; nan, -nan and '.' give NaN."""
    n = starts.size
    if n == 0:
        return np.zeros(0)
    L = ends - starts
    W = int(L.max()) if n else 0
    cols = np.arange(W)
    idx = starts[:, None] + cols[None, :]
    m = a[np.minimum(idx, a.size - 1)]
    m = np.where(cols[None, :] < L[:, None], m, 32).astype(np.int16)
    neg = m[:, 0] == 45
    isdig = (m >= 48) & (m <= 57)
    isdot = m == 46
    ise = (m == 101) | (m == 69)
    has_dot = isdot.any(axis=1)
    dot = np.where(has_dot, isdot.argmax(axis=1), -1)
    has_e = ise.any(axis=1)
    epos = np.where(has_e, ise.argmax(axis=1), W)
    mant_end = np.where(has_dot, dot, epos)                       # the place of the units digit: just before the point, or the token end
    mant_end = np.where(has_dot, dot, np.minimum(epos, L))
    # the place value of column c: units at mant_end-1, tenths at dot+1 ...
    place = np.where(cols[None, :] < mant_end[:, None], mant_end[:, None] - 1 - cols[None, :], mant_end[:, None] - cols[None, :])
    place = np.clip(place, -40, 40)
    use = isdig & (cols[None, :] < epos[:, None])
    v = np.zeros(n, dtype=np.float64)
    for c in range(W):
        u = use[:, c]
        if u.any():
            v += np.where(u, (m[:, c] - 48) * POW10[place[:, c] + 40], 0.0)
    if has_e.any():
        rows = np.flatnonzero(has_e)
        e = np.zeros(rows.size, dtype=np.int64)
        esign = np.ones(rows.size, dtype=np.int64)
        for c in range(W):
            mc = m[rows, c]
            after = c > epos[rows]
            esign = np.where(after & (mc == 45) & (c == epos[rows] + 1), -1, esign)
            d = after & (mc >= 48) & (mc <= 57)
            e = np.where(d, e * 10 + (mc - 48), e)
        v[rows] *= 10.0 ** np.clip(esign * e, -40, 40)
    v = np.where(neg, -v, v)
    bad = ~isdig.any(axis=1) | (m[:, 0] == 110) | (m[:, 1] == 110)          # nan, -nan, '.', empty
    v[bad] = np.nan
    return v


def parse_block(buf, n_fixed, n_samples):
    """A block of complete lines of `bcftools query` output -> (fixed fields as lists of str per column, gt int8 (lines x samples),
    lrr float32, baf float32). Every line carries n_fixed + 3 x n_samples tab-separated fields."""
    a = np.frombuffer(buf, dtype=np.uint8)
    nl = np.flatnonzero(a == 10)
    n = nl.size
    per = n_fixed + 3 * n_samples - 1
    tabs = np.flatnonzero(a == 9)
    if tabs.size != n * per:
        raise ValueError("unexpected field count in the bcftools output (%d tabs for %d lines, %d expected per line): a sample field is missing" % (tabs.size, n, per))
    tabs = tabs.reshape(n, per)
    starts = np.concatenate([np.r_[0, nl[:-1] + 1][:, None], tabs + 1], axis=1)       # field starts (n x fields)
    ends = np.concatenate([tabs, nl[:, None]], axis=1)                               # field ends (exclusive)
    fixed = [[a[starts[i, k]:ends[i, k]].tobytes().decode() for i in range(n)] for k in range(n_fixed)]
    s0 = starts[:, n_fixed::3]
    b0 = a[s0].astype(np.int16) - 48
    b1 = a[np.minimum(s0 + 1, a.size - 1)]
    b2 = a[np.minimum(s0 + 2, a.size - 1)].astype(np.int16) - 48
    haploid = (b1 == 9) | (b1 == 10)
    gt = np.where(haploid, np.where(b0 >= 0, 2 * np.minimum(b0, 1), -1), np.where((b0 < 0) | (b2 < 0), -1, np.minimum(b0, 1) + np.minimum(b2, 1)))
    gt = gt.astype(np.int8)
    lrr = parse_floats(a, starts[:, n_fixed + 1::3].ravel(), ends[:, n_fixed + 1::3].ravel()).reshape(n, n_samples).astype(np.float32)
    baf = parse_floats(a, starts[:, n_fixed + 2::3].ravel(), ends[:, n_fixed + 2::3].ravel()).reshape(n, n_samples).astype(np.float32)
    return fixed, gt, lrr, baf


def query_format(infos):
    fields = ["%POS", "%ID", "%REF", "%ALT"] + ["%INFO/" + k for k in infos]
    return "\t".join(fields) + "[\t%GT\t%LRR\t%BAF]\n", len(fields)


def iter_bcftools(path, contig, infos, samples):
    """Blocks of complete lines from `bcftools query -r contig`."""
    fmt, n_fixed = query_format(infos)
    cmd = [shutil.which("bcftools") or "bcftools", "query", "-r", contig, "-f", fmt]
    if samples is not None:
        cmd += ["-s", ",".join(samples)]
    cmd.append(path)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    rest = b""
    while True:
        chunk = proc.stdout.read(CHUNK_BYTES)
        if not chunk:
            break
        buf = rest + chunk
        cut = buf.rfind(b"\n")
        if cut < 0:
            rest = buf
            continue
        yield buf[:cut + 1], n_fixed
        rest = buf[cut + 1:]
    if rest:
        yield rest + (b"" if rest.endswith(b"\n") else b"\n"), n_fixed
    err = proc.stderr.read().decode()
    if proc.wait() != 0:
        raise RuntimeError("bcftools query failed on %s: %s" % (contig, err.strip()))


def _info_value(tok, kind):
    if kind == "flag":
        return tok not in (".", "0", "")
    try:
        return float(tok)
    except ValueError:
        return np.nan


def extract_bcftools(path, contig, infos, samples, gt, lrr, baf, n_probes):
    """Route (a) for one chromosome into the open arrays; returns the probe rows."""
    rows = dict(pos=[], id=[], ref=[], alt=[], info={k: [] for k in infos})
    i0 = 0
    n_samples = gt.shape[0]
    for buf, n_fixed in iter_bcftools(path, contig, infos, samples):
        fixed, g, l, b = parse_block(buf, n_fixed, n_samples)
        n = g.shape[0]
        if i0 + n > n_probes:
            raise RuntimeError("%s: more records than the index reports (%d+)" % (contig, i0 + n))
        gt.write(i0, g.T)
        lrr.write(i0, l.T)
        baf.write(i0, b.T)
        rows["pos"].extend(fixed[0]), rows["id"].extend(fixed[1]), rows["ref"].extend(fixed[2]), rows["alt"].extend(fixed[3])
        for k, col in zip(infos, fixed[4:]):
            rows["info"][k].extend(col)
        i0 += n
    return rows, i0


# ----------------------------------------------------------------------------------------------------------- route (b)

def extract_cyvcf2(path, contig, infos, samples, gt, lrr, baf, n_probes):
    import cyvcf2
    vf = cyvcf2.VCF(path, gts012=True, lazy=True)
    if samples is not None:
        vf.set_samples(samples)
    rows = dict(pos=[], id=[], ref=[], alt=[], info={k: [] for k in infos})
    n_samples = gt.shape[0]
    i0 = 0
    bg = np.empty((BLOCK, n_samples), dtype=np.int8)
    bl = np.empty((BLOCK, n_samples), dtype=np.float32)
    bb = np.empty((BLOCK, n_samples), dtype=np.float32)
    k = 0

    def flush(k):
        gt.write(i0, bg[:k].T)
        lrr.write(i0, bl[:k].T)
        baf.write(i0, bb[:k].T)

    for v in vf(contig):
        g = v.gt_types
        bg[k] = np.where(g == 3, -1, g)
        fl = v.format("LRR")
        bl[k] = fl[:, 0] if fl is not None else np.nan
        fb = v.format("BAF")
        bb[k] = fb[:, 0] if fb is not None else np.nan
        rows["pos"].append(str(v.POS)), rows["id"].append(v.ID or "."), rows["ref"].append(v.REF), rows["alt"].append(",".join(v.ALT) if v.ALT else ".")
        for key in infos:
            val = v.INFO.get(key)
            rows["info"][key].append("1" if val is True else ("." if val is None or val is False else str(val)))
        k += 1
        if k == BLOCK:
            flush(k)
            i0 += k
            k = 0
    if k:
        flush(k)
        i0 += k
    return rows, i0


# ----------------------------------------------------------------------------------------------------------- the pass

def extract_chrom(args):
    """One chromosome into the cache (a worker; args is a dict). Returns (chrom, n, seconds)."""
    path, out, chrom, contig, infos, samples, n_probes, route, baf_dtype = (args[k] for k in ("path", "out", "chrom", "contig", "infos", "samples", "n_probes", "route", "baf_dtype"))
    t0 = time.time()
    n_samples = len(samples)
    arrays = open_arrays(out, chrom, n_samples, n_probes, np.dtype(baf_dtype))
    fn = extract_cyvcf2 if route == "cyvcf2" else extract_bcftools
    rows, n = fn(path, contig, infos, samples, *arrays, n_probes)
    if n != n_probes:
        raise RuntimeError("%s: %d records read, %d in the index" % (chrom, n, n_probes))
    finish_arrays(out, chrom, arrays)
    np.savez(os.path.join(out, "probes.%s.npz" % chrom), pos=np.asarray(rows["pos"]), id=np.asarray(rows["id"]), ref=np.asarray(rows["ref"]), alt=np.asarray(rows["alt"]),
             **{"info_" + k: np.asarray(v) for k, v in rows["info"].items()})
    write_json(dict(chrom=chrom, n=n, seconds=time.time() - t0, route=route), os.path.join(out, "done", chrom + ".json"))
    return chrom, n, time.time() - t0


def choose_route(route):
    """cyvcf2 when installed (five times faster than bcftools query plus the numpy parser on the mock at cohort scale: METHODS), else bcftools."""
    if route == "auto":
        try:
            import cyvcf2  # noqa: F401
            return "cyvcf2"
        except ImportError:
            return "bcftools"
    if route == "cyvcf2":
        import cyvcf2  # noqa: F401
    return route


def extract(path, out, threads=1, route="auto", baf_dtype="float16", redo=False, samples=None, genome="", chroms=None, log_fn=log):
    """The pass. `samples`: a subset (export order kept); `chroms`: a subset of chromosomes (canonical names)."""
    t0 = time.time()
    os.makedirs(os.path.join(out, "done"), exist_ok=True)
    h = header_info(path, samples)
    G = genome_for(h["contigs"], genome)
    names = {normalize_contig(c): c for c in h["contigs"]}                   # canonical -> as in the file
    counts = index_counts(path)
    if not counts:
        raise SystemExit("no index statistics for %s: it needs a .csi/.tbi index and bcftools on PATH" % path)
    counts = {normalize_contig(c): n for c, n in counts.items()}
    todo = [c for c in G.chroms if c in names and counts.get(c, 0) > 0 and (chroms is None or c in chroms)]
    skipped = [c for c in counts if c not in G.length]
    if skipped:
        log_fn("contigs not in %s, skipped: %s" % (G.label, ", ".join(skipped[:8]) + (" ..." if len(skipped) > 8 else "")), "extract")
    route = choose_route(route)
    with open(os.path.join(out, "samples.txt"), "w") as fh:
        fh.write("\n".join(h["samples"]) + "\n")
    pending = []
    for c in todo:
        if not redo and os.path.exists(os.path.join(out, "done", c + ".json")) and os.path.exists(os.path.join(out, "gt.%s.npy" % c)):
            continue
        pending.append(dict(path=path, out=out, chrom=c, contig=names[c], infos=h["infos"], samples=h["samples"], n_probes=counts[c], route=route, baf_dtype=baf_dtype))
    pending.sort(key=lambda a: -a["n_probes"])
    log_fn("%d samples, %d chromosomes (%d to do), %s, route %s, %d thread(s)" % (len(h["samples"]), len(todo), len(pending), G.label, route, threads), "extract")
    if pending:
        if threads > 1 and len(pending) > 1:
            with ProcessPoolExecutor(max_workers=threads) as ex:
                futs = [ex.submit(extract_chrom, a) for a in pending]
                for f in as_completed(futs):
                    c, n, s = f.result()
                    log_fn("%s: %d records in %.1f s" % (c, n, s), "extract")
        else:
            for a in pending:
                c, n, s = extract_chrom(a)
                log_fn("%s: %d records in %.1f s" % (c, n, s), "extract")
    assemble_probes(out, todo, G, h["infos"])
    st = os.stat(path)
    write_json(dict(export=path, size=st.st_size, mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)), header_md5=h["header_md5"], samples=len(h["samples"]),
                    chroms={c: int(counts[c]) for c in todo}, assembly=G.name, genome_label=G.label, genome_table=genome if (genome and os.path.exists(str(genome))) else "",
                    contig_names={c: names[c] for c in todo}, info_fields=h["infos"], baf_dtype=baf_dtype, route=route, threads=threads, version=__version__,
                    seconds=time.time() - t0, created=time.strftime("%Y-%m-%d %H:%M:%S")), os.path.join(out, "manifest.json"))
    log_fn("cache %s: %d chromosomes, %d probes, %.1f s" % (out, len(todo), sum(counts[c] for c in todo), time.time() - t0), "extract")
    return out


def assemble_probes(out, chroms, G, infos):
    """probes.tsv and probes.npz over all chromosomes from the per-chromosome npz files."""
    cols = dict(pos=[], snv=[], intensity_only=[], par=[], gc=[], gentrain=[], cluster_sep=[])
    tsv = []
    for c in chroms:
        z = np.load(os.path.join(out, "probes.%s.npz" % c))
        pos = z["pos"].astype(np.int64)
        ref, alt = z["ref"], z["alt"]
        rl = np.char.str_len(ref.astype(str))
        al = np.char.str_len(alt.astype(str))
        snv = (rl == 1) & (al == 1) & (alt.astype(str) != ".")
        io = np.zeros(pos.size, dtype=bool)
        if "info_INTENSITY_ONLY" in z.files:
            io = ~np.isin(z["info_INTENSITY_ONLY"].astype(str), (".", "0", ""))
        io |= alt.astype(str) == "."
        gc = _floats(z, "info_GC", pos.size)
        gentrain = _floats(z, "info_GenTrain_Score", pos.size)
        csep = _floats(z, "info_Cluster_Sep", pos.size)
        par = G.par_mask(c, pos)
        cols["pos"].append(pos), cols["snv"].append(snv), cols["intensity_only"].append(io), cols["par"].append(par)
        cols["gc"].append(gc.astype(np.float32)), cols["gentrain"].append(gentrain.astype(np.float32)), cols["cluster_sep"].append(csep.astype(np.float32))
        aa = z["info_ALLELE_A"].astype(str) if "info_ALLELE_A" in z.files else np.full(pos.size, "NA")
        ab = z["info_ALLELE_B"].astype(str) if "info_ALLELE_B" in z.files else np.full(pos.size, "NA")
        ids = z["id"].astype(str)
        for i in range(pos.size):
            tsv.append((c, int(pos[i]), ids[i], str(ref[i]), str(alt[i]), int(snv[i]), int(io[i]), int(par[i]), gc[i], gentrain[i], csep[i], aa[i].replace(".", "NA"), ab[i].replace(".", "NA")))
    np.savez(os.path.join(out, "probes.npz"), **{k: np.concatenate(v) if v else np.zeros(0) for k, v in cols.items()})
    from .cache import PROBE_COLS
    write_tsv(os.path.join(out, "probes.tsv"), list(PROBE_COLS), tsv, nd=5)


def _floats(z, key, n):
    if key not in z.files:
        return np.full(n, np.nan)
    v = z[key].astype(str)
    out = np.full(n, np.nan)
    ok = ~np.isin(v, (".", "", "NA"))
    out[ok] = v[ok].astype(float)
    return out
