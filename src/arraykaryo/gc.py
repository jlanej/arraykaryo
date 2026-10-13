"""The GC content around each probe from a reference FASTA, for an export whose INFO/GC is absent.

`bcftools +gtc2vcf` writes INFO/GC from the reference over +-200 bp of the variant (its default --gc-window-size); Illumina's
own GTCtoVCF and most other converters write none. arraykaryo reads INFO/GC when it is there and, given --fasta, computes the
same quantity here for the chromosomes whose GC is wholly missing: (G+C)/(A+C+G+T) over the 401 bases centred on the probe,
NaN where the contig is absent from the FASTA or the window holds no called base. Without either source the QC step skips the
GC-wave regression and says so (the cohort reference per probe and the per-sample offset still apply, the waviness factor is
unsigned); a probe without GC in a cohort that has GC gets no wave term of its own.
"""
import os

import numpy as np

from .util import log, read_json, write_json

HALF = 200                                   # gtc2vcf's default --gc-window-size
CHUNK = 20_000_000                           # bases of sequence held at a time


def contig_name(fa, chrom):
    """The FASTA's own name for a normalised contig ('chr1' or '1')."""
    refs = set(fa.references)
    for name in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom):
        if name in refs:
            return name
    return None


def gc_from_fasta(fasta, chrom, pos, half=HALF, chunk=CHUNK):
    """GC ratio over +-half bp around each 1-based position, from an indexed FASTA; NaN where it cannot be read."""
    import pysam
    pos = np.asarray(pos, dtype=np.int64)
    out = np.full(pos.size, np.nan, dtype=np.float32)
    with pysam.FastaFile(fasta) as fa:
        name = contig_name(fa, chrom)
        if name is None or pos.size == 0:
            return out
        length = fa.get_reference_length(name)
        order = np.argsort(pos, kind="stable")
        p = pos[order]
        for start in range(0, length, chunk):
            end = min(length, start + chunk)
            i0, i1 = np.searchsorted(p, start + 1, side="left"), np.searchsorted(p, end, side="right")  # positions in start+1..end
            if i1 <= i0:
                continue
            a, b = max(0, start - half), min(length, end + half)
            seq = np.frombuffer(fa.fetch(name, a, b).upper().encode("ascii"), dtype=np.uint8)
            is_gc = ((seq == 71) | (seq == 67)).astype(np.int32)                                       # G, C
            is_base = (is_gc | (seq == 65) | (seq == 84)).astype(np.int32)                             # A, T as well
            is_gc, is_base = np.append(is_gc, 0), np.append(is_base, 0)                                 # a window may end at the sequence's end
            lo = np.clip(p[i0:i1] - 1 - half - a, 0, seq.size)                                         # 0-based [lo, hi)
            hi = np.clip(p[i0:i1] + half - a, 0, seq.size)
            ok = hi > lo
            if not ok.any():
                continue
            idx = np.stack([lo[ok], hi[ok]], axis=1).ravel()
            g = np.add.reduceat(is_gc, idx)[::2]
            n = np.add.reduceat(is_base, idx)[::2]
            with np.errstate(invalid="ignore", divide="ignore"):
                out[order[i0:i1][ok]] = np.where(n > 0, g / np.maximum(n, 1), np.nan).astype(np.float32)
    return out


def gc_source_of(cache):
    """'INFO/GC' when the cache's probes carry GC, 'none' when no probe does, 'partial' when whole chromosomes lack it."""
    have = [np.isfinite(cache.probes(c).gc).any() for c in cache.chroms]
    return "INFO/GC" if all(have) else ("none" if not any(have) else "partial")


def fill_gc(cache, fasta, log_fn=log, force=False):
    """Fill the cache's probe GC from the FASTA for every chromosome whose GC is wholly missing (`force`: every chromosome);
    rewrites probes.npz and probes.tsv and records the source in the manifest. Returns the chromosomes filled."""
    npz = dict(cache.npz)
    gc = npz["gc"].astype(np.float32).copy()
    filled = []
    for c in cache.chroms:
        sl = cache._offsets[c]
        if not force and np.isfinite(gc[sl]).any():
            continue
        v = gc_from_fasta(fasta, c, npz["pos"][sl])
        if np.isfinite(v).any():
            gc[sl] = v
            filled.append(c)
        else:
            log_fn("%s is not in the FASTA (or its windows hold no called base): no GC for its probes" % c, "gc")
    if filled:
        npz["gc"] = gc
        np.savez(os.path.join(cache.path, "probes.npz"), **npz)
        _rewrite_tsv_gc(cache, gc)
        cache._npz = None
        log_fn("GC computed from the FASTA (gtc2vcf's window, +-%d bp) for %d probes on %s" % (HALF, sum(cache.n_probes[c] for c in filled), ", ".join(filled)), "gc")
    m = read_json(os.path.join(cache.path, "manifest.json"))
    m["gc_source"] = gc_source_of(cache) if not filled else ("fasta" if len(filled) == len(cache.chroms) else "INFO/GC + fasta")
    m["gc_fasta"] = os.path.basename(fasta) if filled else m.get("gc_fasta", "")
    write_json(m, os.path.join(cache.path, "manifest.json"))
    return filled


def _rewrite_tsv_gc(cache, gc):
    """The gc column of probes.tsv, streamed (the table of a large array has millions of rows)."""
    path = os.path.join(cache.path, "probes.tsv")
    tmp = path + ".tmp"
    with open(path) as fh, open(tmp, "w") as out:
        header = fh.readline()
        cols = header.rstrip("\n").split("\t")
        j = cols.index("gc")
        out.write(header)
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            f[j] = "NA" if not np.isfinite(gc[i]) else "%.5f" % gc[i]
            out.write("\t".join(f) + "\n")
    os.replace(tmp, path)
