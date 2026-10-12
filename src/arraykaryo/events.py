"""arraykaryo events: near-chromosomal events and runs of homozygosity on every reliable sample, from the corrected LRR,
the BAF and the genotypes; cell fractions; karyotype strings; figures.

Per sample and chromosome (windows.py, segment.py): windows of about 100 probes give three tracks - the median corrected
LRR, the band deviation d at heterozygous calls, the heterozygous rate against the Hardy-Weinberg expectation - which are
standardised by the sample's own window noise and segmented jointly by circular binary segmentation; every boundary is
refined to probe resolution. A segment is typed from its probe-level evidence against the baseline (LRR 0, d 0, heterozygous
rate 1; on the X of an XY sample the sample's own X level): gain (LRR up), loss (LRR down), copy-neutral loss of
heterozygosity (LRR flat, d above 0 at heterozygous calls), or, where the heterozygous rate sits at the error floor with no
split, a run of homozygosity (LRR flat: autozygosity or a constitutional isodisomy; LRR at one copy: a constitutional
deletion, called loss in every cell). A segment covering 90% or more of a chromosome's probed span is whole, of an arm's p
or q, else stretch; floors 2 Mb and 40 probes. Cell fraction from the BAF split, the primary estimate (gain 4d/(1-2d),
loss 4d/(1+2d), copy-neutral LOH 2d), and from the LRR on the cohort's copy scale, the secondary (one copy calibrated from
the male X, three copies from any constitutional trisomy found, else the literature's PennCNV level scaled to the cohort's
compression); both reported with their disagreement. Runs of homozygosity from the genotypes (roh.py: PLINK's sliding
window), the genome fraction in runs per sample, runs of 10 Mb or more as events. The sex chromosomes enter as whole-
chromosome events from the sex step. Every event carries its evidence and a note where the evidence is thin, the sample is
flagged or the region is artefact-prone (genome.artefacts: tagged, never removed). The karyotype string per sample follows
triokaryo's grammar (karyotype.py)."""
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .baf import f_from_d, fit_d, fold
from .cache import Cache, release
from .karyotype import karyotype_string
from .qc import Correction, probe_masks, read_sample_qc, update_sample_qc
from .roh import ROH, runs_of_homozygosity
from .segment import cbs, critical, merge_segments, refine_boundary
from .sex import read_sex
from .util import NA, finite, log, mark_done, num, read_json, robust_sd, write_json, write_tsv
from .windows import save_tracks, tracks_for_block, window_layout

PARAMS = dict(window=100, min_mb=2.0, min_probes=40, lrr_min=0.03, lrr_z=5.0, lrr_alone=0.10, d_min=0.02, d_llr=30.0, n_het_min=20, het_floor=0.2, n_called_min=40, whole_frac=0.9,
              const_f=0.9, alpha=0.01, q_floor=20.0, merge_div=2.0, block=64, l3_default=0.40, l1_literature=-0.66, sample_block=16)
# lrr_min / lrr_z: the smallest LRR shift and its significance (the median's SE from the sample's own noise); lrr_alone: below this |LRR| a
# mosaic gain or loss needs the BAF split too (residual GC waves of a few hundredths survive the correction and must not become events);
# d_min / d_llr: the smallest band deviation and its log-likelihood ratio against d = 0 (the segmentation chose the segment on the same
# data, so the ratio asks more than an unselected test would); het_floor: the heterozygous ratio of a run of homozygosity.
EVENT_COLS = ["sample", "chrom", "start", "end", "size_bp", "span", "bands", "type", "f", "f_lrr", "f_baf", "constitutional", "n_probes", "n_het", "n_called", "lrr", "lrr_se", "d", "d_llr",
              "het_ratio", "source", "artefact", "note", "qc_pass", "sample_flags"]
ROH_COLS = ["sample", "chrom", "start", "end", "size_bp", "n_probes", "n_het", "bands", "lrr", "whole_chromosome", "is_event"]
TYPE_WORD = {"gain": "gain", "loss": "loss", "loh": "copy-neutral LOH", "roh": "run of homozygosity"}


# ------------------------------------------------------------------------------------------------------------- step A

def tracks_chrom(args):
    """Window tracks of every sample on one chromosome (a worker)."""
    out, cache_dir, chrom, w, block = (args[k] for k in ("out", "cache", "chrom", "window", "block"))
    C = Cache(cache_dir)
    corr = Correction(out)
    masks = probe_masks(out, C)[chrom]
    pr = C.probes(chrom)
    keep = masks["keep"] & (~pr.par if chrom == "chrX" else np.ones(pr.n, bool))
    keep_gt = masks["keep_gt"] & (~pr.par if chrom == "chrX" else np.ones(pr.n, bool))
    bounds, mids, spans = window_layout(pr.pos, keep, w)
    shift = np.where(np.isfinite(masks["het_baf_median"]) & (np.abs(masks["het_baf_median"] - 0.5) > 0.01), masks["het_baf_median"] - 0.5, 0.0).astype(np.float32)
    qc = read_sample_qc(out)
    sigma = np.array([max(num(qc[s]["baf_sd"]), 0.02) if finite(num(qc[s]["baf_sd"])) else 0.04 for s in C.samples])
    gt, baf = C.gt(chrom), C.baf(chrom)
    parts = []
    for a in range(0, C.n_samples, block):
        idx = np.arange(a, min(a + block, C.n_samples))
        L = corr.apply(C, chrom, sample_idx=idx)
        parts.append(tracks_for_block(L, np.asarray(baf[idx], dtype=np.float32), np.asarray(gt[idx]), keep, keep_gt, masks["maf"], shift, bounds, sigma[idx]))
        for mm in (gt, baf, C.lrr(chrom)):
            release(mm)
    save_tracks(out, chrom, bounds, mids, spans, parts)
    # the tracks as separate .npy files, so that one sample's row is one contiguous read
    d = os.path.join(out, "events", "windows")
    z = np.load(os.path.join(d, chrom + ".npz"))
    for k in ("lrr_med", "d", "n_het", "het_ratio", "n_called", "n_lrr"):
        np.save(os.path.join(d, "%s.%s.npy" % (chrom, k)), z[k])
    return chrom, bounds.size - 1


# ------------------------------------------------------------------------------------------------------------- step B

class SampleContext:
    """What one sample's calls need, opened once per worker."""

    def __init__(self, out, cache_dir):
        self.out = out
        self.C = Cache(cache_dir)
        self.corr = Correction(out)
        self.masks = probe_masks(self.out, self.C)
        self.qc = read_sample_qc(out)
        self.sex = read_sex(out)
        self.G = self.C.genome
        self.layout = {}
        self.tracks = {}
        for c in self.C.chroms:
            p = os.path.join(out, "events", "windows", c + ".npz")
            if os.path.exists(p):
                z = np.load(p)
                self.layout[c] = dict(bounds=z["bounds"], mids=z["mids"], spans=z["spans"])
                self.tracks[c] = {k: np.load(os.path.join(out, "events", "windows", "%s.%s.npy" % (c, k)), mmap_mode="r") for k in ("lrr_med", "d", "n_het", "het_ratio", "n_called", "n_lrr")}
        self.artefacts = self.G.artefacts()

    def chrom_order(self):
        return self.G.chrom_order()


def sample_noise(ctx, i, chroms):
    """The sample's window-level noise per track over the genome: LRR from first differences (robust to events), d and the
    heterozygous ratio from their robust spread."""
    dl, dd, dh = [], [], []
    for c in chroms:
        t = ctx.tracks[c]
        l = np.asarray(t["lrr_med"][i], dtype=float)
        dl.append(np.diff(l[np.isfinite(l)]))
        dd.append(np.asarray(t["d"][i], dtype=float))
        dh.append(np.asarray(t["het_ratio"][i], dtype=float))
    dl, dd, dh = np.concatenate(dl), np.concatenate(dd), np.concatenate(dh)
    s_l = max(robust_sd(dl) / np.sqrt(2.0), 0.004) if dl.size > 10 else 0.02
    dd = dd[np.isfinite(dd)]
    s_d = max(robust_sd(dd), 0.7 * float(np.mean(dd)) if dd.size else 0.0, 0.004)
    dh = dh[np.isfinite(dh)]
    s_h = max(robust_sd(dh), 0.05)
    return s_l, s_d, s_h


def call_chromosome(ctx, i, chrom, p, noise, x_base=0.0, use_het=True, sample_sd=0.15, baf_sd=0.04):
    """Segments of one chromosome of one sample with their probe-level evidence (not yet typed as events)."""
    t, lay = ctx.tracks[chrom], ctx.layout[chrom]
    bounds, spans = lay["bounds"], lay["spans"]
    n_win = bounds.size - 1
    s_l, s_d, s_h = noise
    lrr_w = np.asarray(t["lrr_med"][i], dtype=float) - x_base
    d_w = np.asarray(t["d"][i], dtype=float)
    h_w = np.asarray(t["het_ratio"][i], dtype=float)
    z = np.stack([np.clip(lrr_w / s_l, -8, 8), np.clip(np.nan_to_num(d_w) / s_d, 0, 10), np.clip((h_w - 1.0) / s_h, -8, 8) if use_het else np.zeros(n_win)])
    bps = cbs(z, alpha=p["alpha"], q_floor=p["q_floor"]) if n_win >= 4 else []
    seg_bounds = merge_segments(z, [0] + bps + [n_win], critical(n_win, p["alpha"], p["q_floor"]) / p["merge_div"]) if bps else [0, n_win]
    # probe-level data of this sample on this chromosome
    pr = ctx.C.probes(chrom)
    m = ctx.masks[chrom]
    keep = m["keep"] & (~pr.par if chrom == "chrX" else True)
    keep_gt = m["keep_gt"] & (~pr.par if chrom == "chrX" else True)
    L = ctx.corr.apply(ctx.C, chrom, sample_idx=[i])[0].astype(float) - x_base
    g = np.asarray(ctx.C.gt(chrom)[i])
    b = np.asarray(ctx.C.baf(chrom)[i], dtype=np.float32)
    shift = np.where(np.isfinite(m["het_baf_median"]) & (np.abs(m["het_baf_median"] - 0.5) > 0.01), m["het_baf_median"] - 0.5, 0.0)
    het = (g == 1) & keep_gt & np.isfinite(b)
    u = fold(b - shift.astype(np.float32))
    exp_het = np.where(np.isfinite(m["maf"]), 2 * m["maf"] * (1 - m["maf"]), np.nan)
    called = (g >= 0) & keep_gt & np.isfinite(exp_het)
    # refine the interior boundaries to probes
    edges = [int(bounds[k]) for k in seg_bounds]
    refined = [0]
    for k in range(1, len(seg_bounds) - 1):
        w = seg_bounds[k]
        lo, hi = int(bounds[max(w - 1, 0)]), int(bounds[min(w + 1, n_win)])
        k_best = refine_boundary([L, np.where(het, u, np.nan)], [keep, het], lo, hi, [sample_sd, baf_sd], pos=pr.pos)
        refined.append(max(refined[-1] + 1, k_best))
    refined.append(pr.n)
    segs = []
    for a, bb in zip(refined[:-1], refined[1:]):
        kp = keep[a:bb] & np.isfinite(L[a:bb])
        n_lrr = int(kp.sum())
        if n_lrr == 0:
            continue
        lrr = float(np.median(L[a:bb][kp]))
        se = 1.2533 * sample_sd / np.sqrt(n_lrr)
        hb = b[a:bb][het[a:bb]] - shift[a:bb][het[a:bb]]
        n_het = int(hb.size)
        d, _, llr = fit_d(hb, sigma=baf_sd) if n_het >= 4 else (np.nan, np.nan, 0.0)
        nc = int(called[a:bb].sum())
        expected = float(np.nansum(exp_het[a:bb][called[a:bb]]))
        observed = int((het[a:bb] & called[a:bb]).sum())
        hr = observed / expected if expected > 2 else np.nan
        pk = pr.pos[a:bb][kp]
        segs.append(dict(chrom=chrom, a=a, b=bb, start=int(pk[0]), end=int(pk[-1]), n_probes=n_lrr, lrr=lrr, lrr_se=se, d=d, d_llr=llr, n_het=n_het, het_ratio=hr, n_called=nc,
                         n_windows=int(seg_bounds[len(segs) + 1] - seg_bounds[len(segs)]) if len(segs) + 1 < len(seg_bounds) else 0))
    return segs, dict(L=L, het=het, u=u, called=called, keep=keep, pos=pr.pos, g=g, b=b, shift=shift, exp_het=exp_het, chrom=chrom)


def type_segment(s, p, L1, L3, use_het=True):
    """The event type of a segment from its evidence, or '' when it is nothing; with f_lrr, f_baf, f and notes."""
    notes = []
    lrr_sig = finite(s["lrr"]) and abs(s["lrr"]) >= max(p["lrr_min"], p["lrr_z"] * s["lrr_se"])
    d_sig = use_het and finite(s["d"]) and s["d"] >= p["d_min"] and s["d_llr"] >= p["d_llr"] and s["n_het"] >= p["n_het_min"]   # no BAF reading on a single X
    het_floor = use_het and finite(s["het_ratio"]) and s["het_ratio"] <= p["het_floor"] and s["n_called"] >= p["n_called_min"]
    kind = ""
    if het_floor:                                                      # no heterozygosity: the few calls left are errors, whatever their BAF says
        d_sig = False
        if lrr_sig and s["lrr"] < 0.5 * L1:
            kind = "loss"
            notes.append("no heterozygous calls and one-copy LRR: a constitutional deletion")
        elif lrr_sig and s["lrr"] >= p["lrr_alone"]:
            kind = "gain"
            notes.append("no heterozygous calls with LRR up: a gain without heterozygosity (an isodisomic duplication, or contamination)")
        else:
            kind = "roh"
    elif lrr_sig and (d_sig or abs(s["lrr"]) >= p["lrr_alone"] or not use_het):
        kind = "gain" if s["lrr"] > 0 else "loss"
    elif d_sig:
        kind = "loh"
    if not kind:
        return "", np.nan, np.nan, np.nan, notes
    f_lrr = np.nan
    if finite(s["lrr"]):
        if kind == "gain":
            f_lrr = float(min(max(s["lrr"] / L3, 0.0), 1.0))
        elif kind == "loss":
            f_lrr = float(min(max(s["lrr"] / L1, 0.0), 1.0))
    f_baf = f_from_d(s["d"], kind if kind != "roh" else "loh") if (d_sig and kind != "roh") else np.nan
    if kind == "roh":
        f = 1.0
    elif kind == "loh":
        f = f_baf
    elif het_floor and kind == "loss":
        f = 1.0
    elif d_sig and s["n_het"] >= p["n_het_min"]:
        f = f_baf
        if finite(f_lrr) and abs(f_lrr - f_baf) > 0.2:
            notes.append("the LRR and the BAF disagree on the cell fraction (%.2f against %.2f)" % (f_lrr, f_baf))
    else:
        f = f_lrr
        if kind in ("gain", "loss") and not d_sig:
            notes.append("cell fraction from the LRR alone (no BAF split: %d heterozygous calls, d %.3f)" % (s["n_het"], s["d"] if finite(s["d"]) else 0.0))
    return kind, f, f_lrr, f_baf, notes


def segment_evidence(pd, a, b, sample_sd, baf_sd):
    """The probe-level evidence of probes [a, b) of one chromosome (the same numbers call_chromosome computes)."""
    L, het, u, called, keep, pos, bb, shift = pd["L"], pd["het"], pd["u"], pd["called"], pd["keep"], pd["pos"], pd["b"], pd["shift"]
    kp = keep[a:b] & np.isfinite(L[a:b])
    n_lrr = int(kp.sum())
    if n_lrr == 0:
        return None
    hb = bb[a:b][het[a:b]] - shift[a:b][het[a:b]]
    d, _, llr = fit_d(hb, sigma=baf_sd) if hb.size >= 4 else (np.nan, np.nan, 0.0)
    exp_het = pd.get("exp_het")
    nc = int(called[a:b].sum())
    expected = float(np.nansum(exp_het[a:b][called[a:b]])) if exp_het is not None else 0.0
    observed = int((het[a:b] & called[a:b]).sum())
    pk = pos[a:b][kp]
    return dict(chrom=pd["chrom"], a=a, b=b, start=int(pk[0]), end=int(pk[-1]), n_probes=n_lrr, lrr=float(np.median(L[a:b][kp])), lrr_se=1.2533 * sample_sd / np.sqrt(n_lrr), d=d, d_llr=llr,
                n_het=int(hb.size), het_ratio=observed / expected if expected > 2 else np.nan, n_called=nc, n_windows=0)


def merge_typed(segs, pd, p, L1, L3, use_het, sample_sd, baf_sd, max_gap_frac=0.15, max_gap_bp=5_000_000):
    """Adjacent segments typed alike (the same kind, levels within reach) are joined, over any short untyped segment between them -
    the segmentation cuts a long event where its windows are noisy, and the pieces must become one event."""
    if len(segs) < 2:
        return segs
    typed = [type_segment(s, p, L1, L3, use_het=use_het)[0] for s in segs]
    changed = True
    while changed and len(segs) > 1:
        changed = False
        for i in range(len(segs) - 1):
            for j in (i + 1, i + 2):
                if j >= len(segs) or not typed[i] or typed[i] != typed[j]:
                    continue
                gap = segs[j]["start"] - segs[i]["end"]
                span = segs[j]["end"] - segs[i]["start"]
                if j == i + 2 and (typed[i + 1] or gap > max(max_gap_bp, max_gap_frac * span)):
                    continue
                se = np.hypot(segs[i]["lrr_se"], segs[j]["lrr_se"])
                if typed[i] in ("gain", "loss") and abs(segs[i]["lrr"] - segs[j]["lrr"]) > max(0.03, 3 * se):
                    continue
                if typed[i] in ("loh",) and finite(segs[i]["d"]) and finite(segs[j]["d"]) and abs(segs[i]["d"] - segs[j]["d"]) > 0.04:
                    continue
                m = segment_evidence(pd, segs[i]["a"], segs[j]["b"], sample_sd, baf_sd)
                if m is None or type_segment(m, p, L1, L3, use_het=use_het)[0] != typed[i]:
                    continue
                segs = segs[:i] + [m] + segs[j + 1:]
                typed = typed[:i] + [typed[i]] + typed[j + 1:]
                changed = True
                break
            if changed:
                break
    return segs


def span_class(G, chrom, start, end, probed, p):
    """whole / p / q / stretch by the share of the chromosome's or arm's probed span covered."""
    c0, c1 = probed["chrom"]
    if c1 > c0 and (min(end, c1) - max(start, c0)) / float(c1 - c0) >= p["whole_frac"]:
        return "whole"
    for arm in ("p", "q"):
        if arm in probed and probed[arm][1] > probed[arm][0]:
            a0, a1 = probed[arm]
            if (min(end, a1) - max(start, a0)) / float(a1 - a0) >= p["whole_frac"] and (end - start) <= 1.25 * (a1 - a0):
                return arm
    return "stretch"


def probed_spans(G, chrom, pos, keep):
    """The probed span of the chromosome and of each arm (first to last kept probe), for the span class."""
    pk = pos[keep]
    out = {"chrom": (int(pk[0]), int(pk[-1])) if pk.size else (0, 0)}
    if G.p_end.get(chrom):
        pa = pk[pk <= G.p_end[chrom]]
        qa = pk[pk > G.p_end[chrom]]
        if pa.size >= 10:
            out["p"] = (int(pa[0]), int(pa[-1]))
        if qa.size >= 10:
            out["q"] = (int(qa[0]), int(qa[-1]))
    return out


def artefact_tags(arts, chrom, start, end):
    tags = []
    for name, c, a, b, why in arts:
        if c != chrom:
            continue
        ov = min(end, b) - max(start, a)
        if ov <= 0:
            continue
        if ov >= 0.5 * (end - start + 1) or (name not in ("pericentromeric", "subtelomeric") and ov >= 0.5 * (b - a + 1)):
            tags.append(name)
    return ",".join(dict.fromkeys(tags))


def call_sample(args):
    """Every event and run of homozygosity of a block of samples (a worker)."""
    out, cache_dir, idxs, p, L1, L3 = (args[k] for k in ("out", "cache", "idxs", "params", "L1", "L3"))
    ctx = SampleContext(out, cache_dir)
    results = []
    for i in idxs:
        results.append(call_one(ctx, i, p, L1, L3))
    return results


def call_one(ctx, i, p, L1, L3):
    C, G = ctx.C, ctx.G
    name = C.samples[i]
    q = ctx.qc[name]
    sx = ctx.sex.get(name, {})
    sample_sd = num(q.get("lrr_sd_diff", q.get("lrr_sd", "")))
    sample_sd = sample_sd if finite(sample_sd) and sample_sd > 0.01 else 0.15
    baf_sd = num(q.get("baf_sd", ""))
    baf_sd = baf_sd if finite(baf_sd) and baf_sd > 0.01 else 0.04
    base = sx.get("base", "")
    xy = base == "XY"
    chroms = [c for c in C.chroms if c in ctx.tracks and c != "chrY"]
    noise = sample_noise(ctx, i, [c for c in chroms if c != "chrX"] or chroms)
    events, rohs = [], []
    roh_bp_total, probed_total = 0, 0
    for chrom in chroms:
        is_x = chrom == "chrX"
        x_base = 0.0
        if is_x and xy:
            l = np.asarray(ctx.tracks[chrom]["lrr_med"][i], dtype=float)
            x_base = float(np.nanmedian(l)) if np.isfinite(l).any() else 0.0
        use_het = not (is_x and xy)
        segs, pd = call_chromosome(ctx, i, chrom, p, noise, x_base=x_base, use_het=use_het, sample_sd=sample_sd, baf_sd=baf_sd)
        probed = probed_spans(G, chrom, pd["pos"], pd["keep"])
        if not is_x:
            probed_total += max(probed["chrom"][1] - probed["chrom"][0], 0)
        min_bp = int(p["min_mb"] * 1e6)
        segs = merge_typed(segs, pd, p, L1, L3, use_het, sample_sd, baf_sd)
        for s in segs:
            kind, f, f_lrr, f_baf, notes = type_segment(s, p, L1, L3, use_het=use_het)
            if not kind:
                continue
            if (s["end"] - s["start"] + 1) < min_bp or s["n_probes"] < p["min_probes"]:
                continue
            span = span_class(G, chrom, s["start"], s["end"], probed, p)
            if is_x and span == "whole" and kind in ("gain", "loss"):
                continue                                          # a whole X gained or lost is the sex step's call
            const = finite(f) and f >= p["const_f"]
            if s["n_het"] < p["n_het_min"] and kind in ("gain", "loss", "loh"):
                notes.append("thin: %d heterozygous calls" % s["n_het"])
            if s["n_probes"] < 100:
                notes.append("thin: %d probes" % s["n_probes"])
            if kind == "roh" and span == "whole":
                notes.append("a whole chromosome without heterozygosity: an isodisomic uniparental disomy, or consanguinity (a trio decides: phase 2)")
            events.append(dict(sample=name, chrom=chrom, start=s["start"], end=s["end"], size_bp=s["end"] - s["start"] + 1, span=span, bands=G.band_range(chrom, s["start"], s["end"]), type=kind,
                               f=f, f_lrr=f_lrr, f_baf=f_baf, constitutional=int(const), n_probes=s["n_probes"], n_het=s["n_het"], n_called=s["n_called"], lrr=s["lrr"], lrr_se=s["lrr_se"],
                               d=s["d"], d_llr=s["d_llr"], het_ratio=s["het_ratio"], source="segmentation", artefact=artefact_tags(ctx.artefacts, chrom, s["start"], s["end"]) or NA,
                               note="; ".join(notes) or NA, qc_pass=q.get("pass", NA), sample_flags=q.get("flags", NA)))
        # runs of homozygosity (autosomes; the X of an XX sample with two X copies - a single X is homozygous by nature)
        x_copies = num(sx.get("x_copies", ""))
        if use_het and not (is_x and finite(x_copies) and x_copies < 1.5):
            runs = runs_of_homozygosity(pd["pos"], pd["het"], pd["called"], ROH)
            for (s0, e0, n_pr, n_h) in runs:
                m = (pd["pos"] >= s0) & (pd["pos"] <= e0) & pd["keep"] & np.isfinite(pd["L"])
                lrr = float(np.median(pd["L"][m])) if m.any() else np.nan
                whole = span_class(G, chrom, s0, e0, probed, p) == "whole"
                is_ev = (e0 - s0 + 1) >= ROH["event_bp"]
                rohs.append(dict(sample=name, chrom=chrom, start=s0, end=e0, size_bp=e0 - s0 + 1, n_probes=n_pr, n_het=n_h, bands=G.band_range(chrom, s0, e0), lrr=lrr, whole_chromosome=int(whole), is_event=int(is_ev)))
                if not is_x:
                    roh_bp_total += e0 - s0 + 1
                if is_ev:
                    overl = [e for e in events if e["chrom"] == chrom and e["type"] in ("roh", "loss") and min(e["end"], e0) - max(e["start"], s0) > 0.5 * min(e["size_bp"], e0 - s0 + 1)]
                    if overl:
                        overl[0]["note"] = (overl[0]["note"] + "; " if overl[0]["note"] != NA else "") + "ROH scan: %.1f Mb, %d probes, %d heterozygous calls" % ((e0 - s0 + 1) / 1e6, n_pr, n_h)
                        continue
                    kind = "loss" if (finite(lrr) and lrr < 0.5 * L1) else "roh"
                    span = span_class(G, chrom, s0, e0, probed, p)
                    note = "from the ROH scan" + ("; one-copy LRR: a constitutional deletion" if kind == "loss" else "")
                    if kind == "roh" and span == "whole":
                        note += "; a whole chromosome without heterozygosity: an isodisomic uniparental disomy, or consanguinity (a trio decides: phase 2)"
                    events.append(dict(sample=name, chrom=chrom, start=s0, end=e0, size_bp=e0 - s0 + 1, span=span, bands=G.band_range(chrom, s0, e0), type=kind, f=1.0, f_lrr=np.nan, f_baf=np.nan,
                                       constitutional=1, n_probes=n_pr, n_het=n_h, n_called=n_pr, lrr=lrr, lrr_se=np.nan, d=np.nan, d_llr=np.nan, het_ratio=np.nan, source="roh_scan",
                                       artefact=artefact_tags(ctx.artefacts, chrom, s0, e0) or NA, note=note, qc_pass=q.get("pass", NA), sample_flags=q.get("flags", NA)))
    # the sex chromosomes as whole-chromosome events
    if sx and sx.get("mosaic", NA) != NA:
        for chrom, kind, f, f_lrr, f_baf in sex_events(sx):
            if chrom not in C.chroms:
                continue
            events.append(dict(sample=name, chrom=chrom, start=1, end=G.length[chrom], size_bp=G.length[chrom], span="whole", bands="", type=kind, f=f, f_lrr=f_lrr, f_baf=f_baf,
                               constitutional=int(finite(f) and f >= p["const_f"]), n_probes=int(num(sx.get("n_x_called" if chrom == "chrX" else "n_y_probes", "0")) or 0), n_het=int(num(sx.get("n_x_het", "0")) or 0) if chrom == "chrX" else 0,
                               n_called=NA, lrr=num(sx.get("x_lrr" if chrom == "chrX" else "y_lrr", "")), lrr_se=np.nan, d=num(sx.get("d_x", "")) if chrom == "chrX" else np.nan, d_llr=np.nan, het_ratio=np.nan,
                               source="sex", artefact=NA, note=sx.get("mosaic", NA), qc_pass=q.get("pass", NA), sample_flags=q.get("flags", NA)))
    f_roh = roh_bp_total / probed_total if probed_total else np.nan
    return dict(sample=name, events=events, rohs=rohs, f_roh=f_roh, complement=sx.get("complement", NA), sex_karyotype=sx.get("karyotype", NA), noise=noise)


def sex_events(sx):
    """(chrom, kind, f, f_lrr, f_baf) rows from a sex.tsv row's mosaic text."""
    out = []
    txt = sx.get("mosaic", "")
    f = num(sx.get("f", "")); f_lrr = num(sx.get("f_lrr", "")); f_baf = num(sx.get("f_baf", ""))
    comp = sx.get("complement", "")
    base = sx.get("base", "")
    if "constitutional" in txt:
        nx, bx = comp.count("X"), base.count("X")
        if nx != bx:
            out.append(("chrX", "gain" if nx > bx else "loss", 1.0, 1.0, f_baf))
        ny, by = comp.count("Y"), base.count("Y")
        if ny != by:
            out.append(("chrY", "gain" if ny > by else "loss", 1.0, 1.0, np.nan))
    for part in txt.split(";"):
        part = part.strip()
        if part.startswith("mosaic") and "of an X" in part:
            out.append(("chrX", "gain" if "gain" in part else "loss", f, f_lrr, f_baf))
        elif "loss of Y" in part and "mosaic" in part:
            fy = _pct(part, f)
            out.append(("chrY", "loss", fy, fy, np.nan))
        elif ("gain of a Y" in part) or ("a Y in" in part):
            fy = _pct(part, f)
            out.append(("chrY", "gain", fy, fy, np.nan))
    return out


def _pct(text, default):
    import re
    m = re.search(r"(\d+)% of cells", text)
    return float(m.group(1)) / 100.0 if m else default


# ------------------------------------------------------------------------------------------------------------- the step

def run_events(out, cache_dir=None, threads=1, params=None, figures="flagged", log_fn=log):
    t0 = time.time()
    cache_dir = cache_dir or os.path.join(out, "cache")
    C = Cache(cache_dir)
    p = dict(PARAMS)
    if params:
        p.update(params)
    edir = os.path.join(out, "events")
    os.makedirs(edir, exist_ok=True)
    qc = read_sample_qc(out)
    sex_summary = read_json(os.path.join(out, "sex", "summary.json")) if os.path.exists(os.path.join(out, "sex", "summary.json")) else {}
    xm = sex_summary.get("x_modes", {})
    L1 = float(xm["one_copy"] - xm["two_copies"]) if xm and xm.get("one_copy") is not None else -0.55
    L3 = p["l3_default"] * (-L1 / -p["l1_literature"])
    calib = dict(one_copy_lrr=L1, one_copy_source="the male X against the female X (sex step)", three_copies_lrr=L3, three_copies_source="PennCNV's +0.40 scaled by the cohort's compression (%.2f / 0.66)" % -L1)
    # step A: tracks
    chroms = [c for c in C.chroms if c != "chrY"]
    jobs = [dict(out=out, cache=cache_dir, chrom=c, window=p["window"], block=p["block"]) for c in chroms]
    t1 = time.time()
    res = _pool(tracks_chrom, jobs, threads)
    log_fn("window tracks: %s (%.1f s)" % (", ".join("%s %d" % r for r in res), time.time() - t1), "events")
    # step B: calls per sample
    passing = [i for i, s in enumerate(C.samples) if qc[s]["hard_fail"] == "0"]
    blocks = [passing[k:k + p["sample_block"]] for k in range(0, len(passing), p["sample_block"])]
    t1 = time.time()
    results = []
    for r in _pool(call_sample, [dict(out=out, cache=cache_dir, idxs=b, params=p, L1=L1, L3=L3) for b in blocks], threads):
        results.extend(r)
    log_fn("calls on %d passing samples (%d not read: hard QC failures) in %.1f s" % (len(passing), C.n_samples - len(passing), time.time() - t1), "events")
    # the three-copy level from any constitutional trisomy found
    tri = [e["lrr"] for r in results for e in r["events"] if e["type"] == "gain" and e["span"] == "whole" and e["source"] == "segmentation" and finite(e["f_baf"]) and e["f_baf"] >= 0.85 and e["chrom"] not in ("chrX", "chrY")]
    if tri:
        L3m = float(np.median(tri))
        calib.update(three_copies_lrr=L3m, three_copies_source="the median LRR of %d constitutional trisom%s found in the cohort" % (len(tri), "y" if len(tri) == 1 else "ies"))
        for r in results:
            for e in r["events"]:
                if e["type"] == "gain" and finite(e["lrr"]) and e["source"] == "segmentation":
                    e["f_lrr"] = float(min(max(e["lrr"] / L3m, 0.0), 1.0))
                    if e.get("note", NA) != NA and "from the LRR alone" in e["note"]:
                        e["f"] = e["f_lrr"]
    # tables
    ev_rows = [e for r in results for e in r["events"]]
    roh_rows = [x for r in results for x in r["rohs"]]
    order = C.genome.chrom_order()
    ev_rows.sort(key=lambda e: (C.samples.index(e["sample"]), order.get(e["chrom"], 99), e["start"]))
    roh_rows.sort(key=lambda e: (C.samples.index(e["sample"]), order.get(e["chrom"], 99), e["start"]))
    write_tsv(os.path.join(edir, "events.all.tsv"), EVENT_COLS, ev_rows, nd=4)
    write_tsv(os.path.join(edir, "roh.all.tsv"), ROH_COLS, roh_rows, nd=4)
    kar_rows = []
    n_ev = {}
    for r in results:
        evs = [e for e in r["events"]]
        kar = karyotype_string([dict(e, constitutional=bool(e["constitutional"])) for e in evs], r["complement"] if r["complement"] != NA else "", r["sex_karyotype"], order)
        n_ev[r["sample"]] = len(evs)
        kar_rows.append(dict(sample=r["sample"], complement=r["complement"], karyotype=kar, n_events=len(evs), n_roh=len(r["rohs"]), f_roh=r["f_roh"],
                             noise_lrr=r["noise"][0], noise_d=r["noise"][1], noise_het=r["noise"][2], qc_pass=qc[r["sample"]]["pass"], flags=qc[r["sample"]]["flags"]))
    for s in C.samples:
        if s not in n_ev:
            kar_rows.append(dict(sample=s, complement=NA, karyotype=NA, n_events=NA, n_roh=NA, f_roh=np.nan, noise_lrr=np.nan, noise_d=np.nan, noise_het=np.nan, qc_pass=qc[s]["pass"], flags=qc[s]["flags"]))
    kar_rows.sort(key=lambda r: C.samples.index(r["sample"]))
    write_tsv(os.path.join(edir, "karyotypes.tsv"), ["sample", "complement", "karyotype", "n_events", "n_roh", "f_roh", "noise_lrr", "noise_d", "noise_het", "qc_pass", "flags"], kar_rows, nd=4)
    # the event-count flag
    counts = np.array([n_ev[s] for s in C.samples if s in n_ev], dtype=float)
    updates = {}
    cut = None
    if counts.size >= 10:
        cut = float(np.median(counts) + 3.0 * max(robust_sd(counts), 1.0))
        for s, n in n_ev.items():
            if n > cut:
                updates[s] = {"add_flag": "event_count", "n_events": n}
            else:
                updates[s] = {"n_events": n}
    else:
        updates = {s: {"n_events": n} for s, n in n_ev.items()}
    update_sample_qc(out, updates, new_cols=("n_events",))
    by_type = {}
    for e in ev_rows:
        k = "%s %s%s" % (e["type"], e["span"], "" if e["constitutional"] else " mosaic")
        by_type[k] = by_type.get(k, 0) + 1
    summary = dict(samples=C.n_samples, called=len(passing), events=len(ev_rows), samples_with_events=len({e["sample"] for e in ev_rows}), by_type=by_type, roh_runs=len(roh_rows),
                   roh_events=sum(1 for x in roh_rows if x["is_event"]), calibration=calib, event_count_cut=cut, params=p, seconds=time.time() - t0)
    write_json(summary, os.path.join(edir, "summary.json"))
    log_fn("%d events on %d samples (%s); %d runs of homozygosity; one-copy LRR %.3f, three-copy %.3f (%s)" % (
        len(ev_rows), summary["samples_with_events"], ", ".join("%s %d" % kv for kv in sorted(by_type.items())) or "none", len(roh_rows), calib["one_copy_lrr"], calib["three_copies_lrr"], calib["three_copies_source"]), "events")
    # figures
    if figures != "none":
        from .plots import draw_sample_figures
        want = set(e["sample"] for e in ev_rows) | {s for s in C.samples if qc[s]["flags"] != NA} if figures == "flagged" else set(C.samples)
        want = [s for s in C.samples if s in want]
        t1 = time.time()
        n = draw_sample_figures(out, cache_dir, want, threads, log_fn)
        log_fn("figures for %d samples in %.1f s -> events/figures/" % (n, time.time() - t1), "events")
    mark_done(out, "events", dict(seconds=time.time() - t0))
    return summary


def _pool(fn, jobs, threads):
    if threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            return list(ex.map(fn, jobs))
    return [fn(j) for j in jobs]


def read_events(out):
    from .util import read_tsv
    p = os.path.join(out, "events", "events.all.tsv")
    return read_tsv(p) if os.path.exists(p) else []


def read_rohs(out):
    from .util import read_tsv
    p = os.path.join(out, "events", "roh.all.tsv")
    return read_tsv(p) if os.path.exists(p) else []


def read_karyotypes(out):
    from .util import read_tsv
    p = os.path.join(out, "events", "karyotypes.tsv")
    return {r["sample"]: r for r in read_tsv(p)} if os.path.exists(p) else {}
