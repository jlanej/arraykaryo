"""arraykaryo qc: probe QC, sample QC, the corrections the calls need (the GC wave per sample, the cohort reference per
probe) and a verdict per sample - pass, or flagged with the reason; a hard failure marks a sample "not to be read". Flags,
not drops: every sample stays in the tables.

Probe QC (cohort-wide, qc/probes.tsv): intensity-only probes keep their LRR and lose BAF and GT; duplicate positions;
non-SNV records; a genotype call rate across samples under 0.95; a per-probe LRR median beyond +-0.3 or an LRR robust SD in
the top 1% (mis-normalised, noisy); a heterozygous-call BAF whose cohort median departs from 0.5 by more than 0.05 (a
shifted cluster); Hardy-Weinberg failure at p < 1e-10 (provisionally over the passing samples; the pedigree step repeats it
over the unrelated ones); GenTrain_Score under 0.6 or Cluster_Sep under 0.3 when present. PAR probes are kept and tagged.
Two verdicts per probe: `keep` (the LRR is usable) and `keep_gt` (the BAF and the genotype are usable).

Sample QC (qc/samples.tsv): call rate (flag under 0.97: Illumina's and MoChA's practice); autosomal LRR SD after correction
(flag above 0.20 for mosaic calling, a hard fail above 0.30: PennCNV's filter default); BAF SD at heterozygous calls (flag
above 0.06; taken as the median over the autosomes of the per-chromosome SD so that one chromosome's split bands do not
flag a sample); BAF drift, the share of homozygous-called probes whose BAF lies in (0.2, 0.25) or (0.75, 0.8) (flag above
0.002; PennCNV's drift counts every probe and filters at 0.01); the wave factor as PennCNV defines it - the median absolute
deviation of the 1-Mb window medians of the raw autosomal LRR from their overall median, signed by the correlation of the
window medians with the window GC (PennCNV signs it against a reference wave) - flag |WF| above 0.05 (PennCNV's qcwf;
detect_cnv warns at 0.04); the LRR mean (flag beyond +-0.1); later the sex check, a duplicate or an unexpected relative, and
an event count beyond 3 robust SD above the median.

Corrections, per sample, on the autosomes and the X: (1) the GC wave - LRR regressed on the probe's GC and on the GC of its
1-Mb window with a quadratic in each (Diskin 2008; MoChA's --adjust-BAF-LRR is the modern reference), fitted on the
autosomes by least squares on the LRR clipped to +-1, and subtracted; (2) the cohort reference - the median corrected LRR
per probe over the passing samples, per batch when a batch holds 20 or more passing samples, else cohort-wide, subtracted
(the array's shared structure, as triokaryo's panel removes sequencing's); on the X the reference comes from the samples
without a Y so that one X reads as one copy; the Y keeps its GC-corrected raw level against the male cluster; (3) BAF:
clipped to [0, 1] at use, and the heterozygous-call BAF re-centred per probe by the cohort's median when its shift exceeds
0.01 (MoChA's cluster adjustment). The corrected LRR is never written: qc/correction.npz holds the parameters and
Correction.apply reproduces it from the cache for any block."""
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .cache import Cache, release
from .util import NA, log, mark_done, read_json, read_sample_sheet, robust_sd, write_json, write_tsv

BLOCK = 4096
THRESHOLDS = dict(
    probe_call_rate=(0.95, "genotype call rate across samples; PennCNV's and MoChA's practice (variant call rate 0.97 in MoChA's plots, 0.95 for keeping a probe)"),
    probe_lrr_median=(0.3, "per-probe cohort median of the raw LRR beyond +-0.3: a mis-normalised probe"),
    probe_lrr_sd_top=(0.01, "per-probe LRR robust SD in the top 1% of probes: a noisy probe"),
    probe_baf_shift=(0.05, "cohort median BAF at heterozygous calls more than 0.05 from 0.5: a shifted cluster (MoChA re-centres shifts above 0.01; beyond 0.05 the probe is dropped)"),
    probe_hwe_p=(1e-10, "Hardy-Weinberg exact-test p among unrelated passing samples (PLINK's --hwe practice for arrays: 1e-10 to 1e-6)"),
    probe_gentrain=(0.6, "Illumina GenTrain score; Illumina's guidance treats probes under 0.6 as poorly clustered"),
    probe_cluster_sep=(0.3, "Illumina cluster separation; under 0.3 the clusters overlap"),
    sample_call_rate=(0.97, "Illumina's and MoChA's recommended minimum sample call rate"),
    sample_lrr_sd_flag=(0.20, "autosomal LRR SD after correction above which mosaic calls are unreliable (MoChA's practice; PennCNV filters at 0.30)"),
    sample_lrr_sd_fail=(0.30, "PennCNV's filter_cnv --qclrrsd default: the sample is not to be read"),
    sample_baf_sd=(0.06, "BAF SD at heterozygous calls (MoChA's baf_auto flags contamination at 0.03; the band SD itself above 0.06 blurs the bands)"),
    sample_baf_drift=(0.002, "share of homozygous-called probes with BAF in (0.2, 0.25) or (0.75, 0.8); PennCNV's drift over every probe filters at 0.01"),
    sample_wf=(0.05, "|wave factor|: PennCNV's waviness factor, the median absolute deviation of the 1-Mb window LRR medians (raw, autosomes), signed by the GC correlation; PennCNV's --qcwf default 0.05"),
    sample_lrr_mean=(0.1, "autosomal raw LRR mean beyond +-0.1: a poorly normalised sample"),
    sample_event_z=(3.0, "event count more than 3 robust SD above the cohort median (after the events step)"),
    batch_min=(20, "passing samples a batch needs for its own reference; else the cohort-wide reference is used"),
    baf_recentre=(0.01, "heterozygous-call BAF re-centred per probe when the cohort median shifts by more than this (MoChA)"),
)


def window_gc(pos, gc, half=500_000):
    """The mean GC of the probes within +-half bp of each probe."""
    ok = np.isfinite(gc)
    cs = np.concatenate([[0.0], np.cumsum(np.where(ok, gc, 0.0))])
    cn = np.concatenate([[0], np.cumsum(ok)])
    a = np.searchsorted(pos, pos - half, side="left")
    b = np.searchsorted(pos, pos + half, side="right")
    n = cn[b] - cn[a]
    return np.where(n > 0, (cs[b] - cs[a]) / np.maximum(n, 1), np.nan)


def design(gc_z, w_z):
    """The GC-wave design matrix: 1, g, g^2, w, w^2 for standardised probe GC g and window GC w (NaN -> 0)."""
    g = np.nan_to_num(gc_z)
    w = np.nan_to_num(w_z)
    return np.stack([np.ones_like(g), g, g * g, w, w * w], axis=1)


def gc_covariates(cache, chrom, stats):
    """Standardised probe GC and 1-Mb window GC of one chromosome (the cohort's autosomal mean and SD)."""
    pr = cache.probes(chrom)
    gc = pr.gc.astype(np.float64)
    w = window_gc(pr.pos, gc)
    return (gc - stats["gc_mean"]) / stats["gc_sd"], (w - stats["w_mean"]) / stats["w_sd"]


def gc_stats(cache):
    """Mean and SD of the probe GC and of the 1-Mb window GC over the autosomes (the standardisation every fit uses)."""
    g, w = [], []
    for c in cache.chroms:
        if c in ("chrX", "chrY"):
            continue
        pr = cache.probes(c)
        gc = pr.gc.astype(np.float64)
        g.append(gc[np.isfinite(gc)])
        ww = window_gc(pr.pos, gc)
        w.append(ww[np.isfinite(ww)])
    g, w = np.concatenate(g), np.concatenate(w)
    return dict(gc_mean=float(g.mean()), gc_sd=float(g.std() + 1e-9), w_mean=float(w.mean()), w_sd=float(w.std() + 1e-9))


# ------------------------------------------------------------------------------------------------------------- pass 1

def pass1_chrom(args):
    """Raw per-probe and per-sample statistics of one chromosome, and the GC regression accumulators."""
    cache_dir, chrom, stats = args["cache"], args["chrom"], args["gc_stats"]
    C = Cache(cache_dir)
    pr = C.probes(chrom)
    S, P = C.n_samples, pr.n
    gt, lrr, baf = C.gt(chrom), C.lrr(chrom), C.baf(chrom)
    gc_z, w_z = gc_covariates(C, chrom, stats)
    X = design(gc_z, w_z)
    auto = chrom not in ("chrX", "chrY")
    p_end = C.genome.p_end.get(chrom, 0)
    arm_of = (pr.pos > p_end).astype(np.int8)                                   # 0 the short arm, 1 the long: the GC fit is on within-arm deviations
    acc = {k: dict(n=0, sx=np.zeros(4), sy=np.zeros(S), sxx=np.zeros((4, 4)), sxy=np.zeros((4, S))) for k in (0, 1)}
    out = dict(chrom=chrom, lrr_median=np.full(P, np.nan), lrr_sd=np.full(P, np.nan), call_rate=np.full(P, np.nan), het_baf_median=np.full(P, np.nan), n_called=np.zeros(P, np.int32),
               n_het=np.zeros(P, np.int32), n_alt=np.zeros(P, np.int32),
               s_called=np.zeros(S, np.int64), s_total=np.zeros(S, np.int64), s_het=np.zeros(S, np.int64), s_lrr_n=np.zeros(S, np.int64), s_lrr_sum=np.zeros(S), s_lrr_sq=np.zeros(S),
               s_hetbaf_n=np.zeros(S, np.int64), s_hetbaf_sq=np.zeros(S), s_hom_n=np.zeros(S, np.int64), s_drift=np.zeros(S, np.int64), s_hist=np.zeros((S, 500), np.int64),
               arms=[], s_wy=np.zeros(S), s_w2=0.0, s_wn=0,
               y_lrr_sum=np.zeros(S), y_lrr_n=np.zeros(S, np.int64), y_called=np.zeros(S, np.int64), y_total=np.zeros(S, np.int64),
               win_med=np.zeros((0, S)), win_gc=np.zeros(0), auto=auto)
    if auto:                                                                 # PennCNV's 1-Mb windows for the waviness factor
        edges = np.arange(0, int(pr.pos[-1]) + 1_000_000, 1_000_000)
        k = np.searchsorted(pr.pos, edges)
        meds, gcs = [], []
        for w in range(len(k) - 1):
            a, b = k[w], k[w + 1]
            if b - a > 10:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    meds.append(np.nanmedian(np.asarray(lrr[:, a:b], dtype=np.float64), axis=1))
                gcs.append(float(np.nanmean(pr.gc[a:b])))
        if meds:
            out["win_med"], out["win_gc"] = np.stack(meds), np.asarray(gcs)
    for a in range(0, P, BLOCK):
        b = min(a + BLOCK, P)
        g = np.asarray(gt[:, a:b])
        L = np.asarray(lrr[:, a:b], dtype=np.float64)
        B = np.asarray(baf[:, a:b], dtype=np.float64)
        called = g >= 0
        het = g == 1
        finiteL = np.isfinite(L)
        snv = pr.snv[a:b] & ~pr.intensity_only[a:b]
        par = pr.par[a:b]
        # per probe
        out["n_called"][a:b] = called.sum(axis=0)
        out["n_het"][a:b] = het.sum(axis=0)
        out["n_alt"][a:b] = np.where(called, g, 0).sum(axis=0)
        out["call_rate"][a:b] = called.sum(axis=0) / float(S)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out["lrr_median"][a:b] = np.nanmedian(L, axis=0)
            q = np.nanpercentile(L, [25, 75], axis=0)
            out["lrr_sd"][a:b] = (q[1] - q[0]) / 1.349
            hb = np.where(het, B, np.nan)
            out["het_baf_median"][a:b] = np.nanmedian(hb, axis=0)
        # per sample (autosomes for the metrics; the X and Y for the sex signals)
        if auto:
            use = snv[None, :]
            out["s_called"] += (called & use).sum(axis=1)
            out["s_total"] += np.broadcast_to(use, g.shape).sum(axis=1)
            out["s_het"] += (het & use).sum(axis=1)
            Lf = np.where(finiteL, L, 0.0)
            out["s_lrr_n"] += finiteL.sum(axis=1)
            out["s_lrr_sum"] += Lf.sum(axis=1)
            out["s_lrr_sq"] += (Lf * Lf).sum(axis=1)
            hu = het & use & np.isfinite(B)
            dev = np.where(hu, np.clip(B, 0, 1) - 0.5, 0.0)
            out["s_hetbaf_n"] += hu.sum(axis=1)
            out["s_hetbaf_sq"] += (dev * dev).sum(axis=1)                 # per chromosome: the median over chromosomes is the sample's BAF SD
            bins = np.minimum((np.abs(dev) * 1000).astype(np.int64), 499)   # the folded values histogrammed at 0.001: the robust band SD comes from their median
            for j in range(S):
                hj = bins[j][hu[j]]
                if hj.size:
                    out["s_hist"][j] += np.bincount(hj, minlength=500)
            hom = (g == 0) | (g == 2)
            homu = hom & use & np.isfinite(B)
            out["s_hom_n"] += homu.sum(axis=1)
            out["s_drift"] += (homu & (((B > 0.2) & (B < 0.25)) | ((B > 0.75) & (B < 0.8)))).sum(axis=1)
            # the GC regression on within-arm deviations: y = clipped LRR, NaN -> 0; the four slope terms of X, shared by every sample;
            # the arm sums let the within-arm products be formed afterwards (sum (x - xbar)(y - ybar) = sum xy - n xbar ybar)
            y = np.where(finiteL, np.clip(L, -1.0, 1.0), 0.0)
            Xb = X[a:b, 1:]
            for arm in (0, 1):
                m = arm_of[a:b] == arm
                if m.any():
                    A = acc[arm]
                    A["n"] += int(m.sum())
                    A["sx"] += Xb[m].sum(axis=0)
                    A["sy"] += y[:, m].sum(axis=1)
                    A["sxx"] += Xb[m].T @ Xb[m]
                    A["sxy"] += Xb[m].T @ y[:, m].T
            wz = np.nan_to_num(w_z[a:b])
            out["s_wy"] += (y * wz[None, :]).sum(axis=1)
            out["s_w2"] += float((wz * wz).sum())
            out["s_wn"] += int(wz.size)
        elif chrom == "chrY":
            use = (~par)[None, :]
            Lf = np.where(finiteL & use, L, 0.0)
            out["y_lrr_sum"] += Lf.sum(axis=1)
            out["y_lrr_n"] += (finiteL & use).sum(axis=1)
            out["y_called"] += (called & use & snv[None, :]).sum(axis=1)
            out["y_total"] += np.broadcast_to(use & snv[None, :], g.shape).sum(axis=1)
        if (a // BLOCK) % 8 == 7:
            for mm in (gt, lrr, baf):
                release(mm)
    if auto:
        for arm in (0, 1):
            A = acc[arm]
            if A["n"] > 10:
                xbar, ybar = A["sx"] / A["n"], A["sy"] / A["n"]
                out["arms"].append((A["sxx"] - A["n"] * np.outer(xbar, xbar), A["sxy"] - A["n"] * np.outer(xbar, ybar)))
    return out


def combine_pass1(parts, cache):
    S = cache.n_samples
    s = dict(called=np.zeros(S, np.int64), total=np.zeros(S, np.int64), het=np.zeros(S, np.int64), lrr_n=np.zeros(S, np.int64), lrr_sum=np.zeros(S), lrr_sq=np.zeros(S),
             hetbaf_n=np.zeros(S, np.int64), hetbaf_sq=np.zeros(S), hom_n=np.zeros(S, np.int64), drift=np.zeros(S, np.int64), hist=np.zeros((S, 500), np.int64), xtx=np.zeros((4, 4)), xty=np.zeros((4, S)),
             wy=np.zeros(S), w2=0.0, wn=0, y_lrr_sum=np.zeros(S), y_lrr_n=np.zeros(S, np.int64), y_called=np.zeros(S, np.int64), y_total=np.zeros(S, np.int64))
    per_probe = {}
    baf_sd_chrom, win_med, win_gc = [], [], []
    for p in parts:
        per_probe[p["chrom"]] = {k: p[k] for k in ("lrr_median", "lrr_sd", "call_rate", "het_baf_median", "n_called", "n_het", "n_alt")}
        for k in s:
            if k in ("xtx", "xty"):
                continue
            s[k] = s[k] + p["s_" + k] if ("s_" + k) in p else s[k] + p[k]
        for xtx_arm, xty_arm in p["arms"]:
            s["xtx"] += xtx_arm
            s["xty"] += xty_arm
        if p["auto"] and p["s_hetbaf_n"].max() > 0:
            baf_sd_chrom.append(np.where(p["s_hetbaf_n"] >= 20, np.sqrt(p["s_hetbaf_sq"] / np.maximum(p["s_hetbaf_n"], 1)), np.nan))
        if p["win_med"].shape[0]:
            win_med.append(p["win_med"]), win_gc.append(p["win_gc"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        s["baf_sd"] = np.nanmedian(np.stack(baf_sd_chrom), axis=0) if baf_sd_chrom else np.full(S, np.nan)
        # the robust band SD: the median of |BAF - 0.5| at heterozygous calls over 0.6745 (a half-normal's median), immune to mis-called genotypes
        cum = np.cumsum(s["hist"], axis=1)
        tot = cum[:, -1]
        med_bin = np.argmax(cum >= (tot / 2.0)[:, None], axis=1)
        s["baf_sd_robust"] = np.where(tot > 50, (med_bin + 0.5) / 1000.0 / 0.6745, np.nan)
        if win_med:
            W, G = np.concatenate(win_med), np.concatenate(win_gc)               # windows x S, windows
            yi = W - np.nanmedian(W, axis=0)[None, :]
            wf = np.nanmedian(np.abs(yi), axis=0)
            Gc = G - np.nanmean(G)
            num = np.nansum(np.nan_to_num(yi) * Gc[:, None], axis=0)
            s["wf"] = np.where(num >= 0, wf, -wf)
        else:
            s["wf"] = np.full(S, np.nan)
    return s, per_probe


# ------------------------------------------------------------------------------------------------------------- pass 2

def pass2_chrom(args):
    """With the GC fit known: the reference per probe (cohort-wide and per batch, over the passing samples; on the X over the
    samples without a Y), and per-sample corrected-LRR accumulators and the GC correlation after correction."""
    cache_dir, chrom, stats, coef, passing, batch_idx, n_batches, ref_samples_x = (args[k] for k in ("cache", "chrom", "gc_stats", "coef", "passing", "batch_idx", "n_batches", "ref_samples_x"))
    C = Cache(cache_dir)
    pr = C.probes(chrom)
    S, P = C.n_samples, pr.n
    lrr = C.lrr(chrom)
    gc_z, w_z = gc_covariates(C, chrom, stats)
    X = design(gc_z, w_z)
    auto = chrom not in ("chrX", "chrY")
    ref_mask = passing if auto else (ref_samples_x if chrom == "chrX" else np.zeros(S, bool))
    out = dict(chrom=chrom, ref=np.full(P, np.nan), ref_batch=np.full((n_batches, P), np.nan), s_n=np.zeros(S, np.int64), s_sum=np.zeros(S), s_sq=np.zeros(S), s_wy=np.zeros(S), s_w2=0.0,
               s_absdiff=np.zeros(S), s_absn=np.zeros(S, np.int64), s_sample=[])
    for a in range(0, P, BLOCK):
        b = min(a + BLOCK, P)
        L = np.asarray(lrr[:, a:b], dtype=np.float64)
        fit = (X[a:b] @ coef).T                                           # S x block
        Lc = L - fit
        if chrom == "chrY":
            continue
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            if ref_mask.any():
                out["ref"][a:b] = np.nanmedian(Lc[ref_mask], axis=0)
            for k in range(n_batches):
                m = ref_mask & (batch_idx == k)
                if m.any():
                    out["ref_batch"][k, a:b] = np.nanmedian(Lc[m], axis=0)
        if auto:
            r = out["ref"][a:b]
            Lr = Lc - np.where(np.isfinite(r), r, 0.0)[None, :]
            ok = np.isfinite(Lr)
            v = np.where(ok, Lr, 0.0)
            out["s_n"] += ok.sum(axis=1)
            out["s_sum"] += v.sum(axis=1)
            out["s_sq"] += (v * v).sum(axis=1)
            wz = np.nan_to_num(w_z[a:b])
            out["s_wy"] += (v * wz[None, :]).sum(axis=1)
            out["s_w2"] += float((wz * wz).sum())
            d = np.abs(np.diff(v, axis=1))
            out["s_absdiff"] += d.sum(axis=1)
            out["s_absn"] += d.shape[1]
            out["s_sample"].append(np.where(ok, Lr, np.nan)[:, ::97].astype(np.float32))      # every 97th probe: the per-sample median offset comes from these
        if (a // BLOCK) % 8 == 7:
            release(lrr)
    out["s_sample"] = np.concatenate(out["s_sample"], axis=1) if out["s_sample"] else np.zeros((S, 0), np.float32)
    return out


class Correction:
    """The GC fit per sample and the reference per probe, applied to any block of the cache."""

    def __init__(self, out):
        z = np.load(os.path.join(out, "qc", "correction.npz"), allow_pickle=True)
        self.coef = z["coef"]                        # 5 x S
        self.batch_idx = z["batch_idx"]              # S (-1: cohort-wide)
        self.gc_stats = {k: float(z["gc_" + k]) for k in ("mean", "sd")}
        self.gc_stats = dict(gc_mean=float(z["gc_mean"]), gc_sd=float(z["gc_sd"]), w_mean=float(z["w_mean"]), w_sd=float(z["w_sd"]))
        self.ref = {str(c): z["ref_" + str(c)] for c in z["chroms"]}
        self.ref_batch = {str(c): z["refb_" + str(c)] for c in z["chroms"]}
        self.use_batch = bool(z["use_batch"])
        self.offset = z["offset"] if "offset" in z.files else np.zeros(self.coef.shape[1])
        self._cov = {}

    def covariates(self, cache, chrom):
        if chrom not in self._cov:
            self._cov[chrom] = gc_covariates(cache, chrom, self.gc_stats)
        return self._cov[chrom]

    def apply(self, cache, chrom, sample_idx=None, probe_slice=None, reference=True):
        """Corrected LRR (float32, samples x probes) of a chromosome for the given samples (indices; None = all) and probe slice."""
        lrr = cache.lrr(chrom)
        sl = probe_slice if probe_slice is not None else slice(0, lrr.shape[1])
        idx = np.arange(cache.n_samples) if sample_idx is None else np.atleast_1d(sample_idx)
        L = np.asarray(lrr[idx, sl], dtype=np.float32)
        gc_z, w_z = self.covariates(cache, chrom)
        X = design(gc_z[sl], w_z[sl])
        L -= (X @ self.coef[:, idx]).T.astype(np.float32)
        if reference and chrom != "chrY" and chrom in self.ref:
            r = self.ref[chrom][sl]
            if self.use_batch:
                rb = self.ref_batch[chrom]
                for j, i in enumerate(idx):
                    k = self.batch_idx[i]
                    rr = rb[k, sl] if k >= 0 else r
                    rr = np.where(np.isfinite(rr), rr, np.where(np.isfinite(r), r, 0.0))
                    L[j] -= rr.astype(np.float32)
            else:
                L -= np.where(np.isfinite(r), r, 0.0).astype(np.float32)[None, :]
            L -= self.offset[idx].astype(np.float32)[:, None]
        return L


# ------------------------------------------------------------------------------------------------------------- the step

def hwe_p(n_aa, n_ab, n_bb):
    """Exact Hardy-Weinberg test p (Wigginton 2005), vectorised over probes by a loop over the het count distribution."""
    n_aa, n_ab, n_bb = (np.asarray(x, dtype=np.int64) for x in (n_aa, n_ab, n_bb))
    n = n_aa + n_ab + n_bb
    p = np.ones(n.shape, dtype=np.float64)
    for i in range(n.size):
        p[i] = _hwe_one(int(n_aa[i]), int(n_ab[i]), int(n_bb[i]))
    return p


def _hwe_one(aa, ab, bb):
    n = aa + ab + bb
    if n == 0:
        return 1.0
    rare = 2 * min(aa, bb) + ab
    if rare == 0:
        return 1.0
    # the distribution of heterozygote counts given the allele counts
    mid = rare * (2 * n - rare) // (2 * n)
    if (mid % 2) != (rare % 2):
        mid += 1
    probs = {mid: 1.0}
    cur_hets, cur_homr = mid, (rare - mid) // 2
    cur_homc = n - cur_hets - cur_homr
    while cur_hets > 1:
        probs[cur_hets - 2] = probs[cur_hets] * cur_hets * (cur_hets - 1.0) / (4.0 * (cur_homr + 1.0) * (cur_homc + 1.0))
        cur_hets -= 2
        cur_homr += 1
        cur_homc += 1
    cur_hets, cur_homr = mid, (rare - mid) // 2
    cur_homc = n - cur_hets - cur_homr
    while cur_hets <= rare - 2:
        probs[cur_hets + 2] = probs[cur_hets] * 4.0 * cur_homr * cur_homc / ((cur_hets + 2.0) * (cur_hets + 1.0))
        cur_hets += 2
        cur_homr -= 1
        cur_homc -= 1
    tot = sum(probs.values())
    target = probs.get(ab, 0.0)
    return float(min(1.0, sum(v for v in probs.values() if v <= target * (1 + 1e-9)) / tot))


def run_qc(out, cache_dir=None, sample_sheet="", exclude_probes="", exclude_samples="", threads=1, params=None, log_fn=log):
    t0 = time.time()
    cache_dir = cache_dir or os.path.join(out, "cache")
    C = Cache(cache_dir)
    S = C.n_samples
    T = {k: v[0] for k, v in THRESHOLDS.items()}
    if params:
        T.update(params)
    qdir = os.path.join(out, "qc")
    os.makedirs(qdir, exist_ok=True)
    sheet, sheet_cols = (read_sample_sheet(sample_sheet) if sample_sheet else ({}, []))
    reported = np.array([sheet.get(s, {}).get("reported_sex", "") for s in C.samples])
    batch = np.array([sheet.get(s, {}).get("batch", "") for s in C.samples])
    excl_s = set(_read_list(exclude_samples))
    excl_p = set(_read_list(exclude_probes))
    stats = gc_stats(C)
    log_fn("%d samples, %d probes on %d chromosomes; GC mean %.3f (SD %.3f), 1-Mb GC SD %.3f" % (S, sum(C.n_probes.values()), len(C.chroms), stats["gc_mean"], stats["gc_sd"], stats["w_sd"]), "qc")
    # pass 1
    jobs = [dict(cache=cache_dir, chrom=c, gc_stats=stats) for c in C.chroms]
    parts = _pool(pass1_chrom, jobs, threads)
    s, per_probe = combine_pass1(parts, C)
    # the GC fit per sample: the four slopes on within-arm deviations (no intercept: the per-sample median offset centres the result)
    xtx = s["xtx"] + 1e-9 * np.eye(4)
    coef = np.vstack([np.zeros((1, S)), np.linalg.solve(xtx, s["xty"])])      # 5 x S, the intercept row 0
    wf = s["wf"]
    # the per-sample raw metrics
    call_rate = s["called"] / np.maximum(s["total"], 1)
    het_rate = s["het"] / np.maximum(s["called"], 1)
    lrr_mean = s["lrr_sum"] / np.maximum(s["lrr_n"], 1)
    lrr_sd_raw = np.sqrt(np.maximum(s["lrr_sq"] / np.maximum(s["lrr_n"], 1) - lrr_mean ** 2, 0))
    baf_sd = s["baf_sd"]
    baf_sd_robust = s["baf_sd_robust"]
    drift = s["drift"] / np.maximum(s["hom_n"], 1)
    y_lrr = np.where(s["y_lrr_n"] > 0, s["y_lrr_sum"] / np.maximum(s["y_lrr_n"], 1), np.nan)
    y_call = s["y_called"] / np.maximum(s["y_total"], 1)
    gc_corr_before = (s["wy"] / max(s["wn"], 1)) / np.maximum(lrr_sd_raw * np.sqrt(s["w2"] / max(s["wn"], 1)), 1e-9)
    # provisional passing (the reference is built from these): call rate, hard LRR SD, not excluded
    prov_pass = (call_rate >= T["sample_call_rate"]) & (lrr_sd_raw < T["sample_lrr_sd_fail"] + 0.1) & ~np.isin(C.samples, list(excl_s))
    if prov_pass.sum() < max(3, S // 10):
        log_fn("only %d of %d samples pass provisionally; the reference uses every sample" % (prov_pass.sum(), S), "qc")
        prov_pass = ~np.isin(C.samples, list(excl_s))
    # provisional sex from the Y (the X reference comes from the samples without a Y)
    y_mid = _y_midpoint(y_lrr)
    no_y = np.isfinite(y_lrr) & (y_lrr < y_mid) if np.isfinite(y_mid) else np.zeros(S, bool)
    ref_x = prov_pass & no_y
    if ref_x.sum() < 10:
        ref_x = prov_pass.copy()
        log_fn("fewer than 10 passing samples without a Y: the X reference is cohort-wide (one X will not read as one copy)", "qc")
    # batches
    batches = sorted({b for b in batch if b})
    batch_idx = np.array([batches.index(b) if b else -1 for b in batch])
    use_batch = False
    if batches:
        counts = {b: int(((batch_idx == i) & prov_pass).sum()) for i, b in enumerate(batches)}
        use_batch = all(v >= T["batch_min"] for v in counts.values())
        log_fn("batches: %s -> %s" % (", ".join("%s %d passing" % kv for kv in counts.items()), "per-batch reference" if use_batch else "cohort-wide reference (a batch has fewer than %d passing samples)" % T["batch_min"]), "qc")
    if not use_batch:
        batch_idx = np.full(S, -1)
    # pass 2
    jobs = [dict(cache=cache_dir, chrom=c, gc_stats=stats, coef=coef, passing=prov_pass, batch_idx=batch_idx, n_batches=len(batches) if use_batch else 0, ref_samples_x=ref_x) for c in C.chroms]
    parts2 = _pool(pass2_chrom, jobs, threads)
    n2 = sum(p["s_n"] for p in parts2)
    sum2 = sum(p["s_sum"] for p in parts2)
    sq2 = sum(p["s_sq"] for p in parts2)
    wy2 = sum(p["s_wy"] for p in parts2)
    w22 = sum(p["s_w2"] for p in parts2)
    absd = sum(p["s_absdiff"] for p in parts2) / np.maximum(sum(p["s_absn"] for p in parts2), 1)
    lrr_mean_c = sum2 / np.maximum(n2, 1)
    lrr_sd = np.sqrt(np.maximum(sq2 / np.maximum(n2, 1) - lrr_mean_c ** 2, 0))
    # the per-sample offset: the median corrected autosomal LRR, subtracted so that two copies sit at 0 in every sample whatever share
    # of its genome an event covers (the regression's intercept is a mean, which a large event pulls)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sub = np.concatenate([p["s_sample"] for p in parts2 if p["s_sample"].shape[1]], axis=1) if any(p["s_sample"].shape[1] for p in parts2) else np.zeros((S, 0))
        offset = np.nanmedian(sub, axis=1) if sub.shape[1] else np.zeros(S)
    offset = np.where(np.isfinite(offset), offset, 0.0)
    gc_corr_after = (wy2 / max(s["wn"], 1)) / np.maximum(lrr_sd * np.sqrt(w22 / max(s["wn"], 1)), 1e-9)
    lrr_sd_diff = absd / np.sqrt(2.0) * np.sqrt(np.pi / 2.0)             # the noise SD from first differences (robust to events)
    np.savez(os.path.join(qdir, "correction.npz"), coef=coef, batch_idx=batch_idx, use_batch=use_batch, gc_mean=stats["gc_mean"], gc_sd=stats["gc_sd"], w_mean=stats["w_mean"], w_sd=stats["w_sd"],
             offset=offset, chroms=np.array(C.chroms), batches=np.array(batches), **{"ref_" + p["chrom"]: p["ref"] for p in parts2}, **{"refb_" + p["chrom"]: p["ref_batch"] for p in parts2})
    # probe QC
    probe_rows, keep_all, keep_gt_all, maf_all, cr_all, hwe_all = [], [], [], [], [], []
    noisy_cut = np.nanpercentile(np.concatenate([per_probe[c]["lrr_sd"] for c in C.chroms if c not in ("chrX", "chrY")]), 100 * (1 - T["probe_lrr_sd_top"]))
    seen = set()
    for c in C.chroms:
        pr = C.probes(c)
        pp = per_probe[c]
        P = pr.n
        keep = np.ones(P, bool)
        keep_gt = pr.snv & ~pr.intensity_only
        reasons = [[] for _ in range(P)]
        ids = pr.ids
        for i in range(P):
            if pr.intensity_only[i]:
                reasons[i].append("intensity_only")
            elif not pr.snv[i]:
                reasons[i].append("non_snv")
            key = (c, int(pr.pos[i]))
            if key in seen:
                keep[i] = False
                keep_gt[i] = False
                reasons[i].append("duplicate_position")
            seen.add(key)
            if ids[i] in excl_p:
                keep[i] = False
                keep_gt[i] = False
                reasons[i].append("excluded")
        med = pp["lrr_median"]
        bad = np.isfinite(med) & (np.abs(med) > T["probe_lrr_median"]) & (c not in ("chrX", "chrY"))
        noisy = np.isfinite(pp["lrr_sd"]) & (pp["lrr_sd"] > noisy_cut) & (c not in ("chrX", "chrY"))
        if c == "chrY" and (~no_y & prov_pass).sum() >= 5:
            cr = np.where(pr.par, pp["call_rate"], s_y_call_rate(C, c, prov_pass & ~no_y))      # the Y itself among the samples with a Y; its PAR probes among everyone
        else:
            cr = pp["call_rate"]
        lowcr = keep_gt & (cr < T["probe_call_rate"])
        shift = keep_gt & np.isfinite(pp["het_baf_median"]) & (np.abs(pp["het_baf_median"] - 0.5) > T["probe_baf_shift"]) & (pp["n_het"] >= 10)
        gt_low = keep_gt & np.isfinite(pr.gentrain) & (pr.gentrain < T["probe_gentrain"])
        cs_low = keep_gt & np.isfinite(pr.cluster_sep) & (pr.cluster_sep < T["probe_cluster_sep"])
        n_called = pp["n_called"]
        n_ab = pp["n_het"]
        n_bb = (pp["n_alt"] - n_ab) // 2
        n_aa = n_called - n_ab - n_bb
        hwe = np.ones(P)
        if c not in ("chrX", "chrY"):
            cand = keep_gt & (n_called >= 10) & (n_ab > 0)
            hwe[cand] = hwe_p(n_aa[cand], n_ab[cand], n_bb[cand])
        hwe_fail = keep_gt & (hwe < T["probe_hwe_p"])
        for i in range(P):
            for flag, name in ((bad[i], "lrr_median"), (noisy[i], "lrr_sd"), (lowcr[i], "call_rate"), (shift[i], "baf_shift"), (gt_low[i], "gentrain"), (cs_low[i], "cluster_sep"), (hwe_fail[i], "hwe")):
                if flag:
                    reasons[i].append(name)
        keep &= ~(bad | noisy)
        keep_gt &= ~(lowcr | shift | gt_low | cs_low | hwe_fail | bad | noisy)
        maf = np.where(n_called > 0, pp["n_alt"] / np.maximum(2 * n_called, 1), np.nan)
        maf = np.minimum(maf, 1 - maf)
        for i in range(P):
            probe_rows.append((c, int(pr.pos[i]), ids[i], int(keep[i]), int(keep_gt[i]), ",".join(reasons[i]) or NA, int(pr.par[i]), int(pr.intensity_only[i]), int(pr.snv[i]),
                               cr[i], med[i], pp["lrr_sd"][i], pp["het_baf_median"][i], hwe[i] if c not in ("chrX", "chrY") else NA, pr.gentrain[i], pr.cluster_sep[i], maf[i], int(n_called[i])))
        keep_all.append(keep), keep_gt_all.append(keep_gt), maf_all.append(maf), cr_all.append(cr), hwe_all.append(hwe)
    write_tsv(os.path.join(qdir, "probes.tsv"), ["chrom", "pos", "id", "keep", "keep_gt", "reason", "par", "intensity_only", "snv", "call_rate", "lrr_median", "lrr_sd", "het_baf_median",
                                                 "hwe_p", "gentrain", "cluster_sep", "maf", "n_called"], probe_rows, nd=4)
    het_shift = {c: per_probe[c]["het_baf_median"] for c in C.chroms}
    np.savez(os.path.join(qdir, "probes.npz"), keep=np.concatenate(keep_all), keep_gt=np.concatenate(keep_gt_all), maf=np.concatenate(maf_all), call_rate=np.concatenate(cr_all),
             hwe=np.concatenate(hwe_all), het_baf_median=np.concatenate([het_shift[c] for c in C.chroms]), n_het=np.concatenate([per_probe[c]["n_het"] for c in C.chroms]), chroms=np.array(C.chroms))
    # sample QC verdicts
    rows = []
    cohort = {}
    flags_all = []
    for i, name in enumerate(C.samples):
        fl, hard = [], False
        if name in excl_s:
            fl.append("excluded")
            hard = True
        if call_rate[i] < T["sample_call_rate"]:
            fl.append("call_rate")
        if lrr_sd[i] > T["sample_lrr_sd_fail"]:
            fl.append("lrr_sd_fail")
            hard = True
        elif lrr_sd[i] > T["sample_lrr_sd_flag"]:
            fl.append("lrr_sd")
        if baf_sd[i] > T["sample_baf_sd"]:
            fl.append("baf_sd")
        if drift[i] > T["sample_baf_drift"]:
            fl.append("baf_drift")
        if abs(wf[i]) > T["sample_wf"]:
            fl.append("wave")
        if abs(lrr_mean[i]) > T["sample_lrr_mean"]:
            fl.append("lrr_mean")
        flags_all.append(fl)
        rows.append(dict(sample=name, reported_sex=reported[i] or NA, batch=batch[i] or NA, call_rate=call_rate[i], het_rate=het_rate[i], lrr_mean=lrr_mean[i], lrr_offset=offset[i], lrr_sd_raw=lrr_sd_raw[i],
                         lrr_sd=lrr_sd[i], lrr_sd_diff=lrr_sd_diff[i], baf_sd=baf_sd[i], baf_sd_robust=baf_sd_robust[i], baf_drift=drift[i], wf=wf[i], gc_corr_before=gc_corr_before[i], gc_corr_after=gc_corr_after[i],
                         y_lrr_raw=y_lrr[i], y_call_rate=y_call[i], **{"pass": int(not fl), "hard_fail": int(hard), "flags": ",".join(fl) or NA}))
    for k in ("call_rate", "het_rate", "lrr_mean", "lrr_sd_raw", "lrr_sd", "baf_sd", "baf_drift", "wf", "gc_corr_before", "gc_corr_after"):
        v = np.array([r[k] for r in rows], dtype=float)
        cohort[k] = dict(median=float(np.nanmedian(v)), robust_sd=robust_sd(v))
    cols = ["sample", "reported_sex", "batch", "call_rate", "het_rate", "lrr_mean", "lrr_offset", "lrr_sd_raw", "lrr_sd", "lrr_sd_diff", "baf_sd", "baf_sd_robust", "baf_drift", "wf", "gc_corr_before", "gc_corr_after",
            "y_lrr_raw", "y_call_rate", "pass", "hard_fail", "flags"]
    write_tsv(os.path.join(qdir, "samples.tsv"), cols, rows, nd=4)
    n_pass = sum(1 for r in rows if r["pass"])
    n_hard = sum(1 for r in rows if r["hard_fail"])
    summary = dict(samples=S, passing=n_pass, flagged=S - n_pass, hard_fail=n_hard, probes=len(probe_rows), probes_keep=int(np.concatenate(keep_all).sum()), probes_keep_gt=int(np.concatenate(keep_gt_all).sum()),
                   reasons=_count([r[5] for r in probe_rows]), flags=_count([",".join(f) or NA for f in flags_all]), cohort=cohort, thresholds={k: dict(value=T[k], source=THRESHOLDS[k][1]) for k in THRESHOLDS},
                   batches=batches, use_batch=use_batch, gc=stats, y_midpoint=float(y_mid) if np.isfinite(y_mid) else None, x_reference_samples=int(ref_x.sum()), seconds=time.time() - t0)
    write_json(summary, os.path.join(qdir, "summary.json"))
    log_fn("probes: %d kept for LRR, %d for genotypes of %d; samples: %d pass, %d flagged (%d hard) of %d; LRR SD median %.3f raw -> %.3f corrected; GC correlation %.3f -> %.3f" % (
        summary["probes_keep"], summary["probes_keep_gt"], len(probe_rows), n_pass, S - n_pass, n_hard, S, cohort["lrr_sd_raw"]["median"], cohort["lrr_sd"]["median"],
        float(np.nanmedian(np.abs(gc_corr_before))), float(np.nanmedian(np.abs(gc_corr_after)))), "qc")
    try:
        from .plots import fig_qc
        fig_qc(out, C, rows, T, Correction(out), log_fn)
    except Exception as e:                                                 # figures never block the tables
        log_fn("figures skipped: %s" % e, "qc")
    mark_done(out, "qc", dict(seconds=time.time() - t0))
    return summary


def s_y_call_rate(C, chrom, males):
    gt = C.gt(chrom)
    g = np.asarray(gt[males])
    return (g >= 0).sum(axis=0) / float(max(males.sum(), 1))


def _y_midpoint(y_lrr):
    """The midpoint between the two modes of the raw Y level (present near 0, absent far below); NaN when one mode."""
    v = y_lrr[np.isfinite(y_lrr)]
    if v.size < 4:
        return np.nan
    lo, hi = v.min(), v.max()
    if hi - lo < 1.0:
        return np.nan
    c1, c2 = np.percentile(v, 10), np.percentile(v, 90)
    for _ in range(20):
        mid = (c1 + c2) / 2
        a, b = v[v < mid], v[v >= mid]
        if a.size == 0 or b.size == 0:
            return np.nan
        c1, c2 = np.median(a), np.median(b)
    return (c1 + c2) / 2 if c2 - c1 > 1.0 else np.nan


def _count(items):
    out = {}
    for x in items:
        for k in x.split(","):
            out[k] = out.get(k, 0) + 1
    return out


def _read_list(path):
    if not path:
        return []
    with open(path) as fh:
        return [l.strip().split()[0] for l in fh if l.strip()]


def _pool(fn, jobs, threads):
    if threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            return list(ex.map(fn, jobs))
    return [fn(j) for j in jobs]


def read_sample_qc(out):
    """qc/samples.tsv as {sample: row} (strings)."""
    from .util import read_tsv
    return {r["sample"]: r for r in read_tsv(os.path.join(out, "qc", "samples.tsv"))}


def update_sample_qc(out, updates, new_cols=()):
    """Add or change columns of qc/samples.tsv: updates = {sample: {col: value}}; a flag added to `flags` and `pass` recomputed."""
    from .util import read_tsv
    path = os.path.join(out, "qc", "samples.tsv")
    rows = read_tsv(path)
    cols = list(rows[0].keys()) if rows else []
    for c in new_cols:
        if c not in cols:
            cols.insert(cols.index("pass"), c)
    for r in rows:
        u = updates.get(r["sample"], {})
        for k, v in u.items():
            if k == "add_flag":
                fl = [f for f in r["flags"].split(",") if f and f != NA]
                if v and v not in fl:
                    fl.append(v)
                r["flags"] = ",".join(fl) or NA
                r["pass"] = "0" if fl else "1"
            else:
                r[k] = v
        for c in new_cols:
            r.setdefault(c, NA)
    write_tsv(path, cols, rows)


def probe_masks(out, cache):
    """keep, keep_gt, maf, het_baf_median per chromosome from qc/probes.npz."""
    z = np.load(os.path.join(out, "qc", "probes.npz"))
    off, n = {}, 0
    for c in cache.chroms:
        off[c] = slice(n, n + cache.n_probes[c])
        n += cache.n_probes[c]
    n_het = z["n_het"] if "n_het" in z.files else np.full(z["keep"].size, 10 ** 6)
    return {c: dict(keep=z["keep"][off[c]], keep_gt=z["keep_gt"][off[c]], maf=z["maf"][off[c]], het_baf_median=z["het_baf_median"][off[c]], call_rate=z["call_rate"][off[c]], n_het=n_het[off[c]],
                    baf_shift=baf_shift(z["het_baf_median"][off[c]], n_het[off[c]])) for c in cache.chroms}


RECENTRE_MIN_HET = 30


def baf_shift(het_median, n_het, min_shift=0.01, min_het=RECENTRE_MIN_HET):
    """The per-probe re-centring of heterozygous-call BAF (MoChA's cluster adjustment): the cohort's median at the probe minus 0.5
    when it shifts by more than 0.01 and at least 30 heterozygous calls support it (the median of fewer calls is noise)."""
    return np.where(np.isfinite(het_median) & (np.abs(het_median - 0.5) > min_shift) & (n_het >= min_het), het_median - 0.5, 0.0).astype(np.float32)


def recompute_hwe(out, cache, unrelated):
    """Hardy-Weinberg over the unrelated passing samples (after the pedigree step): qc/probes.tsv and probes.npz updated (the
    `hwe` reason and keep_gt). Returns the number of probes failing."""
    from .util import read_tsv
    T = THRESHOLDS["probe_hwe_p"][0]
    path = os.path.join(out, "qc", "probes.tsv")
    rows = read_tsv(path)
    z = np.load(os.path.join(out, "qc", "probes.npz"))
    keep_gt, hwe_all = z["keep_gt"].copy(), z["hwe"].copy()
    off, n = {}, 0
    for c in cache.chroms:
        off[c] = n
        n += cache.n_probes[c]
    n_fail = 0
    for c in cache.chroms:
        if c in ("chrX", "chrY"):
            continue
        g = np.asarray(cache.gt(c)[unrelated])
        n_ab = (g == 1).sum(axis=0)
        n_bb = (g == 2).sum(axis=0)
        n_aa = (g == 0).sum(axis=0)
        n_called = n_aa + n_ab + n_bb
        sl = slice(off[c], off[c] + cache.n_probes[c])
        cand = (n_called >= 10) & (n_ab > 0)
        p = np.ones(cache.n_probes[c])
        p[cand] = hwe_p(n_aa[cand], n_ab[cand], n_bb[cand])
        hwe_all[sl] = p
        fail = p < T
        for i in np.flatnonzero(fail | (np.asarray([("hwe" in r["reason"]) for r in rows[sl]]))):
            r = rows[off[c] + i]
            reasons = [x for x in r["reason"].split(",") if x not in ("NA", "hwe", "")]
            if fail[i]:
                reasons.append("hwe")
            r["reason"] = ",".join(reasons) or NA
            r["hwe_p"] = fmt_p(p[i])
            was_only_hwe = (not fail[i]) and r["keep_gt"] == "0" and reasons == []
            r["keep_gt"] = "0" if (fail[i] or (reasons and r["keep_gt"] == "0")) else ("1" if was_only_hwe or r["keep_gt"] == "1" else r["keep_gt"])
            keep_gt[off[c] + i] = r["keep_gt"] == "1"
        n_fail += int(fail.sum())
    cols = list(rows[0].keys())
    write_tsv(path, cols, rows)
    np.savez(os.path.join(out, "qc", "probes.npz"), keep=z["keep"], keep_gt=keep_gt, maf=z["maf"], call_rate=z["call_rate"], hwe=hwe_all, het_baf_median=z["het_baf_median"], n_het=z["n_het"], chroms=z["chroms"])
    return n_fail


def fmt_p(p):
    return "%.3g" % p
