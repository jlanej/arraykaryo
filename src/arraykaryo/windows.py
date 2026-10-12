"""Window tracks per sample and chromosome: windows of about W probes pooled by probe count (a similar amount of evidence
everywhere); per window the median corrected LRR, the band deviation d at heterozygous calls (EM fit of the symmetric
two-band mixture with the sample's own band SD; baf.d_em_windows) with its log-likelihood ratio against no split, and the
heterozygous rate relative to the Hardy-Weinberg expectation from the cohort's allele frequencies at those probes. Computed for a block of samples at once and stored per
chromosome (events/windows/<chrom>.npz) for the calls, the figures and the report."""
import os

import numpy as np

from .baf import d_em_windows, fold

W_DEFAULT = 100


def window_layout(pos, keep, w=W_DEFAULT):
    """Windows over the kept probes of a chromosome: (bounds, mids) with bounds an array of probe-index boundaries [b0, b1, ...]
    (into the full probe array; the last is n_probes) such that each window holds about w kept probes, and mids their median
    positions. The last window is merged into the previous one when it is small."""
    kept = np.flatnonzero(keep)
    if kept.size == 0:
        return np.array([0, pos.size]), np.array([float(pos[0]) if pos.size else 0.0]), np.zeros(0, dtype=int)
    n_win = max(1, int(round(kept.size / float(w))))
    edges = np.linspace(0, kept.size, n_win + 1).astype(int)
    starts = kept[edges[:-1]]
    bounds = np.concatenate([starts, [pos.size]])
    bounds[0] = 0
    mids = np.array([np.median(pos[kept[edges[k]:edges[k + 1]]]) for k in range(n_win)])
    spans = np.array([(pos[kept[edges[k]]], pos[kept[edges[k + 1] - 1]]) for k in range(n_win)])
    return bounds, mids, spans


def gather_matrix(bounds):
    """An index matrix (n_win x width) into the probe array for every window, -1 padded."""
    widths = np.diff(bounds)
    width = int(widths.max())
    idx = np.full((bounds.size - 1, width), -1, dtype=np.int64)
    for k in range(bounds.size - 1):
        idx[k, :widths[k]] = np.arange(bounds[k], bounds[k + 1])
    return idx


def tracks_for_block(lrr, baf, gt, keep, keep_gt, maf, baf_shift, bounds, sigma):
    """For a block of samples (arrays B x P): lrr_med (B x n_win), d (B x n_win), n_het, het_ratio, n_called, n_lrr.
    sigma: the band SD per sample (B). baf_shift: per-probe re-centring (subtracted from BAF; 0 where none)."""
    B, P = lrr.shape
    idx = gather_matrix(bounds)
    n_win, width = idx.shape
    valid = idx >= 0
    safe = np.where(valid, idx, 0)
    # LRR: median over kept probes
    L = np.where(keep[None, :] & np.isfinite(lrr), lrr, np.nan)
    Lw = L[:, safe]                                     # B x n_win x width
    Lw = np.where(valid[None, :, :], Lw, np.nan)
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            lrr_med = np.nanmedian(Lw, axis=2)
    n_lrr = np.isfinite(Lw).sum(axis=2)
    # BAF: d at heterozygous calls
    het = (gt == 1) & keep_gt[None, :]
    u = fold(baf.astype(np.float32) - baf_shift[None, :].astype(np.float32))
    u = np.where(np.isfinite(u), u, 0.0)
    uw = u[:, safe]
    hw = het[:, safe] & valid[None, :, :] & np.isfinite(baf[:, safe])
    d, n_het, llr = d_em_windows(uw.reshape(B * n_win, width), hw.reshape(B * n_win, width), np.repeat(sigma, n_win))
    d, n_het, llr = d.reshape(B, n_win), n_het.reshape(B, n_win), llr.reshape(B, n_win)
    # heterozygous rate against the Hardy-Weinberg expectation at the called probes
    called = (gt >= 0) & keep_gt[None, :] & np.isfinite(maf)[None, :]
    exp_het = 2 * maf * (1 - maf)
    cw = called[:, safe] & valid[None, :, :]
    expected = (np.where(cw, np.nan_to_num(exp_het)[safe][None, :, :], 0.0)).sum(axis=2)
    observed = (het[:, safe] & cw).sum(axis=2)
    with np.errstate(all="ignore"):
        het_ratio = np.where(expected > 2, observed / expected, np.nan)
    return dict(lrr_med=lrr_med.astype(np.float32), d=d.astype(np.float32), d_llr=llr.astype(np.float32), n_het=n_het.astype(np.int16), het_ratio=het_ratio.astype(np.float32),
                n_called=cw.sum(axis=2).astype(np.int16), n_lrr=n_lrr.astype(np.int16))


def save_tracks(out, chrom, bounds, mids, spans, parts):
    os.makedirs(os.path.join(out, "events", "windows"), exist_ok=True)
    np.savez(os.path.join(out, "events", "windows", chrom + ".npz"), bounds=bounds, mids=mids, spans=spans,
             **{k: np.concatenate([p[k] for p in parts], axis=0) for k in parts[0]})


def load_tracks(out, chrom):
    z = np.load(os.path.join(out, "events", "windows", chrom + ".npz"))
    return {k: z[k] for k in z.files}
