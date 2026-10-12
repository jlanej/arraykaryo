"""Circular binary segmentation (Olshen 2004) on the window tracks of one chromosome, jointly over the tracks: for every
candidate segment (i, j] of the n windows the statistic is the sum over tracks of the squared t-statistic of the segment's
mean against the rest of the chromosome, each track standardised to unit noise; the best segment is accepted when its
statistic exceeds a chi-square quantile (3 degrees of freedom) at a level divided by the number of candidate pairs (a
Bonferroni-style bound in place of permutations, which the cohort scale forbids; calibrated on the mock's null samples,
docs/calibration), and the pieces are scanned again. Adjacent segments whose means do not differ are merged. Boundaries are
then refined to probe resolution by the least-squares split of the probe-level signal around each window boundary."""
import numpy as np

MIN_WINDOWS = 2


def chi2_3_sf(q):
    """P(chi-square with 3 df > q), closed form."""
    from math import erfc, exp, pi, sqrt
    q = max(float(q), 0.0)
    return erfc(sqrt(q / 2.0)) + sqrt(2.0 * q / pi) * exp(-q / 2.0)


def chi2_3_isf(p):
    """The q with P(chi-square 3 df > q) = p, by bisection."""
    lo, hi = 0.0, 400.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if chi2_3_sf(mid) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


_CRIT = {}


def critical(n, alpha, q_floor=20.0, df_tracks=3):
    """The statistic a segment must exceed on a chromosome of n windows: the 3-df chi-square quantile at alpha over the number
    of candidate (i, j) pairs, never below q_floor."""
    key = (int(n), float(alpha), float(q_floor))
    if key not in _CRIT:
        pairs = max(n * (n - 1) / 2.0, 1.0)
        _CRIT[key] = max(q_floor, chi2_3_isf(min(alpha / pairs, 0.5)))
    return _CRIT[key]


def best_segment(z, lo, hi, min_len=MIN_WINDOWS):
    """The (i, j, Q) of the segment z[:, lo+i : lo+j] whose means differ most from the rest of z[:, lo:hi] (tracks x windows,
    standardised; NaN treated as the track's mean, i.e. no evidence). Q is the sum over tracks of squared t-statistics."""
    x = z[:, lo:hi]
    x = np.where(np.isfinite(x), x, 0.0)
    T, n = x.shape
    if n < 2 * min_len:
        return 0, n, 0.0
    S = np.concatenate([np.zeros((T, 1)), np.cumsum(x, axis=1)], axis=1)           # T x (n+1)
    total = S[:, -1]
    if n <= 2600:                                                                   # every (i, j) at once: n^2 floats per track
        idx = np.arange(n + 1)
        k = idx[None, :] - idx[:, None]                                             # (i, j) -> j - i
        valid = (k >= min_len) & ((n - k) >= min_len)
        kf = np.where(valid, k, 1).astype(np.float64)
        mf = np.where(valid, n - k, 1).astype(np.float64)
        Q = np.zeros((n + 1, n + 1))
        for t in range(T):
            s_in = S[t][None, :] - S[t][:, None]
            d = s_in / kf - (total[t] - s_in) / mf
            Q += d * d / (1.0 / kf + 1.0 / mf)
        Q[~valid] = -1.0
        a = int(np.argmax(Q))
        i, j = divmod(a, n + 1)
        return int(i), int(j), float(Q[i, j])
    best = (0, n, 0.0)
    i_all = np.arange(n + 1)
    for j in range(min_len, n + 1):
        i = i_all[:j - min_len + 1]                                                 # segment (i, j] with length k >= min_len
        k = j - i
        m = n - k
        ok = m >= min_len
        if not ok.any():
            continue
        i, k, m = i[ok], k[ok], m[ok]
        s_in = S[:, j][:, None] - S[:, i]
        mean_in = s_in / k
        mean_out = (total[:, None] - s_in) / m
        t = (mean_in - mean_out) / np.sqrt(1.0 / k + 1.0 / m)
        Q = (t * t).sum(axis=0)
        a = int(np.argmax(Q))
        if Q[a] > best[2]:
            best = (int(i[a]), j, float(Q[a]))
    return best


def cbs(z, alpha=0.01, min_len=MIN_WINDOWS, q_floor=20.0, max_depth=40):
    """Breakpoints (window indices where a new segment starts, 0 excluded, n excluded) of the standardised tracks z (T x n)."""
    T, n = z.shape
    bps = set()
    stack = [(0, n, 0)]
    while stack:
        lo, hi, depth = stack.pop()
        if hi - lo < 2 * min_len or depth > max_depth:
            continue
        i, j, Q = best_segment(z, lo, hi, min_len)
        if Q < critical(hi - lo, alpha, q_floor):
            continue
        a, b = lo + i, lo + j
        pieces = []
        if a > lo:
            bps.add(a)
            pieces.append((lo, a))
        if b < hi:
            bps.add(b)
            pieces.append((b, hi))
        pieces.append((a, b))
        for p in pieces:
            stack.append((p[0], p[1], depth + 1))
    return sorted(bps)


def merge_segments(z, bounds, q_merge):
    """Merge adjacent segments (bounds: [b0, b1, ...] with b0 = 0 and the last = n) whose means do not differ: the joint statistic of
    the two segments' means (sum over tracks of squared t) under q_merge."""
    x = np.where(np.isfinite(z), z, 0.0)
    S = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(x, axis=1)], axis=1)
    bounds = list(bounds)
    while len(bounds) > 2:
        b = np.asarray(bounds)
        a, m, c = b[:-2], b[1:-1], b[2:]
        n1, n2 = (m - a).astype(float), (c - m).astype(float)
        m1 = (S[:, m] - S[:, a]) / n1
        m2 = (S[:, c] - S[:, m]) / n2
        Q = (((m1 - m2) ** 2) / (1.0 / n1 + 1.0 / n2)).sum(axis=0)
        k = int(np.argmin(Q))
        if Q[k] >= q_merge:
            break
        del bounds[k + 1]
    return bounds


def refine_boundary(values, masks, lo, hi, scales, pos=None, tie=1.0):
    """The probe index k in [lo, hi] at which splitting the probe-level tracks (each a 1-D array over probes with a mask of usable
    probes and a noise scale) into [lo, k) and [k, hi) reduces the residual sum of squares most, summed over the tracks;
    lo when nothing is usable. Among splits within `tie` (standardised units) of the best, the one at the largest gap between
    consecutive probes wins when `pos` is given: a boundary falling in a centromere or a probe desert is placed at the gap."""
    best_k, best_gain = lo, -np.inf
    n = hi - lo
    if n < 4:
        return lo
    gain = np.zeros(n + 1)
    any_track = False
    for v, m, s in zip(values, masks, scales):
        x = np.where(m[lo:hi] & np.isfinite(v[lo:hi]), v[lo:hi], np.nan)
        ok = np.isfinite(x)
        if ok.sum() < 4:
            continue
        any_track = True
        xv = np.where(ok, x, 0.0)
        c = np.concatenate([[0.0], np.cumsum(xv)])
        c2 = np.concatenate([[0.0], np.cumsum(xv * xv)])
        cn = np.concatenate([[0], np.cumsum(ok)])
        nl, nr = cn, cn[-1] - cn
        sl, sr = c, c[-1] - c
        with np.errstate(all="ignore"):
            rss_l = c2 - sl * sl / np.maximum(nl, 1)
            rss_r = (c2[-1] - c2) - sr * sr / np.maximum(nr, 1)
            total = c2[-1] - c[-1] ** 2 / max(cn[-1], 1)
            g = (total - rss_l - rss_r) / (s * s)
        g[(nl < 2) | (nr < 2)] = -np.inf
        gain += np.where(np.isfinite(g), g, -np.inf)
    if not any_track:
        return lo
    k = int(np.argmax(gain))
    if not np.isfinite(gain[k]):
        return lo
    if pos is not None:
        cand = np.flatnonzero(gain >= gain[k] - tie)
        cand = cand[(cand > 0) & (cand < n)]
        if cand.size:
            gaps = pos[lo + cand] - pos[lo + cand - 1]
            k = int(cand[np.argmax(gaps)])
    return lo + k
