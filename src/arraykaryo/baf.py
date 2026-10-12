"""The B-allele frequency's band deviation d: at heterozygous calls the BAF sits at 0.5 +- d (one band per homologue), and d
is the location of that symmetric two-band mixture, fitted by maximum likelihood with a shared band SD. The per-window mean
of |BAF - 0.5| is biased upward by noise (its expectation is sigma x sqrt(2/pi) when d = 0) and must not be used for a cell
fraction. The cell fraction follows from d: gain f = 4d / (1 - 2d), loss f = 4d / (1 + 2d), copy-neutral loss of
heterozygosity f = 2d (the band positions (1 + f) / (2 + f), 1 / (2 - f) and (1 + f) / 2 solved for f)."""
import numpy as np

D_GRID = np.linspace(0.0, 0.49, 246)       # step 0.002 in d: 0.008 in f for a gain, 0.004 for a copy-neutral LOH


def fold(baf):
    """|BAF - 0.5| after clipping to [0, 1]."""
    return np.abs(np.clip(baf, 0.0, 1.0) - 0.5)


def loglik_grid(u, sigma, grid=D_GRID):
    """Log-likelihood of the folded values u (finite) under 0.5 N(u; d, sigma) + 0.5 N(u; -d, sigma) for every d of the grid
    (a folded normal: u >= 0 carries both signs). Returns an array over the grid (the constant sqrt(2 pi) dropped; the
    -n log sigma term kept so that fits with different sigma compare)."""
    u = np.asarray(u, dtype=np.float64)
    if u.size == 0:
        return np.zeros(grid.size)
    s2 = 2.0 * sigma * sigma
    a = -((u[:, None] - grid[None, :]) ** 2) / s2
    b = -((u[:, None] + grid[None, :]) ** 2) / s2
    m = np.maximum(a, b)
    return (m + np.log(np.exp(a - m) + np.exp(b - m))).sum(axis=0) - u.size * np.log(sigma)


MAX_FIT = 4000      # heterozygous calls used by a fit (an evenly spaced subsample beyond: the SE of d at 4,000 calls is under 0.001)
COARSE = D_GRID[::5]


def fit_d(baf, sigma=None, grid=D_GRID):
    """The maximum-likelihood d of one set of heterozygous-call BAFs with a shared band SD. sigma given: the sample's band SD
    (from its own heterozygous calls); sigma None: fitted jointly on a grid. Returns (d, sigma, llr) with llr the log-likelihood
    ratio against d = 0 (0 when there are fewer than 4 values; scaled to the full count when a subsample was fitted). The
    search is coarse (every fifth grid point) then fine around the best."""
    u = fold(np.asarray(baf, dtype=np.float64))
    u = u[np.isfinite(u)]
    n_all = u.size
    if n_all < 4:
        return float("nan"), float("nan"), 0.0
    if n_all > MAX_FIT:
        u = u[:: int(np.ceil(n_all / float(MAX_FIT)))]
    scale = n_all / float(u.size)
    if sigma is None:
        best = (-np.inf, 0.0, 0.05)
        for s in (0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.13):
            ll = loglik_grid(u, s, COARSE)
            k = int(np.argmax(ll))
            if ll[k] > best[0]:
                best = (ll[k], float(COARSE[k]), s)
        _, d0, sigma = best
    else:
        ll = loglik_grid(u, sigma, COARSE)
        d0 = float(COARSE[int(np.argmax(ll))])
    fine = grid[(grid >= d0 - 0.012) & (grid <= d0 + 0.012)]
    if fine.size == 0 or fine[0] > 0.0:
        fine = np.concatenate([[0.0], fine])
    ll = loglik_grid(u, sigma, fine)
    k = int(np.argmax(ll))
    return float(fine[k]), float(sigma), float((ll[k] - ll[0]) * scale)


def d_em_windows(u, valid, sigma, iters=12):
    """d per window for many windows at once by EM with a known band SD: u (windows x width) folded values, valid the mask of
    heterozygous calls; sigma per window (broadcast). The E-step gives each value's probability of the upper band,
    r = sigmoid(2 d u / sigma^2); the M-step d = mean((2r - 1) u). Returns (d, n) with d NaN where fewer than 4 calls."""
    u = np.where(valid, u, 0.0)
    n = valid.sum(axis=1)
    sig2 = np.broadcast_to(np.asarray(sigma, dtype=np.float64)[..., None] if np.ndim(sigma) else sigma, u.shape) ** 2 if np.ndim(sigma) else sigma ** 2
    d = np.sqrt(np.maximum((u * u).sum(axis=1) / np.maximum(n, 1) - (sig2[:, 0] if np.ndim(sig2) else sig2), 0.0))      # the moment start
    for _ in range(iters):
        z = 2.0 * d[:, None] * u / sig2
        r = 1.0 / (1.0 + np.exp(-np.clip(z, -50, 50)))
        d = np.where(n > 0, ((2.0 * r - 1.0) * u * valid).sum(axis=1) / np.maximum(n, 1), 0.0)
    d = np.where(n >= 4, d, np.nan)
    return d, n


def f_gain(d):
    return 4.0 * d / np.maximum(1.0 - 2.0 * d, 1e-6)


def f_loss(d):
    return 4.0 * d / (1.0 + 2.0 * d)


def f_loh(d):
    return 2.0 * d


def f_from_d(d, kind):
    """Cell fraction from the band deviation for an event kind (gain, loss, loh); capped at 1."""
    if not np.isfinite(d):
        return float("nan")
    f = {"gain": f_gain, "loss": f_loss, "loh": f_loh}[kind](float(d))
    return float(min(f, 1.0))


def d_from_f(f, kind):
    """The band deviation an event in a cell fraction f produces (the inverse of f_from_d)."""
    f = float(f)
    if kind == "gain":
        return f / (2.0 * (2.0 + f))
    if kind == "loss":
        return f / (2.0 * (2.0 - f))
    return f / 2.0
