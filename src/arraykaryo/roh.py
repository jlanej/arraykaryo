"""Runs of homozygosity from the genotypes with PLINK's sliding window (Purcell 2007, --homozyg defaults as the starting
point): over the called SNV probes of a chromosome, a window of 50 probes is homozygous when it holds at most one
heterozygous call (genotyping errors allowed); a probe is in a run when at least 5% of the windows covering it are
homozygous; runs of consecutive such probes are kept when they span at least 1.5 Mb (the brief's floor; PLINK's 1 Mb)
with at least 100 probes, no gap above 1 Mb between consecutive probes and a density of at most one probe per 50 kb.
Per sample: the genome fraction in runs (F_ROH, over the autosomal probed span) and every run; runs of 10 Mb or more are
events (a whole chromosome without heterozygosity is an isodisomic uniparental disomy or consanguinity: both are said)."""
import numpy as np

ROH = dict(window=50, window_het=1, window_threshold=0.05, min_bp=1_500_000, min_probes=100, max_gap=1_000_000, max_kb_per_probe=50.0, event_bp=10_000_000)


def runs_of_homozygosity(pos, het, called, p=ROH):
    """Runs over one chromosome of one sample: pos (1-based positions of the probes used), het and called boolean arrays over them.
    Returns a list of (start, end, n_probes, n_het) with start/end positions."""
    sel = np.flatnonzero(called)
    if sel.size < p["min_probes"]:
        return []
    h = het[sel].astype(np.int32)
    x = pos[sel]
    m = sel.size
    w = p["window"]
    if m < w:
        return []
    cs = np.concatenate([[0], np.cumsum(h)])
    wh = cs[w:] - cs[:-w]                                            # het count of window k = probes k..k+w-1 (m-w+1 windows)
    hom = (wh <= p["window_het"]).astype(np.int32)
    ch = np.concatenate([[0], np.cumsum(hom)])
    i = np.arange(m)
    k0 = np.maximum(i - w + 1, 0)
    k1 = np.minimum(i, m - w) + 1
    n_cov = np.maximum(k1 - k0, 1)
    share = (ch[k1] - ch[k0]) / n_cov
    in_run = share >= p["window_threshold"]
    # runs of consecutive in-run probes, split at large gaps
    edges = np.flatnonzero(np.diff(np.concatenate([[0], in_run.astype(np.int8), [0]])))
    out = []
    for a, b in zip(edges[::2], edges[1::2]):
        xs = x[a:b]
        gaps = np.flatnonzero(np.diff(xs) > p["max_gap"])
        pieces = np.split(np.arange(a, b), gaps + 1)
        for pc in pieces:
            if pc.size < p["min_probes"]:
                continue
            s, e = int(x[pc[0]]), int(x[pc[-1]])
            if e - s + 1 < p["min_bp"]:
                continue
            if (e - s + 1) / 1000.0 / pc.size > p["max_kb_per_probe"]:
                continue
            out.append((s, e, int(pc.size), int(h[pc].sum())))
    return out
