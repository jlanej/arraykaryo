"""Figures. Self-explanatory: each image carries its title, axis labels and legend; a text sidecar (<figure>.txt: title,
caption, key) goes beside it for a manuscript. Okabe-Ito colours throughout (colour-blind safe), the same hue for the same
thing on every page: gain vermillion, loss blue, copy-neutral LOH reddish purple, run of homozygosity orange, LRR points dark
grey, BAF points light grey, fits black, thresholds dotted."""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PAL = dict(lrr="#444444", baf="#9a9a9a", fit="#000000", gain="#D55E00", loss="#0072B2", loh="#CC79A7", roh="#E69F00", het="#009E73", ref="#888888", band="#efefef",
           male="#0072B2", female="#D55E00", other="#009E73", flag="#CC79A7", pass_="#56B4E9", thr="#000000", mosaic="#E69F00")
TYPE_COLOUR = {"gain": PAL["gain"], "loss": PAL["loss"], "loh": PAL["loh"], "roh": PAL["roh"]}
TYPE_NAME = {"gain": "gain", "loss": "loss", "loh": "copy-neutral LOH", "roh": "run of homozygosity"}
MM = 1 / 25.4
DPI = 160
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "figure.dpi": DPI, "savefig.dpi": DPI,
                     "axes.spines.top": False, "axes.spines.right": False, "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"]})


def save(fig, path, title, caption, key=()):
    """PNG plus the sidecar: `<path>.png` and `<path>.txt`."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path + ".png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with open(path + ".txt", "w") as fh:
        fh.write("title: %s\ncaption: %s\n" % (title, caption))
        for k, v in key:
            fh.write("key: %s = %s\n" % (k, v))
    return path + ".png"


def _hist(ax, v, thr=None, thr2=None, label="", colour=PAL["pass_"], bins=30, log=False):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        ax.text(0.5, 0.5, "no values", ha="center", va="center", transform=ax.transAxes)
        return
    ax.hist(v, bins=bins, color=colour, edgecolor="white", linewidth=0.3)
    med = np.median(v)
    ax.axvline(med, color=PAL["fit"], lw=0.8, ls="-", label="median %.3g" % med)
    if thr is not None:
        ax.axvline(thr, color=PAL["gain"], lw=1.0, ls=":", label="flag %g" % thr)
    if thr2 is not None:
        ax.axvline(thr2, color=PAL["loh"], lw=1.0, ls="--", label="hard fail %g" % thr2)
    ax.set_xlabel(label)
    ax.set_ylabel("samples")
    ax.legend(frameon=False, loc="upper right")


def fig_qc(out, cache, rows, T, corr, log_fn=None):
    """qc/figures/qc_distributions: the sample metrics as histograms with the thresholds; qc/figures/gc_wave: LRR against the
    1-Mb GC before and after correction for a typical sample (the median wave factor) and the worst."""
    fdir = os.path.join(out, "qc", "figures")
    os.makedirs(fdir, exist_ok=True)
    get = lambda k: np.array([r[k] for r in rows], dtype=float)  # noqa: E731
    fig, axes = plt.subplots(2, 4, figsize=(260 * MM, 120 * MM))
    axes = axes.ravel()
    _hist(axes[0], get("call_rate"), thr=T["sample_call_rate"], label="genotype call rate (autosomal SNV probes)")
    _hist(axes[1], get("lrr_sd"), thr=T["sample_lrr_sd_flag"], thr2=T["sample_lrr_sd_fail"], label="autosomal LRR SD after correction")
    _hist(axes[2], get("lrr_sd_raw"), label="autosomal LRR SD before correction", colour=PAL["ref"])
    _hist(axes[3], get("baf_sd"), thr=T["sample_baf_sd"], label="BAF SD at heterozygous calls")
    drift = get("baf_drift")
    _hist(axes[4], drift, thr=T["sample_baf_drift"], label="BAF drift (share of homozygous calls\nwith BAF in 0.2-0.25 or 0.75-0.8)", bins=np.linspace(0, max(2 * T["sample_baf_drift"], np.nanmax(drift) * 1.1 if np.isfinite(np.nanmax(drift)) else 0.01), 30))
    _hist(axes[5], get("wf"), thr=T["sample_wf"], label="wave factor (PennCNV's: MAD of the 1-Mb\nwindow LRR medians, signed by GC)")
    axes[5].axvline(-T["sample_wf"], color=PAL["gain"], lw=1.0, ls=":")
    _hist(axes[6], get("lrr_mean"), thr=T["sample_lrr_mean"], label="autosomal LRR mean (raw)")
    axes[6].axvline(-T["sample_lrr_mean"], color=PAL["gain"], lw=1.0, ls=":")
    before, after = np.abs(get("gc_corr_before")), np.abs(get("gc_corr_after"))
    edges = np.linspace(0, max(np.nanmax(before) * 1.05, 0.05), 30)
    axes[7].hist(before[np.isfinite(before)], bins=edges, color=PAL["ref"], edgecolor="white", linewidth=0.3, label="before correction")
    axes[7].hist(after[np.isfinite(after)], bins=edges, color=PAL["pass_"], alpha=0.85, edgecolor="white", linewidth=0.3, label="after correction")
    axes[7].set_xlabel("|correlation of LRR with the 1-Mb GC|")
    axes[7].set_ylabel("samples")
    axes[7].legend(frameon=False, loc="upper right")
    fig.suptitle("Sample QC: %d samples, %d passing (dotted: flag threshold; dashed: hard fail)" % (len(rows), sum(int(r["pass"]) for r in rows)), y=1.01)
    fig.tight_layout()
    p1 = save(fig, os.path.join(fdir, "qc_distributions"), "Sample QC distributions",
              "Per-sample quality metrics over the autosomes with the thresholds that flag a sample (dotted) or fail it (dashed); the black line is the cohort median.",
              [("blue bars", "samples"), ("dotted red line", "flag threshold"), ("dashed purple line", "hard-fail threshold")])
    # the GC wave, before and after, for the median sample and the worst
    wf = get("wf")
    order = np.argsort(np.abs(np.nan_to_num(wf)))
    picks = [("typical", int(order[len(order) // 2])), ("strongest wave", int(order[-1]))]
    chrom = next((c for c in cache.chroms if c not in ("chrX", "chrY")), cache.chroms[0])
    gc_z, w_z = corr.covariates(cache, chrom)
    w = w_z * corr.gc_stats["w_sd"] + corr.gc_stats["w_mean"]
    fig, axes = plt.subplots(2, 2, figsize=(180 * MM, 120 * MM), sharex=True)
    for col, (what, i) in enumerate(picks):
        raw = np.asarray(cache.lrr(chrom)[i], dtype=float)
        cor = corr.apply(cache, chrom, sample_idx=[i])[0]
        for row, (y, lab) in enumerate(((raw, "before"), (cor, "after"))):
            ax = axes[row, col]
            ok = np.isfinite(y) & np.isfinite(w)
            step = max(1, ok.sum() // 20000)
            ax.scatter(w[ok][::step], y[ok][::step], s=1.5, c=PAL["lrr"], lw=0, alpha=0.5, rasterized=True)
            bins = np.linspace(np.nanpercentile(w[ok], 1), np.nanpercentile(w[ok], 99), 25)
            k = np.digitize(w[ok], bins)
            med = [np.median(y[ok][k == j]) if (k == j).sum() > 20 else np.nan for j in range(1, len(bins))]
            ax.plot((bins[:-1] + bins[1:]) / 2, med, color=PAL["gain"], lw=1.5, label="median per GC bin")
            ax.axhline(0, color=PAL["ref"], lw=0.6, ls=":")
            ax.set_ylim(-1, 1)
            ax.set_title("%s sample (%s): LRR %s correction, %s" % (what, rows[i]["sample"], lab, chrom))
            ax.set_ylabel("LRR")
            if row == 1:
                ax.set_xlabel("GC fraction of the 1-Mb window around the probe")
            ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    p2 = save(fig, os.path.join(fdir, "gc_wave"), "The GC wave before and after correction",
              "Raw and corrected LRR of two samples against the GC fraction of the 1-Mb window around each probe; the red line is the median per GC bin. The correction regresses LRR on the probe's GC and its window GC (quadratic) per sample and subtracts the fit, then the cohort reference.",
              [("grey points", "probes"), ("red line", "median LRR per GC bin")])
    if log_fn:
        log_fn("figures: %s, %s" % (os.path.relpath(p1, out), os.path.relpath(p2, out)), "qc")
    return [p1, p2]


def fig_sex(out, rows, summary, log_fn=None, max_labels=40):
    """sex/figures/x_vs_y: the X level against the Y level, coloured by the X heterozygosity relative to the autosomes, every
    sample not plainly XX or XY ringed and labelled (up to max_labels; the rest are in sex.tsv)."""
    fdir = os.path.join(out, "sex", "figures")
    os.makedirs(fdir, exist_ok=True)
    x = np.array([r["x_lrr"] for r in rows], dtype=float)
    y = np.array([r["y_lrr"] for r in rows], dtype=float)
    h = np.array([r["x_het_rel"] for r in rows], dtype=float)
    comp = [r["complement"] for r in rows]
    odd = np.array([(r["mosaic"] != "NA") or (r["complement"] not in ("XX", "XY")) or (r["check"] not in ("agrees", "NA", "not read", "")) for r in rows])
    fig, ax = plt.subplots(figsize=(160 * MM, 115 * MM))
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() == 0:
        ax.text(0.5, 0.5, "no X and Y levels", ha="center", va="center", transform=ax.transAxes)
    else:
        sc = ax.scatter(x[ok], y[ok], c=np.clip(np.nan_to_num(h[ok]), 0, 1.2), cmap="viridis", s=18, lw=0.3, edgecolor="white", vmin=0, vmax=1.2, zorder=3)
        cb = fig.colorbar(sc, ax=ax, pad=0.02)
        cb.set_label("X heterozygosity / autosomal heterozygosity")
        if (odd & ok).any():
            ax.scatter(x[odd & ok], y[odd & ok], s=80, facecolor="none", edgecolor=PAL["gain"], lw=1.2, label="not plainly XX or XY (ringed)", zorder=4)
        xm, ym = summary["x_modes"], summary["y_modes"]
        xlo, xhi = np.nanmin(x[ok]), np.nanmax(x[ok])
        ylo, yhi = np.nanmin(y[ok]), np.nanmax(y[ok])
        ax.set_xlim(xlo - 0.12, xhi + 0.12)
        ax.set_ylim(ylo - 0.5, yhi + 0.6)
        for v, lab in ((xm["one_copy"], "one X"), (xm["two_copies"], "two X's")):
            if v is not None and np.isfinite(v):
                ax.axvline(v, color=PAL["ref"], lw=0.6, ls=":", zorder=1)
                ax.text(v, ylo - 0.45, lab, ha="center", va="bottom", fontsize=7, color=PAL["ref"], bbox=dict(fc="white", ec="none", pad=1))
        for v, lab in ((ym["present"], "Y present"), (ym["absent"], "Y absent")):
            if v is not None and np.isfinite(v):
                ax.axhline(v, color=PAL["ref"], lw=0.6, ls=":", zorder=1)
                ax.text(xlo - 0.11, v, lab, ha="left", va="bottom", fontsize=7, color=PAL["ref"], bbox=dict(fc="white", ec="none", pad=1))
        n_lab = 0
        for i in np.flatnonzero(odd & ok):
            if n_lab >= max_labels:
                break
            kar = rows[i]["karyotype"] if rows[i]["karyotype"] != "NA" else comp[i]
            rep = rows[i]["reported_sex"]
            lab = "%s %s" % (rows[i]["sample"], kar) + ((" (reported %s)" % rep) if rep != "NA" and rows[i]["check"] != "agrees" else "")
            ax.annotate(lab, (x[i], y[i]), xytext=(5, 5 if (n_lab % 2 == 0) else -9), textcoords="offset points", fontsize=6, color=PAL["lrr"], zorder=5)
            n_lab += 1
        if odd.sum() > max_labels:
            ax.text(0.01, 0.99, "%d of %d ringed samples labelled; all are in sex.tsv" % (max_labels, odd.sum()), transform=ax.transAxes, fontsize=6, va="top")
        ax.legend(frameon=False, loc="center right")
    xm, ym = summary["x_modes"], summary["y_modes"]
    ax.set_xlabel("X level: median corrected LRR outside the PARs")
    ax.set_ylabel("Y level: median LRR over the Y-specific probes")
    ax.set_title("Sex chromosomes of %d samples (X modes %.2f, %.2f; Y modes %s, %s)" % (len(rows), xm["one_copy"], xm["two_copies"], _f2(ym["present"]), _f2(ym["absent"])))
    p = save(fig, os.path.join(fdir, "x_vs_y"), "X level against Y level",
             "Each point is a sample: the median corrected LRR on the X outside the pseudoautosomal regions against the median LRR over the Y-specific probes, coloured by the X heterozygosity relative to the autosomes. Dotted lines mark the cohort's modes (one and two X copies; Y present and absent). Samples not plainly XX or XY are ringed and labelled with their complement, cell fraction and, where it disagrees, the reported sex.",
             [("colour", "X heterozygosity / autosomal heterozygosity (0 one X, 1 two X's)"), ("orange ring", "a mosaic, an aneuploidy, or a complement against the reported sex"), ("dotted lines", "the cohort's X and Y modes")])
    if log_fn:
        log_fn("figure: " + os.path.relpath(p, out), "sex")
    return p


def _f2(v):
    return "%.2f" % v if (v is not None and np.isfinite(v)) else "NA"


def fig_sex_tracks(out, C, corr, masks, row):
    """sex/figures/tracks_<sample>: the sample's X (corrected LRR, BAF) and Y (LRR) along the chromosomes."""
    fdir = os.path.join(out, "sex", "figures")
    name = row["sample"]
    i = C.sample_index(name)
    have_y = "chrY" in C.chroms
    fig, axes = plt.subplots(3 if have_y else 2, 1, figsize=(180 * MM, (105 if have_y else 75) * MM))
    axes = np.atleast_1d(axes)
    prx = C.probes("chrX")
    L = corr.apply(C, "chrX", sample_idx=[i])[0]
    B = np.asarray(C.baf("chrX")[i], dtype=float)
    g = np.asarray(C.gt("chrX")[i])
    keep = masks["chrX"]["keep"]
    pos = prx.pos / 1e6
    ax = axes[0]
    ax.scatter(pos[keep], L[keep], s=1.2, c=PAL["lrr"], lw=0, rasterized=True)
    med = np.nanmedian(L[keep & ~prx.par])
    ax.axhline(0, color=PAL["ref"], lw=0.6, ls=":")
    ax.axhline(med, color=PAL["gain"], lw=1.0, label="median X level %.2f" % med)
    for a, b in C.genome.par.get("chrX", ()):
        ax.axvspan(a / 1e6, b / 1e6, color=PAL["band"], lw=0, label="PAR" if a < 1e7 else None)
    ax.set_ylim(-1.5, 1.0)
    ax.set_ylabel("X LRR (corrected)")
    ax.set_title("%s: %s; %s; %s" % (name, row["karyotype"] if row["karyotype"] != "NA" else row["complement"], row["mosaic"] if row["mosaic"] != "NA" else "no mosaic", row["check"]))
    ax.legend(frameon=False, loc="upper right", ncol=3)
    ax = axes[1]
    het = (g == 1) & keep
    ax.scatter(pos[keep & (g != 1)], np.clip(B[keep & (g != 1)], -0.05, 1.05), s=1.0, c=PAL["baf"], lw=0, rasterized=True, label="homozygous / no call")
    ax.scatter(pos[het], np.clip(B[het], -0.05, 1.05), s=2.0, c=PAL["het"], lw=0, rasterized=True, label="heterozygous call")
    for a, b in C.genome.par.get("chrX", ()):
        ax.axvspan(a / 1e6, b / 1e6, color=PAL["band"], lw=0)
    ax.axhline(0.5, color=PAL["ref"], lw=0.6, ls=":")
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("X BAF")
    ax.set_xlabel("position on chrX (Mb)")
    ax.legend(frameon=False, loc="upper right", ncol=2, markerscale=4)
    if have_y:
        pry = C.probes("chrY")
        Ly = corr.apply(C, "chrY", sample_idx=[i], reference=False)[0]
        ky = masks["chrY"]["keep"] & ~pry.par
        ax = axes[2]
        ax.scatter(pry.pos[ky] / 1e6, np.clip(Ly[ky], -6, 2), s=1.5, c=PAL["lrr"], lw=0, rasterized=True)
        ax.axhline(0, color=PAL["ref"], lw=0.6, ls=":")
        ax.axhline(row["y_lrr"], color=PAL["loss"], lw=1.0, label="Y level %.2f (%.2f copies)" % (float(row["y_lrr"]), float(row["y_copies"])))
        ax.set_ylim(-6, 2)
        ax.set_ylabel("Y LRR (GC-corrected)")
        ax.set_xlabel("position on chrY (Mb)")
        ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    return save(fig, os.path.join(fdir, "tracks_%s" % name), "X and Y tracks of %s" % name,
                "The corrected LRR and the BAF along the X (pseudoautosomal regions shaded) and the LRR along the Y of one sample, with its X level, its Y level and its complement.",
                [("dark grey", "LRR per probe"), ("green", "BAF at heterozygous calls"), ("light grey", "BAF at homozygous calls and no-calls"), ("orange line", "the median X level"), ("blue line", "the median Y level")])


def fig_kinship(out, kin, ibs0, cls, p, po_cut, log_fn=None):
    """pedigree/figures/kinship_vs_ibs0: every pair, the related ones coloured by class, the boundaries drawn."""
    fdir = os.path.join(out, "pedigree", "figures")
    os.makedirs(fdir, exist_ok=True)
    kin, ibs0 = np.asarray(kin, dtype=float), np.asarray(ibs0, dtype=float)
    cls = np.asarray(cls)
    ok = np.isfinite(kin) & np.isfinite(ibs0)
    fig, ax = plt.subplots(figsize=(150 * MM, 105 * MM))
    colours = {"": PAL["baf"], "3rd": "#999999", "2nd": PAL["het"], "fs": PAL["loss"], "po": PAL["gain"], "dup": PAL["loh"]}
    names = {"": "unrelated", "3rd": "third degree", "2nd": "second degree", "fs": "full siblings", "po": "parent-offspring", "dup": "duplicate / monozygotic"}
    un = ok & (cls == "")
    step = max(1, un.sum() // 50000)
    ax.scatter(ibs0[un][::step], kin[un][::step], s=3, c=colours[""], lw=0, rasterized=True, label="unrelated (%d pairs)" % un.sum())
    for c in ("3rd", "2nd", "fs", "po", "dup"):
        m = ok & (cls == c)
        if m.any():
            ax.scatter(ibs0[m], kin[m], s=14, c=colours[c], lw=0.3, edgecolor="white", label="%s (%d)" % (names[c], m.sum()), zorder=3)
    for v, lab in ((p["dup"], "duplicate"), (p["first"], "first degree"), (p["second"], "second degree"), (p["third"], "third degree")):
        ax.axhline(v, color=PAL["ref"], lw=0.6, ls=":")
        ax.text(ax.get_xlim()[1] if False else np.nanmax(ibs0[ok]) * 1.02, v, lab, fontsize=6, va="bottom", ha="right", color=PAL["ref"])
    ax.axvline(po_cut, color=PAL["ref"], lw=0.6, ls="--")
    ax.text(po_cut, p["dup"] + 0.02, "parent-offspring | full siblings", fontsize=6, ha="left", va="bottom", color=PAL["ref"])
    ax.set_xlabel("IBS0: share of probes at which the pair is opposite homozygous")
    ax.set_ylabel("KING-robust kinship")
    ax.set_ylim(min(-0.1, np.nanmin(kin[ok]) - 0.02), 0.56)
    ax.set_title("Relatedness of every pair of samples")
    ax.legend(frameon=False, loc="upper right", markerscale=1.5)
    pth = save(fig, os.path.join(fdir, "kinship_vs_ibs0"), "Kinship against IBS0",
               "Every pair of samples: the KING-robust kinship against the share of probes at which the two are opposite homozygotes (IBS0). Dotted lines are the class boundaries (0.354 duplicate or monozygotic; 0.177 first degree; 0.0884 second; 0.0442 third); the dashed line splits first-degree pairs into parent-offspring (IBS0 at the genotyping-error floor) and full siblings.",
               [(names[c], colours[c]) for c in ("", "3rd", "2nd", "fs", "po", "dup")])
    if log_fn:
        log_fn("figure: " + os.path.relpath(pth, out), "pedigree")
    return pth


# ------------------------------------------------------------------------------------------------------- sample figures

STAIN = {"gneg": "#ffffff", "gpos25": "#c8c8c8", "gpos50": "#9a9a9a", "gpos75": "#6a6a6a", "gpos100": "#2a2a2a", "acen": "#b03a2e", "gvar": "#d7d7e0", "stalk": "#8da0cb"}


def _ideogram(ax, G, chrom, x0=0.0, scale=1e6, height=1.0):
    """The cytogenetic bands of one chromosome as a strip on ax (x in Mb from x0)."""
    for s, e, name, stain in G.bands.get(chrom, []):
        ax.add_patch(plt.Rectangle((x0 + s / scale, 0), (e - s) / scale, height, facecolor=STAIN.get(stain, "#ffffff"), edgecolor="#555555", lw=0.3))
    ax.set_ylim(0, height)
    ax.set_yticks([])


def _event_label(e):
    f = float(e["f"]) if e["f"] not in ("NA", "", None) else float("nan")
    cells = "" if (not np.isfinite(f) or f >= 0.9) else " %.0f%%" % (100 * f)
    where = e["chrom"][3:] + (e["span"] if e["span"] in ("p", "q") else "")
    return "%s %s%s" % (TYPE_NAME.get(e["type"], e["type"]), where, cells)


def _thin(n, max_points):
    return max(1, int(np.ceil(n / float(max_points))))


def fig_genome(out, ctx, name, events, rohs, kar, max_points=40000):
    """events/figures/<sample>/genome: every chromosome in one image - the corrected LRR with the segment levels, the BAF, the
    heterozygous ratio per window with the runs of homozygosity - events drawn as coloured bars with labels."""
    C, G, corr, masks = ctx.C, ctx.G, ctx.corr, ctx.masks
    i = C.sample_index(name)
    chroms = [c for c in C.chroms if c in ctx.tracks or c == "chrY"]
    off, total = {}, 0.0
    for c in chroms:
        off[c] = total
        total += G.length[c] / 1e6
    fig, axes = plt.subplots(4, 1, figsize=(260 * MM, 120 * MM), sharex=True, gridspec_kw=dict(height_ratios=[0.25, 1.1, 1.0, 0.7], hspace=0.08))
    ax_ev, ax_l, ax_b, ax_h = axes
    xs_l, ys_l, xs_b, ys_b, hb = [], [], [], [], []
    q = ctx.qc[name]
    sx = ctx.sex.get(name, {})
    xy = sx.get("base", "") == "XY"
    for c in chroms:
        pr = C.probes(c)
        m = masks[c]
        keep = m["keep"] & (~pr.par if c in ("chrX", "chrY") else True)
        if c == "chrY":
            L = corr.apply(C, c, sample_idx=[i], reference=False)[0]
            keep &= ~pr.par
        else:
            L = corr.apply(C, c, sample_idx=[i])[0]
        step = _thin(int(keep.sum()), max_points // len(chroms))
        idx = np.flatnonzero(keep)[::step]
        xs_l.append(off[c] + pr.pos[idx] / 1e6)
        ys_l.append(L[idx])
        g = np.asarray(C.gt(c)[i])
        b = np.asarray(C.baf(c)[i], dtype=np.float32)
        kg = m["keep_gt"] & (g >= 0) & np.isfinite(b)
        idb = np.flatnonzero(kg)[::_thin(int(kg.sum()), max_points // len(chroms))]
        xs_b.append(off[c] + pr.pos[idb] / 1e6)
        ys_b.append(np.clip(b[idb], -0.02, 1.02))
        hb.append(g[idb] == 1)
    x_l, y_l = np.concatenate(xs_l), np.concatenate(ys_l)
    ax_l.scatter(x_l, np.clip(y_l, -1.5, 1.2), s=1.0, c=PAL["lrr"], lw=0, rasterized=True, alpha=0.6)
    x_b, y_b, het = np.concatenate(xs_b), np.concatenate(ys_b), np.concatenate(hb)
    ax_b.scatter(x_b[~het], y_b[~het], s=0.8, c=PAL["baf"], lw=0, rasterized=True, alpha=0.5)
    ax_b.scatter(x_b[het], y_b[het], s=1.2, c=PAL["het"], lw=0, rasterized=True, alpha=0.8)
    # window heterozygous ratio
    for c in chroms:
        if c not in ctx.tracks:
            continue
        lay = ctx.layout[c]
        hr = np.asarray(ctx.tracks[c]["het_ratio"][i], dtype=float)
        ax_h.plot(off[c] + lay["mids"] / 1e6, np.clip(hr, 0, 1.6), color=PAL["het"], lw=0.7)
    for k, c in enumerate(chroms):
        for ax in (ax_l, ax_b, ax_h):
            if k % 2:
                ax.axvspan(off[c], off[c] + G.length[c] / 1e6, color=PAL["band"], lw=0, zorder=0)
            cen = G.centromere[c]
            if cen[1]:
                ax.axvline(off[c] + cen[0] / 1e6, color="#bbbbbb", lw=0.4, ls=":", zorder=0)
    one = ctx.levels[0]
    ax_l.axhline(0, color=PAL["ref"], lw=0.6, ls=":")
    ax_l.axhline(one, color=PAL["loss"], lw=0.5, ls=":")
    ax_l.axhline(ctx.levels[1], color=PAL["gain"], lw=0.5, ls=":")
    ax_l.text(total, one, " 1 copy", fontsize=6, color=PAL["loss"], va="center")
    ax_l.text(total, ctx.levels[1], " 3 copies", fontsize=6, color=PAL["gain"], va="center")
    for v in (0.5,):
        ax_b.axhline(v, color=PAL["ref"], lw=0.5, ls=":")
    ax_h.axhline(1.0, color=PAL["ref"], lw=0.5, ls=":")
    # events and runs
    for e in events:
        c = e["chrom"]
        if c not in off:
            continue
        x0, x1 = off[c] + int(e["start"]) / 1e6, off[c] + int(e["end"]) / 1e6
        col = TYPE_COLOUR.get(e["type"], PAL["lrr"])
        y = {"gain": 0.78, "loss": 0.55, "loh": 0.32, "roh": 0.1}.get(e["type"], 0.5)
        ax_ev.plot([x0, x1], [y, y], color=col, lw=4, solid_capstyle="butt")
        ax_ev.text((x0 + x1) / 2, y + 0.1, _event_label(e), fontsize=5.5, ha="center", va="bottom", color=col)
        lv = float(e["lrr"]) if e["lrr"] not in ("NA", "") and e["source"] != "sex" else np.nan
        if np.isfinite(lv) and e["type"] in ("gain", "loss"):
            ax_l.plot([x0, x1], [lv, lv], color=col, lw=2.0, solid_capstyle="butt")
        elif e["type"] in ("loh", "roh"):
            ax_l.plot([x0, x1], [0, 0], color=col, lw=2.0, solid_capstyle="butt")
        for ax in (ax_l, ax_b, ax_h):
            ax.axvspan(x0, x1, color=col, alpha=0.08, lw=0, zorder=0)
    for r in rohs:
        c = r["chrom"]
        if c in off:
            ax_h.plot([off[c] + int(r["start"]) / 1e6, off[c] + int(r["end"]) / 1e6], [0.03, 0.03], color=PAL["roh"], lw=3, solid_capstyle="butt")
    ax_ev.set_ylim(0, 1.05)
    ax_ev.set_yticks([])
    for sp in ("left", "bottom"):
        ax_ev.spines[sp].set_visible(False)
    ax_l.set_ylim(-1.5, 1.2)
    ax_l.set_ylabel("LRR (corrected)")
    ax_b.set_ylim(-0.02, 1.02)
    ax_b.set_ylabel("BAF")
    ax_h.set_ylim(0, 1.6)
    ax_h.set_ylabel("het ratio")
    ax_h.set_xlim(0, total)
    ax_h.set_xticks([off[c] + G.length[c] / 2e6 for c in chroms])
    ax_h.set_xticklabels([c[3:] for c in chroms], fontsize=6)
    ax_h.set_xlabel("chromosome")
    title = "%s: %s%s" % (name, kar.get("karyotype", "") if kar else "", ("; QC flags: " + q["flags"]) if q.get("flags", "NA") != "NA" else "")
    ax_ev.set_title(title, fontsize=9, loc="left")
    handles = [plt.Line2D([], [], color=TYPE_COLOUR[t], lw=4, label=TYPE_NAME[t]) for t in ("gain", "loss", "loh", "roh")]
    handles += [plt.Line2D([], [], marker="o", ls="", color=PAL["het"], markersize=3, label="BAF at heterozygous calls"), plt.Line2D([], [], marker="o", ls="", color=PAL["baf"], markersize=3, label="BAF at homozygous calls"),
                plt.Line2D([], [], color=PAL["het"], lw=1, label="heterozygous rate / expected (per window)"), plt.Line2D([], [], color=PAL["roh"], lw=3, label="run of homozygosity")]
    ax_h.legend(handles=handles, frameon=False, fontsize=5.5, loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.45))
    fdir = os.path.join(out, "events", "figures", name)
    return save(fig, os.path.join(fdir, "genome"), "%s, every chromosome" % name,
                "Corrected LRR (dotted lines: two copies, one copy, three copies on the cohort's scale), BAF (green at heterozygous calls) and the heterozygous rate against expectation per window; events as bars above with their type and cell fraction, shaded across the rows; runs of homozygosity as orange bars. %s" % title,
                [("vermillion", "gain"), ("blue", "loss"), ("reddish purple", "copy-neutral LOH"), ("orange", "run of homozygosity"), ("green", "heterozygous calls / rate")])


def fig_chrom(out, ctx, name, chrom, events, rohs, kar, max_points=12000):
    """events/figures/<sample>/chr_<chrom>: one chromosome - ideogram, LRR with the segment levels, BAF with the expected bands of the
    called cell fraction, heterozygous ratio with the runs of homozygosity."""
    from .baf import d_from_f
    C, G, corr, masks = ctx.C, ctx.G, ctx.corr, ctx.masks
    i = C.sample_index(name)
    pr = C.probes(chrom)
    m = masks[chrom]
    is_y = chrom == "chrY"
    keep = m["keep"] & (~pr.par if chrom in ("chrX", "chrY") else True)
    L = corr.apply(C, chrom, sample_idx=[i], reference=not is_y)[0]
    g = np.asarray(C.gt(chrom)[i])
    b = np.asarray(C.baf(chrom)[i], dtype=np.float32)
    fig, axes = plt.subplots(4, 1, figsize=(200 * MM, 125 * MM), sharex=True, gridspec_kw=dict(height_ratios=[0.18, 1.1, 1.0, 0.7], hspace=0.1))
    ax_i, ax_l, ax_b, ax_h = axes
    _ideogram(ax_i, G, chrom)
    ax_i.set_title("%s, %s: %s" % (name, chrom, "; ".join(_event_label(e) + (" (%s)" % e["bands"] if e.get("bands") not in ("", "NA", None) else "") for e in events) or "no event"), fontsize=9, loc="left")
    idx = np.flatnonzero(keep)[::_thin(int(keep.sum()), max_points)]
    pos = pr.pos / 1e6
    ax_l.scatter(pos[idx], np.clip(L[idx], -1.5 if not is_y else -6, 1.2), s=1.5, c=PAL["lrr"], lw=0, rasterized=True, alpha=0.7)
    if not is_y:
        ax_l.axhline(0, color=PAL["ref"], lw=0.6, ls=":")
        ax_l.axhline(ctx.levels[0], color=PAL["loss"], lw=0.5, ls=":")
        ax_l.axhline(ctx.levels[1], color=PAL["gain"], lw=0.5, ls=":")
        ax_l.set_ylim(-1.5, 1.2)
    else:
        ax_l.set_ylim(-6, 1.5)
    kg = m["keep_gt"] & (g >= 0) & np.isfinite(b)
    het = kg & (g == 1)
    idb = np.flatnonzero(kg)[::_thin(int(kg.sum()), max_points)]
    hh = g[idb] == 1
    ax_b.scatter(pos[idb][~hh], np.clip(b[idb][~hh], -0.02, 1.02), s=1.2, c=PAL["baf"], lw=0, rasterized=True, alpha=0.5)
    ax_b.scatter(pos[idb][hh], np.clip(b[idb][hh], -0.02, 1.02), s=2.5, c=PAL["het"], lw=0, rasterized=True, alpha=0.9)
    ax_b.axhline(0.5, color=PAL["ref"], lw=0.5, ls=":")
    if chrom in ctx.tracks:
        lay = ctx.layout[chrom]
        hr = np.asarray(ctx.tracks[chrom]["het_ratio"][i], dtype=float)
        ax_h.plot(lay["mids"] / 1e6, np.clip(hr, 0, 1.6), color=PAL["het"], lw=0.9)
        ax_h.axhline(1.0, color=PAL["ref"], lw=0.5, ls=":")
    for e in events:
        x0, x1 = int(e["start"]) / 1e6, int(e["end"]) / 1e6
        col = TYPE_COLOUR.get(e["type"], PAL["lrr"])
        f = float(e["f"]) if e["f"] not in ("NA", "") else np.nan
        lv = float(e["lrr"]) if e["lrr"] not in ("NA", "") and e["source"] != "sex" else np.nan
        if np.isfinite(lv) and e["type"] in ("gain", "loss"):
            ax_l.plot([x0, x1], [lv, lv], color=col, lw=2.5, solid_capstyle="butt", label="segment median LRR %.3f" % lv)
        elif e["type"] in ("loh", "roh") and not is_y:
            ax_l.plot([x0, x1], [0, 0], color=col, lw=2.5, solid_capstyle="butt")
        if np.isfinite(f) and e["type"] in ("gain", "loss", "loh") and f < 1.0 and not is_y:
            d = d_from_f(f, e["type"])
            for y in (0.5 - d, 0.5 + d):
                ax_b.plot([x0, x1], [y, y], color=col, lw=1.2, ls="--", label="expected bands at f=%.2f" % f if y > 0.5 else None)
        elif e["type"] == "gain" and np.isfinite(f) and f >= 1.0 and not is_y:
            for y in (1 / 3, 2 / 3):
                ax_b.plot([x0, x1], [y, y], color=col, lw=1.2, ls="--", label="expected bands at 3 copies" if y > 0.5 else None)
        for ax in (ax_l, ax_b, ax_h):
            ax.axvspan(x0, x1, color=col, alpha=0.1, lw=0, zorder=0)
    for r in rohs:
        ax_h.plot([int(r["start"]) / 1e6, int(r["end"]) / 1e6], [0.03, 0.03], color=PAL["roh"], lw=4, solid_capstyle="butt", label="run of homozygosity (%.1f Mb)" % (int(r["size_bp"]) / 1e6))
    for ax in (ax_l, ax_b, ax_h):
        for a, bb in G.par.get(chrom, ()):
            ax.axvspan(a / 1e6, bb / 1e6, color=PAL["band"], lw=0, zorder=0)
    ax_l.set_ylabel("LRR" + (" (GC-corrected, raw level)" if is_y else " (corrected)"))
    ax_b.set_ylim(-0.02, 1.02)
    ax_b.set_ylabel("BAF")
    ax_h.set_ylim(0, 1.6)
    ax_h.set_ylabel("het ratio")
    ax_h.set_xlabel("position on %s (Mb)" % chrom)
    ax_h.set_xlim(0, G.length[chrom] / 1e6)
    for ax in (ax_l, ax_b, ax_h):
        h, l = ax.get_legend_handles_labels()
        if h:
            ax.legend(h[:4], l[:4], frameon=False, fontsize=6, loc="upper right")
    fdir = os.path.join(out, "events", "figures", name)
    return save(fig, os.path.join(fdir, "chr_%s" % chrom), "%s, %s" % (name, chrom),
                "Ideogram, corrected LRR with the called segment's median level (dotted: two, one and three copies on the cohort's scale), BAF with the bands expected at the called cell fraction (dashed), the heterozygous rate per window against expectation with runs of homozygosity. Events shaded; pseudoautosomal regions in grey.",
                [("vermillion", "gain"), ("blue", "loss"), ("reddish purple", "copy-neutral LOH"), ("orange", "run of homozygosity"), ("green", "heterozygous calls / rate")])


def _figures_worker(args):
    out, cache_dir, names = args["out"], args["cache"], args["names"]
    from .events import SampleContext, read_events, read_karyotypes, read_rohs
    from .util import read_json
    ctx = SampleContext(out, cache_dir)
    summ = read_json(os.path.join(out, "events", "summary.json"))
    cal = summ["calibration"]
    ctx.levels = (cal["one_copy_lrr"], cal["three_copies_lrr"])
    ev_all, roh_all, kar_all = read_events(out), read_rohs(out), read_karyotypes(out)
    n = 0
    for name in names:
        evs = [e for e in ev_all if e["sample"] == name]
        rohs = [r for r in roh_all if r["sample"] == name]
        kar = kar_all.get(name, {})
        fig_genome(out, ctx, name, evs, rohs, kar)
        for c in sorted({e["chrom"] for e in evs}, key=lambda c: ctx.G.chrom_order().get(c, 99)):
            fig_chrom(out, ctx, name, c, [e for e in evs if e["chrom"] == c], [r for r in rohs if r["chrom"] == c], kar)
        n += 1
    return n


def draw_sample_figures(out, cache_dir, names, threads=1, log_fn=None):
    """The genome figure and one figure per chromosome with an event for every named sample; samples in parallel."""
    from concurrent.futures import ProcessPoolExecutor
    if not names:
        return 0
    k = max(1, int(np.ceil(len(names) / float(max(threads, 1)))))
    jobs = [dict(out=out, cache=cache_dir, names=names[a:a + k]) for a in range(0, len(names), k)]
    if threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            return sum(ex.map(_figures_worker, jobs))
    return sum(_figures_worker(j) for j in jobs)


# -------------------------------------------------------------------------------------------------------- cohort figures

def fig_landscape(out, events, G, samples_order, max_rows=120):
    """report/figures/landscape: every chromosome across, one row per sample with an event; each event a bar coloured by its
    type, thick when constitutional, thin when mosaic; above, the number of events per chromosome by type."""
    fdir = os.path.join(out, "report", "figures")
    os.makedirs(fdir, exist_ok=True)
    chroms = [c for c in G.chroms if c != "chrM"]
    off, total = {}, 0.0
    for c in chroms:
        off[c] = total
        total += G.length[c] / 1e6
    with_ev = [s for s in samples_order if any(e["sample"] == s for e in events)]
    shown = with_ev[:max_rows]
    fig, (ax_c, ax) = plt.subplots(2, 1, figsize=(260 * MM, (45 + 2.2 * max(len(shown), 8)) * MM), gridspec_kw=dict(height_ratios=[1.2, max(len(shown), 8) / 10.0], hspace=0.12), sharex=True)
    counts = {t: np.zeros(len(chroms)) for t in TYPE_COLOUR}
    for e in events:
        if e["chrom"] in off:
            counts[e["type"]][chroms.index(e["chrom"])] += 1
    width = np.array([G.length[c] / 1e6 for c in chroms]) * 0.8
    base = np.zeros(len(chroms))
    centres = np.array([off[c] + G.length[c] / 2e6 for c in chroms])
    for t in ("gain", "loss", "loh", "roh"):
        ax_c.bar(centres, counts[t], width=width, bottom=base, color=TYPE_COLOUR[t], label=TYPE_NAME[t], lw=0)
        base += counts[t]
    ax_c.set_ylabel("events")
    ax_c.legend(frameon=False, ncol=4, fontsize=6, loc="upper right")
    ax_c.set_title("Events per chromosome, and one row per sample with an event (%d of %d samples with events shown)" % (len(shown), len(with_ev)), fontsize=9, loc="left")
    for k, c in enumerate(chroms):
        if k % 2:
            for a in (ax_c, ax):
                a.axvspan(off[c], off[c] + G.length[c] / 1e6, color=PAL["band"], lw=0, zorder=0)
    for r, s in enumerate(shown):
        y = len(shown) - r
        for e in events:
            if e["sample"] != s or e["chrom"] not in off:
                continue
            x0, x1 = off[e["chrom"]] + int(e["start"]) / 1e6, off[e["chrom"]] + int(e["end"]) / 1e6
            ax.plot([x0, max(x1, x0 + 1.5)], [y, y], color=TYPE_COLOUR.get(e["type"], PAL["lrr"]), lw=3.5 if e["constitutional"] == "1" else 1.6, solid_capstyle="butt")
    ax.set_yticks(range(1, len(shown) + 1))
    ax.set_yticklabels(shown[::-1], fontsize=5)
    ax.set_ylim(0.3, len(shown) + 0.7)
    ax.set_xlim(0, total)
    ax.set_xticks(centres)
    ax.set_xticklabels([c[3:] for c in chroms], fontsize=6)
    ax.set_xlabel("chromosome")
    return save(fig, os.path.join(fdir, "landscape"), "The landscape of events",
                "Top: the number of events per chromosome, stacked by type. Below: one row per sample with an event, each event a bar from its start to its end, coloured by type, thick when constitutional and thin when mosaic.",
                [("vermillion", "gain"), ("blue", "loss"), ("reddish purple", "copy-neutral LOH"), ("orange", "run of homozygosity"), ("thick bar", "constitutional"), ("thin bar", "mosaic")])


def fig_patterns(out):
    """report/figures/patterns: what each event kind looks like on an array - idealised LRR and BAF tracks with array noise
    (LRR SD 0.15, BAF SD 0.04) for a constitutional trisomy, a mosaic gain, a mosaic loss, a mosaic copy-neutral LOH, a run of
    homozygosity, 47,XXY and a mosaic loss of Y. No data: every point is simulated from the expected values."""
    from .baf import d_from_f
    fdir = os.path.join(out, "report", "figures")
    os.makedirs(fdir, exist_ok=True)
    rng = np.random.default_rng(7)
    n = 1200
    x = np.linspace(0, 100, n)
    inside = (x > 35) & (x < 65)
    levels = dict(one=-0.55, three=0.35)
    cards = [("constitutional trisomy", "gain", 1.0), ("mosaic gain, 30% of cells", "gain", 0.3), ("mosaic loss, 30% of cells", "loss", 0.3),
             ("mosaic copy-neutral LOH, 40% of cells", "loh", 0.4), ("run of homozygosity (constitutional)", "roh", 1.0), ("47,XXY (whole X; two X's with a Y)", "xxy", 1.0),
             ("mosaic loss of Y, 50% of cells", "loy", 0.5)]
    fig, axes = plt.subplots(2, len(cards), figsize=(275 * MM, 70 * MM), sharex=True, gridspec_kw=dict(hspace=0.15, wspace=0.25))
    for k, (title, kind, f) in enumerate(cards):
        ax_l, ax_b = axes[0, k], axes[1, k]
        lrr = rng.normal(0, 0.15, n)
        het = rng.random(n) < 0.3
        baf = np.where(het, 0.5, rng.choice([0.0, 1.0], n)) + rng.normal(0, 0.04, n)
        if kind == "gain":
            lrr[inside] += levels["three"] * (f if f < 1 else 1.0)
            d = d_from_f(f, "gain") if f < 1 else 1 / 6
            baf[inside & het] = 0.5 + d * rng.choice([-1, 1], (inside & het).sum()) + rng.normal(0, 0.04, (inside & het).sum())
        elif kind == "loss":
            lrr[inside] += levels["one"] * f
            d = d_from_f(f, "loss")
            baf[inside & het] = 0.5 + d * rng.choice([-1, 1], (inside & het).sum()) + rng.normal(0, 0.04, (inside & het).sum())
        elif kind == "loh":
            d = d_from_f(f, "loh")
            baf[inside & het] = 0.5 + d * rng.choice([-1, 1], (inside & het).sum()) + rng.normal(0, 0.04, (inside & het).sum())
        elif kind == "roh":
            baf[inside & het] = rng.choice([0.0, 1.0], (inside & het).sum()) + rng.normal(0, 0.04, (inside & het).sum())
            het = het & ~inside
        elif kind == "xxy":
            lrr = rng.normal(0, 0.15, n)                                         # two X's: LRR at the two-copy level, heterozygous calls present
            ax_l.text(50, 0.95, "X: two copies, heterozygous\nY: present (LRR near 0)", ha="center", va="top", fontsize=6)
        elif kind == "loy":
            lrr = rng.normal(-0.55 * 1.0, 0.25, n)                              # the Y at half the male level on the compressed scale
            het = np.zeros(n, bool)
            baf = rng.choice([0.0, 1.0], n) + rng.normal(0, 0.04, n)
            ax_l.text(50, 0.95, "Y level halfway between\nthe male and the absent mode", ha="center", va="top", fontsize=6)
        ax_l.scatter(x, lrr, s=1, c=PAL["lrr"], lw=0, rasterized=True)
        for v, col in ((0, PAL["ref"]), (levels["one"], PAL["loss"]), (levels["three"], PAL["gain"])):
            ax_l.axhline(v, color=col, lw=0.5, ls=":")
        ax_b.scatter(x[~het], np.clip(baf[~het], -0.02, 1.02), s=0.8, c=PAL["baf"], lw=0, rasterized=True)
        ax_b.scatter(x[het], np.clip(baf[het], -0.02, 1.02), s=1.5, c=PAL["het"], lw=0, rasterized=True)
        ax_b.axhline(0.5, color=PAL["ref"], lw=0.5, ls=":")
        if kind in ("gain", "loss", "loh", "roh"):
            for a in (ax_l, ax_b):
                a.axvspan(35, 65, color=TYPE_COLOUR[kind], alpha=0.1, lw=0)
        ax_l.set_ylim(-1.2, 1.0)
        ax_b.set_ylim(-0.02, 1.02)
        ax_l.set_title(title, fontsize=7)
        if k == 0:
            ax_l.set_ylabel("LRR")
            ax_b.set_ylabel("BAF")
        ax_b.set_xlabel("position (shaded: the event)", fontsize=6)
        ax_l.tick_params(labelsize=6), ax_b.tick_params(labelsize=6)
    return save(fig, os.path.join(fdir, "patterns"), "Pattern cards: what each event looks like on an array",
                "Idealised LRR (dotted: two, one and three copies) and BAF (green at heterozygous calls) for each kind of event, simulated with typical array noise. A gain lifts the LRR and splits the heterozygous band towards 1/3 and 2/3 (fully for a trisomy, less for a mosaic); a loss lowers the LRR and splits the band outward; a copy-neutral LOH leaves the LRR flat and splits the band; a run of homozygosity has no heterozygous calls and a flat LRR; 47,XXY shows an X at two copies with heterozygous calls and a Y present; a mosaic loss of Y sits between the male and the absent Y level.",
                [("dark grey", "LRR"), ("green", "BAF at heterozygous calls"), ("light grey", "BAF at homozygous calls")])
