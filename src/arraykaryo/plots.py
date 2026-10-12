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
