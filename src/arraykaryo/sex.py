"""arraykaryo sex: the sex-chromosome complement of every sample from the X outside the pseudoautosomal regions and the Y,
constitutional or mosaic, against the reported sex.

Per sample: the heterozygous rate on the X relative to the sample's autosomal rate; the median corrected LRR on the X (the
reference built from the samples without a Y, so one X reads as one copy); the median GC-corrected LRR on the Y over the
probes that are empirically Y-specific (those whose level is near the male cluster in the samples with a Y and far below it
in the samples without: the X-transposed and PAR-like probes of the Y are excluded by this test, not by a table); the Y call
rate. The cohort's two X modes (one copy, two copies) and two Y modes (present, absent) are found by mixture fitting on the
passing samples and give the copy scale - arrays compress LRR, a one-copy X sits near -0.5, and the scale is the cohort's
own: copies of X = 1 + (LRR_X - m1) / (m2 - m1); copies of Y = 2^((LRR_Y - mY) / (m2 - m1)), the Y level read on the same
compression (Forsberg 2014's mLRR-Y reading, calibrated on the cohort). The complement the Y implies is the anchor (XY with
a Y, XX without; the reported sex plays no part, as in triokaryo's sexchrom): an X deviating from the anchor's count by a
cell fraction f of 0.10 or more is a mosaic, by 0.9 or more a constitutional aneuploidy (47,XXY, 45,X, 47,XXX); a Y below one
copy in a man is a mosaic loss of Y; a Y above is 47,XYY; a Y in part of the cells of an XX genome is named as such. The BAF
at X heterozygous calls gives the second estimate of f (the band deviation d: loss f = 4d/(1+2d), gain f = 4d/(1-2d); an
extra X in part of the cells of an XY sample shows heterozygous-looking calls with d = (1-f)/(2(1+f))). The check against
the reported sex says what a disagreement means (a swap, an XX male, an XY female, or an aneuploidy a heterozygosity-called
sex mislabels: peddy's reads 47,XXY as female and 45,X as male)."""
import os
import time

import numpy as np

from .baf import f_gain, f_loss, fit_d
from .cache import Cache
from .qc import Correction, probe_masks, read_sample_qc, update_sample_qc
from .util import NA, finite, log, mark_done, num, robust_sd, write_json, write_tsv

MIN_F = 0.10           # a deviation in copies this large is a mosaic - or 3 robust SD of the cluster's levels when that is larger (the X and Y
                       # levels of men spread by about 0.07 copies on the public example; a fixed floor would read that spread as mosaics)
CLUSTER_Z = 3.0
D_X_MIN = 0.03         # the band deviation at X heterozygous calls that confirms an X mosaic in an XX sample (the null's 95th percentile on the example)
LLR_X_MIN = 30.0
CONST_F = 0.90         # at or above: constitutional
Y_PRESENT = 0.25       # Y copies at or above: a Y, in all or part of the cells
ONE_COPY_DEFAULT = -0.55   # the one-copy X level relative to two copies when the cohort has no male (or no female) to fit it
KARYOTYPE_NAME = {"XXY": "47,XXY", "XYY": "47,XYY", "XXX": "47,XXX", "X": "45,X", "XXYY": "48,XXYY", "XXXY": "48,XXXY", "XXXX": "48,XXXX", "XY": "46,XY", "XX": "46,XX"}


def two_modes(v, init_lo, init_hi, iters=50):
    """A two-component Gaussian mixture on one dimension by EM from two starting means; returns (m_lo, m_hi, sd_lo, sd_hi, w_lo)."""
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 4 or not (np.isfinite(init_lo) and np.isfinite(init_hi)):
        return np.nan, np.nan, np.nan, np.nan, np.nan
    m = np.array([init_lo, init_hi], dtype=float)
    s = np.array([0.1, 0.1])
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        ll = -0.5 * ((v[:, None] - m[None, :]) / s[None, :]) ** 2 - np.log(s[None, :]) + np.log(w[None, :])
        ll -= ll.max(axis=1, keepdims=True)
        r = np.exp(ll)
        r /= r.sum(axis=1, keepdims=True)
        n = r.sum(axis=0) + 1e-9
        m = (r * v[:, None]).sum(axis=0) / n
        s = np.sqrt((r * (v[:, None] - m[None, :]) ** 2).sum(axis=0) / n) + 0.01
        w = n / n.sum()
    if m[0] > m[1]:
        m, s, w = m[::-1], s[::-1], w[::-1]
    return float(m[0]), float(m[1]), float(s[0]), float(s[1]), float(w[0])


def y_specific_probes(C, keep, par, has_y, no_y, min_group=5):
    """Y probes whose level sits near the male cluster in the samples with a Y and far below in those without."""
    base = keep & ~par
    if has_y.sum() < min_group or no_y.sum() < min_group:
        return base, "all non-PAR Y probes (fewer than %d samples in one of the Y groups: no empirical selection)" % min_group
    lrr = C.lrr("chrY")
    with np.errstate(all="ignore"):
        m_yes = np.nanmedian(np.asarray(lrr[has_y], dtype=np.float32), axis=0)
        m_no = np.nanmedian(np.asarray(lrr[no_y], dtype=np.float32), axis=0)
    spec = base & np.isfinite(m_yes) & np.isfinite(m_no) & (m_yes > -0.6) & (m_no < -2.0) & (m_yes - m_no > 2.0)
    if spec.sum() < 20:
        return base, "all non-PAR Y probes (fewer than 20 passed the empirical test)"
    return spec, "%d of %d non-PAR Y probes are Y-specific (near the male level in the samples with a Y, below -2 in those without)" % (spec.sum(), base.sum())


def anchor(x_copies, y_copies, x_het_rel=np.nan):
    """XY ('M') with a Y, XX ('F') without; a Y in part of the cells goes with the nearer complement; one X without
    heterozygosity and little or no Y is read against XY too (a man whose cells lost the Y, or a 45,X: both are one X, and
    the Y's share of cells says which)."""
    if not np.isfinite(y_copies):
        return ""
    if y_copies >= 0.5:
        return "M"
    if np.isfinite(x_copies) and x_copies < 1.5 and np.isfinite(x_het_rel) and x_het_rel < 0.15:
        return "M"
    if y_copies < Y_PRESENT:
        return "F"
    if not np.isfinite(x_copies):
        return "M"
    return "M" if abs(x_copies - 1.0) + (1.0 - y_copies) <= abs(x_copies - 2.0) + y_copies else "F"


def complement_of(x_copies, y_copies):
    nx = int(np.clip(np.rint(x_copies), 1, 4)) if np.isfinite(x_copies) else 0
    ny = int(np.clip(np.rint(y_copies), 0, 2)) if np.isfinite(y_copies) else 0
    if nx + ny == 0:
        return ""
    return "X" * nx + "Y" * ny


def sex_check(reported, comp):
    """The complement against the reported sex, in triokaryo's words."""
    if not reported:
        return ""
    if not comp:
        return "not read"
    expected = "XY" if reported == "M" else "XX"
    if comp == expected:
        return "agrees"
    who = "male" if reported == "M" else "female"
    if reported == "M" and comp == "XX":
        why = "an XX male, or a sample swap"
    elif reported == "F" and comp == "XY":
        why = "a sample swap, or an XY female"
    elif reported == "F" and comp == "XXY":
        why = "47,XXY; a sex called from X heterozygosity reads two X's as female"
    elif reported == "M" and comp == "X":
        why = "45,X; a sex called from X heterozygosity reads one X as male"
    elif reported == "M" and comp == "XXY":
        why = "47,XXY"
    elif reported == "F" and comp == "X":
        why = "45,X"
    else:
        why = KARYOTYPE_NAME.get(comp, comp)
    return "%s in a reported %s (%s)" % (comp, who, why)


def run_sex(out, cache_dir=None, threads=1, params=None, log_fn=log):
    t0 = time.time()
    cache_dir = cache_dir or os.path.join(out, "cache")
    C = Cache(cache_dir)
    S = C.n_samples
    corr = Correction(out)
    masks = probe_masks(out, C)
    qc = read_sample_qc(out)
    sdir = os.path.join(out, "sex")
    os.makedirs(sdir, exist_ok=True)
    p = dict(min_f=MIN_F, const_f=CONST_F, y_present=Y_PRESENT)
    if params:
        p.update(params)
    passing = np.array([qc[s]["hard_fail"] == "0" for s in C.samples])
    reported = np.array([qc[s]["reported_sex"] if qc[s]["reported_sex"] != NA else "" for s in C.samples])
    auto_het = np.array([num(qc[s]["het_rate"]) for s in C.samples])
    baf_sd = np.array([num(qc[s].get("baf_sd_robust", qc[s]["baf_sd"])) for s in C.samples])
    # the X
    have_x, have_y = "chrX" in C.chroms, "chrY" in C.chroms
    x_lrr = np.full(S, np.nan)
    x_het = np.full(S, np.nan)
    n_x = n_xhet = np.zeros(S, int)
    d_x = np.full(S, np.nan)
    x_het_bafs = [None] * S
    if have_x:
        prx = C.probes("chrX")
        mx = masks["chrX"]
        use_l = mx["keep"] & ~prx.par
        use_g = mx["keep_gt"] & ~prx.par
        Lx = corr.apply(C, "chrX")
        with np.errstate(all="ignore"):
            x_lrr = np.nanmedian(np.where(use_l[None, :], Lx, np.nan), axis=1)
        gx = np.asarray(C.gt("chrX"))
        bx = np.asarray(C.baf("chrX"), dtype=np.float32)
        called = (gx >= 0) & use_g[None, :]
        het = (gx == 1) & use_g[None, :]
        n_x = called.sum(axis=1)
        n_xhet = het.sum(axis=1)
        x_het = n_xhet / np.maximum(n_x, 1)
        for i in range(S):
            if n_xhet[i] >= 20:
                x_het_bafs[i] = bx[i][het[i]]
    x_het_rel = x_het / np.maximum(auto_het, 1e-6)
    # the Y
    y_lrr = np.full(S, np.nan)
    y_call = np.full(S, np.nan)
    n_y = np.zeros(S, int)
    y_note = "no Y in the export"
    if have_y:
        pry = C.probes("chrY")
        my = masks["chrY"]
        Ly = corr.apply(C, "chrY", reference=False)
        base = my["keep"] & ~pry.par
        with np.errstate(all="ignore"):
            y0 = np.nanmedian(np.where(base[None, :], Ly, np.nan), axis=1)
        mid = _midpoint(y0[passing])
        has_y0 = np.isfinite(y0) & (y0 >= mid) if np.isfinite(mid) else np.isfinite(y0) & (y0 > -1.5)
        no_y0 = np.isfinite(y0) & ~has_y0
        spec, y_note = y_specific_probes(C, my["keep"], pry.par, has_y0 & passing, no_y0 & passing)
        with np.errstate(all="ignore"):
            y_lrr = np.nanmedian(np.where(spec[None, :], Ly, np.nan), axis=1)
        gy = np.asarray(C.gt("chrY"))
        n_y = np.full(S, int(spec.sum()))
        y_call = ((gy >= 0) & spec[None, :] & my["keep_gt"][None, :]).sum(axis=1) / max(int((spec & my["keep_gt"]).sum()), 1)
    # the modes, on the passing samples
    pv = passing & np.isfinite(x_lrr)
    if have_y and np.isfinite(y_lrr).any():
        mid = _midpoint(y_lrr[passing & np.isfinite(y_lrr)])
        yes = np.isfinite(y_lrr) & (y_lrr >= mid) if np.isfinite(mid) else np.zeros(S, bool)
        no = np.isfinite(y_lrr) & ~yes
    else:
        yes, no = np.zeros(S, bool), np.ones(S, bool)
    m1_init = float(np.nanmedian(x_lrr[pv & yes])) if (pv & yes).sum() >= 3 else np.nan
    m2_init = float(np.nanmedian(x_lrr[pv & no])) if (pv & no).sum() >= 3 else np.nan
    note_modes = ""
    if np.isfinite(m1_init) and np.isfinite(m2_init):
        m1, m2, s1, s2, w1 = two_modes(x_lrr[pv], m1_init, m2_init)
        if not (np.isfinite(m1) and m2 - m1 > 0.15):
            m1, m2 = m1_init, m2_init
    elif np.isfinite(m2_init):
        m2, m1 = m2_init, m2_init + ONE_COPY_DEFAULT
        note_modes = "no sample with a Y: the one-copy X level is the literature default (%.2f below two copies)" % -ONE_COPY_DEFAULT
    elif np.isfinite(m1_init):
        m1, m2 = m1_init, m1_init - ONE_COPY_DEFAULT
        note_modes = "no sample without a Y: the two-copy X level is the literature default (%.2f above one copy)" % -ONE_COPY_DEFAULT
    else:
        m1, m2 = ONE_COPY_DEFAULT, 0.0
        note_modes = "no X levels to fit"
    k = m2 - m1
    mY1 = float(np.nanmedian(y_lrr[passing & yes])) if (passing & yes).sum() >= 3 else np.nan
    mY0 = float(np.nanmedian(y_lrr[passing & no])) if (passing & no).sum() >= 3 else np.nan
    x_copies = 1.0 + (x_lrr - m1) / k
    if np.isfinite(mY1):
        y_copies = np.where(np.isfinite(y_lrr), np.minimum(2.0 ** ((y_lrr - mY1) / k), 3.0), np.nan)
    elif have_y and np.isfinite(mY0):
        y_copies = np.where(np.isfinite(y_lrr), np.minimum(2.0 ** ((y_lrr - (mY0 + 4.0)) / k), 3.0), np.nan)       # no sample with a Y: the male level assumed 4 above the absent one
    else:
        y_copies = np.full(S, np.nan)
    # the spread of the levels within each cluster: a mosaic must stand out from it
    anchors = np.array([anchor(float(x_copies[i]), float(y_copies[i]), float(x_het_rel[i])) for i in range(S)])
    def spread(v, m):
        vv = v[m & passing & np.isfinite(v)]
        return robust_sd(vv) if vv.size >= 10 else 0.0
    s_x = {"M": spread(x_copies, anchors == "M"), "F": spread(x_copies, anchors == "F")}
    s_y = {"M": spread(y_copies, anchors == "M"), "F": spread(y_copies, anchors == "F")}
    floor_x = {b: max(p["min_f"], CLUSTER_Z * s_x[b]) for b in ("M", "F")}
    floor_y = {b: max(p["min_f"], CLUSTER_Z * s_y[b]) for b in ("M", "F")}
    # per sample
    rows, counts = [], {}
    updates = {}
    for i, name in enumerate(C.samples):
        xc, yc = float(x_copies[i]), float(y_copies[i])
        base = anchor(xc, yc, float(x_het_rel[i]))
        comp = complement_of(xc, yc)
        exp_x = 1 if base == "M" else 2
        mosaic, f, f_lrr, f_baf, kar = "", np.nan, np.nan, np.nan, ""
        notes = []
        # the BAF at X heterozygous calls
        dd, hx, llr = np.nan, x_het_bafs[i], 0.0
        if hx is not None and hx.size >= 20:
            dd, _, llr = fit_d(hx, sigma=baf_sd[i] if finite(baf_sd[i]) and baf_sd[i] > 0.01 else None)
        d_x[i] = dd
        d_sig = finite(dd) and dd >= D_X_MIN and llr >= LLR_X_MIN
        if base and np.isfinite(xc) and xc < 0.5:
            notes.append("X at %.2f copies: below one copy, not read as a complement" % xc)
            comp = ""
        elif base and np.isfinite(xc):
            dev = xc - exp_x
            kind = "gain" if dev > 0 else "loss"
            f_lrr = float(min(abs(dev), 1.0))
            if base == "F" and d_sig:
                f_baf = float(min(f_loss(dd) if kind == "loss" else f_gain(dd), 1.0))
            elif base == "M" and kind == "gain" and d_sig and x_het_rel[i] > 0.15:
                f_baf = float(np.clip((1 - 2 * dd) / (1 + 2 * dd), 0, 1))
            # a mosaic stands out from its cluster's spread; in an XX sample the BAF split confirms it (or a large level change stands alone);
            # in an XY sample an extra X shows heterozygous-looking calls, and a lost single X is not a viable state
            big = abs(dev) >= max(0.25, 5.0 * s_x[base])
            if base == "F":
                called = abs(dev) >= floor_x["F"] and (finite(f_baf) or big)
            elif kind == "gain":
                called = abs(dev) >= floor_x["M"] and (x_het_rel[i] > 0.15 or big)
            else:
                called = False
                if abs(dev) >= floor_x["M"]:
                    notes.append("X at %.2f copies in an XY sample: a single X lost in part of the cells is not a viable state (an artefact, or a mixture)" % xc)
            if called:
                f = f_baf if finite(f_baf) else f_lrr
                if base == "M" and kind == "gain" and x_het_rel[i] >= 0.7:
                    f = f_lrr                                              # two fully heterozygous X's: the extra X is in every cell; the band formula is for a part
                if finite(f_baf) and finite(f_lrr) and abs(f_baf - f_lrr) > 0.2:
                    notes.append("the LRR and the BAF disagree on the cell fraction (%.2f against %.2f)" % (f_lrr, f_baf))
                comp2 = comp if f >= p["const_f"] else ("X" * (exp_x + (1 if kind == "gain" else -1)) + ("Y" if base == "M" else ""))
                if f >= p["const_f"]:
                    mosaic = "constitutional %s" % KARYOTYPE_NAME.get(comp, comp)
                    kar = KARYOTYPE_NAME.get(comp, "%d,%s" % (44 + len(comp), comp))
                else:
                    mosaic = "mosaic %s of an X in %.0f%% of cells (%s/%s)" % (kind, 100 * f, KARYOTYPE_NAME.get(comp2, comp2), KARYOTYPE_NAME.get(base_comp(base), base_comp(base)))
                    kar = "mos %s[%.2f]/%s" % (KARYOTYPE_NAME.get(comp2, "%d,%s" % (44 + len(comp2), comp2)), f, KARYOTYPE_NAME[base_comp(base)])
                    comp = base_comp(base)
            else:
                f_lrr, f_baf = np.nan, np.nan
                comp = base_comp(base)
                kar = KARYOTYPE_NAME[comp]
        if base == "M" and np.isfinite(yc):
            dev_y = yc - 1.0
            if dev_y <= -floor_y["M"]:
                fy = float(min(-dev_y, 1.0))
                if yc < 0.05:
                    txt = "no Y: 45,X (in a man, the Y lost in every cell)"
                elif fy >= p["const_f"]:
                    txt = "loss of Y in %.0f%% of cells" % (100 * fy)
                else:
                    txt = "mosaic loss of Y in %.0f%% of cells" % (100 * fy)
                mosaic = (mosaic + "; " if mosaic else "") + txt
                f = fy if not np.isfinite(f) else f
                f_lrr = fy if not np.isfinite(f_lrr) else f_lrr
                if fy < p["const_f"]:
                    kar = "mos 45,X[%.2f]/%s" % (fy, kar or "46,XY")
                else:                                                      # the Y gone from (nearly) every cell: the complement is X
                    comp = complement_of(xc, 0.0) if comp in ("XY", "") else comp
                    kar = KARYOTYPE_NAME.get(comp, "%d,%s" % (44 + len(comp), comp))
                notes.append("Y at %.2f copies against 1 (%s)" % (yc, txt))
            elif dev_y >= floor_y["M"]:
                fy = float(min(dev_y, 2.0))
                if fy >= p["const_f"]:
                    txt = "constitutional 47,XYY" if comp.count("Y") >= 2 else "an extra Y"
                else:
                    txt = "mosaic gain of a Y in %.0f%% of cells (46,XY/47,XYY)" % (100 * fy)
                mosaic = (mosaic + "; " if mosaic else "") + txt
                f = fy if not np.isfinite(f) else f
                f_lrr = fy if not np.isfinite(f_lrr) else f_lrr
                if fy < p["const_f"]:
                    kar = "mos 47,XYY[%.2f]/%s" % (fy, kar or "46,XY")
        elif base == "F" and np.isfinite(yc) and yc >= p["y_present"]:
            fy = float(min(yc, 1.0))
            mosaic = (mosaic + "; " if mosaic else "") + "a Y in %.0f%% of cells (46,XX/47,XXY mosaic, or an XX/XY mixture)" % (100 * fy)
            f = fy if not np.isfinite(f) else f
            kar = "mos 47,XXY[%.2f]/%s" % (fy, kar or "46,XX")
        called_sex = "M" if base == "M" else ("F" if base == "F" else "")
        chk = sex_check(reported[i], comp) if comp else ("not read" if reported[i] else "")
        if not passing[i]:
            notes.append("sample fails QC: not to be read")
        if have_x and n_x[i] < 200:
            notes.append("thin: %d X probes called" % n_x[i])
        if have_y and n_y[i] < 20:
            notes.append("thin: %d Y probes" % n_y[i])
        if base == "M" and x_het_rel[i] > 0.3 and not mosaic:
            notes.append("X heterozygosity %.2f of the autosomal rate in an XY sample: an X in part of the cells, or contamination" % x_het_rel[i])
        rows.append(dict(sample=name, reported_sex=reported[i] or NA, qc_pass=int(passing[i]), x_lrr=x_lrr[i], y_lrr=y_lrr[i], x_het_rate=x_het[i], x_het_rel=x_het_rel[i], y_call_rate=y_call[i],
                         n_x_called=int(n_x[i]), n_x_het=int(n_xhet[i]), n_y_probes=int(n_y[i]), x_copies=xc, y_copies=yc, base=base_comp(base) if base else NA, complement=comp or NA,
                         sex=called_sex or NA, karyotype=kar or NA, mosaic=mosaic or NA, f=f, f_lrr=f_lrr, f_baf=f_baf, d_x=dd, check=chk or NA, note="; ".join(notes) or NA))
        key = (comp or "?") + (" mosaic" if mosaic and "constitutional" not in mosaic else "")
        counts[key] = counts.get(key, 0) + 1
        upd = dict(sex=called_sex or NA, sex_check=chk or NA)
        if chk and chk != "agrees" and comp in ("XX", "XY"):
            upd["add_flag"] = "sex_mismatch"
        updates[name] = upd
    cols = ["sample", "reported_sex", "qc_pass", "x_lrr", "y_lrr", "x_het_rate", "x_het_rel", "y_call_rate", "n_x_called", "n_x_het", "n_y_probes", "x_copies", "y_copies", "base", "complement",
            "sex", "karyotype", "mosaic", "f", "f_lrr", "f_baf", "d_x", "check", "note"]
    write_tsv(os.path.join(sdir, "sex.tsv"), cols, rows, nd=4)
    update_sample_qc(out, updates, new_cols=("sex", "sex_check"))
    checks = {}
    for r in rows:
        checks[r["check"]] = checks.get(r["check"], 0) + 1
    summary = dict(samples=S, x_modes=dict(one_copy=m1, two_copies=m2, scale_per_copy=k), y_modes=dict(present=mY1, absent=mY0), y_probes=y_note, note=note_modes, complements=counts, checks=checks,
                   cluster_spread=dict(x_in_xy=s_x["M"], x_in_xx=s_x["F"], y_in_xy=s_y["M"], y_in_xx=s_y["F"]), mosaic_floor=dict(x_in_xy=floor_x["M"], x_in_xx=floor_x["F"], y_in_xy=floor_y["M"]),
                   params=p, seconds=time.time() - t0)
    write_json(summary, os.path.join(sdir, "summary.json"))
    log_fn("X modes %.3f (one copy) and %.3f (two copies), scale %.3f per copy; Y modes %.2f present, %.2f absent; %s" % (m1, m2, k, mY1 if np.isfinite(mY1) else np.nan, mY0 if np.isfinite(mY0) else np.nan, y_note), "sex")
    log_fn("cluster spread (robust SD in copies): X in XY %.3f, X in XX %.3f, Y in XY %.3f; mosaic floors: X in XY %.2f, X in XX %.2f, Y in XY %.2f" % (s_x["M"], s_x["F"], s_y["M"], floor_x["M"], floor_x["F"], floor_y["M"]), "sex")
    log_fn("complements: %s; checks: %s" % (", ".join("%s %d" % kv for kv in sorted(counts.items())), ", ".join("%s %d" % kv for kv in sorted(checks.items()))), "sex")
    try:
        from .plots import fig_sex, fig_sex_tracks
        fig_sex(out, rows, summary, log_fn)
        odd = [r for r in rows if r["mosaic"] != NA or r["complement"] not in ("XX", "XY") or (r["check"] not in ("agrees", NA, "not read"))]
        for r in odd:
            fig_sex_tracks(out, C, corr, masks, r)
        if odd:
            log_fn("X and Y tracks drawn for %d samples not plainly XX or XY: sex/figures/" % len(odd), "sex")
    except Exception as e:
        log_fn("figures skipped: %s" % e, "sex")
    mark_done(out, "sex", dict(seconds=time.time() - t0))
    return summary


def base_comp(base):
    return "XY" if base == "M" else "XX"


def _midpoint(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 4 or v.max() - v.min() < 1.0:
        return np.nan
    c1, c2 = np.percentile(v, 10), np.percentile(v, 90)
    for _ in range(30):
        mid = (c1 + c2) / 2
        a, b = v[v < mid], v[v >= mid]
        if a.size == 0 or b.size == 0:
            return np.nan
        c1, c2 = np.median(a), np.median(b)
    return (c1 + c2) / 2 if c2 - c1 > 1.0 else np.nan


def read_sex(out):
    from .util import read_tsv
    p = os.path.join(out, "sex", "sex.tsv")
    return {r["sample"]: r for r in read_tsv(p)} if os.path.exists(p) else {}
