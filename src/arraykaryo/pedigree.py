"""arraykaryo pedigree: relatedness between every pair of samples from the genotypes alone, the pedigree it implies, and the
comparison with any pedigree given - which is checked against the genotypes and never trusted over them.

Probes: autosomal SNV probes passing QC with a cohort call rate of 0.98 or more and a minor allele frequency of 0.05 or
more, thinned to one probe per 50 kb (about 50,000 on a 2.4-million-probe array; a distance thinning stands in for LD pruning
at this density). Estimator: KING-robust kinship (Manichaikul 2010, the between-family estimator PLINK 2's --make-king
implements, equation 11 of the paper): phi = (2 N_AaAa - 4 N_AAaa - N_Aa(i) - N_Aa(j)) / (4 min(N_Aa(i), N_Aa(j))) + 1/2,
where N_AaAa is the number of probes at which both are heterozygous, N_AAaa the number at which they are opposite
homozygotes (IBS0) and N_Aa(i) the number of heterozygous calls of sample i, each counted over the probes called in both
(tested against PLINK 2 --make-king-table on the mock: kinship within 0.0001, IBS0 exact; the simplified form
(N_AaAa - 2 N_AAaa) / (2 min N_Aa) differs by up to 0.03 when the two samples' heterozygote counts differ, as a child with a
uniparental disomy does). Every pair at once by matrix products of the 0/1 indicator
matrices (heterozygous, homozygous-reference, homozygous-alternate, called), accumulated over probe blocks. IBS0 (the share
of probes at which the pair is opposite homozygous) separates parent-offspring (near zero: the genotyping-error floor)
from full siblings (a quarter of the unrelated expectation, sum 2 p^2 q^2 / n).

Classes by kinship (KING's boundaries, each a parameter): above 0.354 duplicate or monozygotic; 0.177-0.354 first degree,
split by IBS0 into parent-offspring and full sibling at the midpoint between the error floor and the full-sibling
expectation; 0.0884-0.177 second degree; 0.0442-0.0884 third degree. Families from the first-degree graph: a trio is a
sample with two parent-offspring relatives who are unrelated to each other, of opposite sex (section sex), and whose three
genotypes show a Mendelian-error rate at the error floor in that orientation and far above it in any other - the orientation
test resolves who is the child without ages. A parent-offspring pair with no third member is a duo, oriented by the given
pedigree or a sample sheet age when available, else reported as an unoriented pair. The check against a given pedigree says,
per claimed relationship and per discovered one: agrees; a claimed parent unrelated (non-paternity, or a sample swap when
the sex check or the duplicate table says so); an undisclosed first- or second-degree pair; an undisclosed duplicate; a sex
contradiction. Hardy-Weinberg is then recomputed over the unrelated passing samples and the probe table updated."""
import os
import time

import numpy as np

from .cache import Cache
from .qc import probe_masks, read_sample_qc, update_sample_qc
from .sex import read_sex
from .util import NA, log, mark_done, num, read_pedigree, read_sample_sheet, write_json, write_pedigree, write_tsv

PARAMS = dict(call_rate=0.98, maf=0.05, thin_bp=50_000, max_probes=60_000, dup=0.354, first=0.177, second=0.0884, third=0.0442, mie_max=0.02, mie_ratio=3.0, mie_floor_ratio=5.0, mie_floor_min=0.003)
CLASS_NAME = {"dup": "duplicate or monozygotic", "po": "parent-offspring", "fs": "full sibling", "1st": "first degree", "2nd": "second degree", "3rd": "third degree"}
BLOCK = 8000


def select_probes(C, masks, p):
    """(chrom, index array) of the probes used, thinned to one per thin_bp."""
    out = []
    n = 0
    for c in C.chroms:
        if c in ("chrX", "chrY"):
            continue
        pr = C.probes(c)
        m = masks[c]
        ok = m["keep_gt"] & pr.snv & ~pr.intensity_only & (m["call_rate"] >= p["call_rate"]) & np.isfinite(m["maf"]) & (m["maf"] >= p["maf"])
        idx = np.flatnonzero(ok)
        if idx.size == 0:
            continue
        pos = pr.pos[idx]
        keep = [0]
        last = pos[0]
        for k in range(1, idx.size):
            if pos[k] - last >= p["thin_bp"]:
                keep.append(k)
                last = pos[k]
        sel = idx[np.asarray(keep)]
        out.append((c, sel))
        n += sel.size
    if n > p["max_probes"]:
        step = int(np.ceil(n / p["max_probes"]))
        out = [(c, sel[::step]) for c, sel in out]
    return out


def kinship_matrices(C, selected):
    """N_hethet, N_ibs0, N_both (probes called in both), N_het_i (i's heterozygous calls over the probes called in both) as S x S
    arrays, plus the expected IBS0 of unrelated pairs (mean 2 p^2 q^2 over the probes) - accumulated over probe blocks."""
    S = C.n_samples
    hh = np.zeros((S, S)); ibs0 = np.zeros((S, S)); both = np.zeros((S, S)); het_i = np.zeros((S, S))
    e_ibs0, n_used = 0.0, 0
    for chrom, idx in selected:
        gt = C.gt(chrom)
        for a in range(0, idx.size, BLOCK):
            cols = idx[a:a + BLOCK]
            g = np.asarray(gt[:, cols])
            H = (g == 1).astype(np.float32)
            A = (g == 0).astype(np.float32)
            B = (g == 2).astype(np.float32)
            M = (g >= 0).astype(np.float32)
            hh += H @ H.T
            ibs0 += A @ B.T + B @ A.T
            both += M @ M.T
            het_i += H @ M.T
            n = M.sum(axis=0)
            pf = np.where(n > 0, (A.sum(axis=0) * 2 + H.sum(axis=0)) / np.maximum(2 * n, 1), 0.5)
            e_ibs0 += float((2 * pf ** 2 * (1 - pf) ** 2).sum())
            n_used += cols.size
    return hh, ibs0, both, het_i, (e_ibs0 / max(n_used, 1)), n_used


def king_robust(hh, ibs0, het_i):
    """KING-robust kinship for every pair (Manichaikul 2010, eq. 11, as PLINK 2 implements it):
    (2 N_AaAa - 4 N_AAaa - N_Aa(i) - N_Aa(j)) / (4 min(N_Aa(i), N_Aa(j))) + 1/2."""
    m = np.minimum(het_i, het_i.T)
    with np.errstate(all="ignore"):
        k = (2 * hh - 4 * ibs0 - het_i - het_i.T) / np.maximum(4 * m, 1) + 0.5
    return np.where(m > 0, k, np.nan)


def mendelian_errors(kid, dad, mom):
    """Probes, called in all three, at which the child's genotype is impossible given the parents' (a boolean array) and the mask of
    probes called in all three."""
    ok = (kid >= 0) & (dad >= 0) & (mom >= 0)
    k, d, m = kid, dad, mom
    err = ((d == 0) & (m == 0) & (k != 0)) | ((d == 2) & (m == 2) & (k != 2)) | ((d == 0) & (m == 2) & (k != 1)) | ((d == 2) & (m == 0) & (k != 1)) | \
          ((d == 0) & (k == 2)) | ((d == 2) & (k == 0)) | ((m == 0) & (k == 2)) | ((m == 2) & (k == 0))
    return err & ok, ok


def mendelian_error_rate(kid, dad, mom, chrom_of=None):
    """(the genome-wide rate, the rate per chromosome as a dict, the number of probes called in all three). With chrom_of (the
    chromosome of every probe) the per-chromosome rates are returned; the orientation test uses their median, so that one
    chromosome's errors (a uniparental disomy, a deletion in the child) do not hide a true trio."""
    err, ok = mendelian_errors(kid, dad, mom)
    n = int(ok.sum())
    if n == 0:
        return np.nan, {}, 0
    per = {}
    if chrom_of is not None:
        for c in np.unique(chrom_of):
            m = (chrom_of == c) & ok
            if m.sum() >= 50:
                per[str(c)] = float(err[m].sum() / m.sum())
    return float(err.sum() / n), per, n


def robust_rate(per, overall):
    """The median over chromosomes (the overall rate with fewer than three chromosomes)."""
    return float(np.median(list(per.values()))) if len(per) >= 3 else overall


def outlier_chroms(per, floor):
    """Chromosomes whose Mendelian-error rate is far above the others: a uniparental disomy, a deletion, or a chromosome the child did
    not inherit as expected."""
    if len(per) < 3:
        return []
    med = np.median(list(per.values()))
    return sorted((c for c, v in per.items() if v > max(5 * med, floor + 0.01, 0.01)), key=lambda c: -per[c])


def genotypes_of(C, selected, idx):
    """The selected probes' genotypes of a few samples (indices) as (len(idx) x n) int8, and the chromosome of every probe."""
    parts, chroms = [], []
    for chrom, cols in selected:
        gt = C.gt(chrom)
        parts.append(np.asarray(gt[np.ix_(idx, cols)]))
        chroms.append(np.full(cols.size, chrom))
    return np.concatenate(parts, axis=1), np.concatenate(chroms)


def run_pedigree(out, cache_dir=None, pedigree="", sample_sheet="", threads=1, params=None, log_fn=log):
    t0 = time.time()
    cache_dir = cache_dir or os.path.join(out, "cache")
    C = Cache(cache_dir)
    S = C.n_samples
    p = dict(PARAMS)
    if params:
        p.update(params)
    pdir = os.path.join(out, "pedigree")
    os.makedirs(pdir, exist_ok=True)
    masks = probe_masks(out, C)
    qc = read_sample_qc(out)
    sexes = read_sex(out)
    passing = np.array([qc[s]["hard_fail"] == "0" for s in C.samples])
    sex_of = {s: (sexes[s]["sex"] if s in sexes and sexes[s]["sex"] != NA else (qc[s].get("reported_sex", "") if qc[s].get("reported_sex", NA) != NA else "")) for s in C.samples}
    selected = select_probes(C, masks, p)
    n_sel = sum(sel.size for _, sel in selected)
    log_fn("%d probes (autosomal SNV, call rate >= %.2f, MAF >= %.2f, one per %d kb) on %d chromosomes" % (n_sel, p["call_rate"], p["maf"], p["thin_bp"] // 1000, len(selected)), "pedigree")
    hh, ibs0, both, het_i, e_ibs0, n_used = kinship_matrices(C, selected)
    kin = king_robust(hh, ibs0, het_i)
    ibs0_rate = np.where(both > 0, ibs0 / np.maximum(both, 1), np.nan)
    np.savez(os.path.join(pdir, "kinship.npz"), kinship=kin.astype(np.float32), ibs0=ibs0_rate.astype(np.float32), n=both.astype(np.int32), samples=np.array(C.samples), e_ibs0=e_ibs0)
    # the error floor for IBS0: the lowest IBS0 among first-degree pairs, else a tenth of the full-sibling expectation
    iu = np.triu_indices(S, 1)
    kv, iv = kin[iu], ibs0_rate[iu]
    first = np.isfinite(kv) & (kv >= p["first"]) & (kv < p["dup"])
    fs_expect = e_ibs0 / 4.0
    floor = float(np.nanmin(iv[first])) if first.any() else fs_expect / 10.0
    po_cut = (floor + fs_expect) / 2.0
    # classes
    def classify(k, i0):
        if not np.isfinite(k):
            return ""
        if k >= p["dup"]:
            return "dup"
        if k >= p["first"]:
            return "po" if (np.isfinite(i0) and i0 < po_cut) else "fs"
        if k >= p["second"]:
            return "2nd"
        if k >= p["third"]:
            return "3rd"
        return ""
    rel_rows = []
    cls = {}
    for a, b, k, i0 in zip(iu[0], iu[1], kv, iv):
        c = classify(k, i0)
        if c:
            cls[(int(a), int(b))] = c
            note = ""
            if not (passing[a] and passing[b]):
                note = "a sample fails QC"
            rel_rows.append(dict(sample_a=C.samples[a], sample_b=C.samples[b], n_probes=int(both[a, b]), kinship=k, ibs0=i0, hethet=int(hh[a, b]), **{"class": CLASS_NAME[c]}, note=note or NA))
    rel_rows.sort(key=lambda r: -r["kinship"])
    write_tsv(os.path.join(pdir, "relatedness.tsv"), ["sample_a", "sample_b", "n_probes", "kinship", "ibs0", "hethet", "class", "note"], rel_rows, nd=5)
    # families from the relatedness graph (passing samples only): parent-offspring pairs first, but full-sibling- and second-degree-
    # classed pairs too, since one chromosome's missing parental copy (an isodisomy, a deletion) lifts a parent-offspring pair's IBS0
    # and lowers its kinship - the Mendelian-error orientation test decides, and the pair is reclassified when it passes
    po = {}
    for (a, b), c in cls.items():
        if c in ("po", "fs", "2nd") and passing[a] and passing[b]:
            po.setdefault(a, set()).add(b)
            po.setdefault(b, set()).add(a)
    related = lambda a, b: cls.get((min(a, b), max(a, b)), "") != ""  # noqa: E731
    families, used_in_trio, fam_rows = [], set(), []
    mie_floor = None
    reclassified = {}
    for c in sorted(po):
        rels = sorted(po[c])
        for i in range(len(rels)):
            for j in range(i + 1, len(rels)):
                a, b = rels[i], rels[j]
                if related(a, b):
                    continue
                sa, sb = sex_of[C.samples[a]], sex_of[C.samples[b]]
                if sa and sb and sa == sb:
                    continue
                g, chrom_of = genotypes_of(C, selected, [c, a, b])
                orient, per_c = {}, {}
                for key, (x, y, z) in (("c", (0, 1, 2)), ("a", (1, 0, 2)), ("b", (2, 0, 1))):
                    overall, per, _ = mendelian_error_rate(g[x], g[y], g[z], chrom_of)
                    orient[key] = robust_rate(per, overall)
                    per_c[key] = (overall, per)
                best = min(orient, key=lambda k: orient[k])
                others = [v for k, v in orient.items() if k != best]
                if best != "c" or not (orient["c"] <= p["mie_max"] and min(others) >= p["mie_ratio"] * max(orient["c"], 1e-4)):
                    continue
                if len(outlier_chroms(per_c["c"][1], orient["c"])) > 2:
                    continue                                              # errors standing out on several chromosomes: not a trio (a grandparent, an uncle)
                dad, mom = (a, b) if (sa == "M" or sb == "F") else (b, a)
                if not sa and not sb:
                    dad, mom = a, b
                key = (c, dad, mom)
                if key in [(f["kid"], f["dad"], f["mom"]) for f in families]:
                    continue
                notes = [] if (sa and sb) else ["a parent's sex unknown: father and mother assigned by what is known"]
                odd = outlier_chroms(per_c["c"][1], orient["c"])
                if odd:
                    notes.append("Mendelian errors concentrated on %s (%s against %.4f elsewhere): a uniparental disomy, a deletion in the child, or a chromosome not inherited as expected" % (
                        ", ".join(odd), ", ".join("%.4f" % per_c["c"][1][x] for x in odd), orient["c"]))
                for x in (a, b):
                    pair = (min(c, x), max(c, x))
                    if cls.get(pair) in ("fs", "2nd"):
                        reclassified[pair] = "%s by kinship and IBS0 (kinship %.3f, IBS0 %.4f); parent-offspring by the Mendelian-error orientation test" % (CLASS_NAME[cls[pair]], kin[pair], ibs0_rate[pair])
                        cls[pair] = "po"
                families.append(dict(kid=c, dad=dad, mom=mom, kind="trio", mie=orient["c"], mie_overall=per_c["c"][0], mie_other=min(others), note="; ".join(notes)))
                used_in_trio |= {c, dad, mom}
                mie_floor = orient["c"] if mie_floor is None else min(mie_floor, orient["c"])
    # the cohort's Mendelian-error floor from the trios accepted so far: a trio must sit near it, not merely under the absolute cap
    # (a grandparent or an uncle in a parent's place gives 0.01-0.05 where true trios give 0.0003 on a good array)
    if families:
        floor = float(np.median([f["mie"] for f in families]))
        cut = max(p["mie_floor_ratio"] * floor, p["mie_floor_min"])
        dropped = [f for f in families if f["mie"] > cut]
        for f in dropped:
            for x in (f["dad"], f["mom"]):
                pair = (min(f["kid"], x), max(f["kid"], x))
                if pair in reclassified:
                    cls[pair] = "2nd" if kin[pair] < p["first"] else "fs"
                    del reclassified[pair]
        families = [f for f in families if f["mie"] <= cut]
        mie_floor = floor
        if dropped and log_fn:
            log_fn("%d candidate trio(s) dropped: Mendelian-error rate above %.4f (%.0f x the cohort floor %.4f)" % (len(dropped), cut, p["mie_floor_ratio"], floor), "pedigree")
    for r in rel_rows:
        pair = (index_of(C, r["sample_a"]), index_of(C, r["sample_b"]))
        if pair in reclassified:
            r["class"] = CLASS_NAME["po"]
            r["note"] = reclassified[pair]
    if reclassified:
        write_tsv(os.path.join(pdir, "relatedness.tsv"), ["sample_a", "sample_b", "n_probes", "kinship", "ibs0", "hethet", "class", "note"], rel_rows, nd=5)
    # the given pedigree and the sample sheet (ages) for duos
    given_rows, given_sex, _ = read_pedigree(pedigree) if pedigree else ([], {}, {})
    from .util import read_relative_columns
    disclosed = read_relative_columns(pedigree) if pedigree else {}
    given_parent_of = {}
    for r in given_rows:
        for par in (r["dad"], r["mom"]):
            if par:
                given_parent_of[(r["kid"], par)] = True
    ages = {}
    if sample_sheet:
        sheet, cols = read_sample_sheet(sample_sheet)
        age_col = next((c for c in cols if c.strip().lower() in ("age", "age_years", "years")), None)
        if age_col:
            ages = {s: num(r.get(age_col, "")) for s, r in sheet.items()}
    trio_pairs = set()
    for f in families:
        trio_pairs |= {(min(f["kid"], f["dad"]), max(f["kid"], f["dad"])), (min(f["kid"], f["mom"]), max(f["kid"], f["mom"]))}
    for (a, b), c in sorted(cls.items()):
        if c != "po" or (a, b) in trio_pairs or not (passing[a] and passing[b]):
            continue
        na, nb = C.samples[a], C.samples[b]
        kid, par, how = None, None, ""
        if (na, nb) in given_parent_of:
            kid, par, how = a, b, "oriented by the given pedigree"
        elif (nb, na) in given_parent_of:
            kid, par, how = b, a, "oriented by the given pedigree"
        elif np.isfinite(ages.get(na, np.nan)) and np.isfinite(ages.get(nb, np.nan)) and ages[na] != ages[nb]:
            kid, par = (a, b) if ages[na] < ages[nb] else (b, a)
            how = "oriented by the sample sheet ages"
        if kid is None:
            families.append(dict(kid=a, dad=None, mom=None, kind="parent-offspring pair, unoriented", other=b, mie=np.nan, mie_other=np.nan, note="no pedigree or age says who is the parent"))
        else:
            ps = sex_of[C.samples[par]]
            families.append(dict(kid=kid, dad=par if ps == "M" else None, mom=par if ps == "F" else None, kind="duo", other=par if not ps else None, mie=np.nan, mie_other=np.nan,
                                 note=how + ("" if ps else "; the parent's sex is unknown")))
    for (a, b), c in sorted(cls.items()):
        if c == "fs" and passing[a] and passing[b]:
            families.append(dict(kid=a, dad=None, mom=None, kind="full siblings", other=b, mie=np.nan, mie_other=np.nan, note=""))
        if c == "dup":
            families.append(dict(kid=a, dad=None, mom=None, kind="duplicate or monozygotic twins", other=b, mie=np.nan, mie_other=np.nan, note="" if (passing[a] and passing[b]) else "a sample fails QC"))
    # families.tsv and the inferred pedigree
    for k, f in enumerate(families):
        kid = C.samples[f["kid"]]
        members = [kid] + [C.samples[x] for x in (f.get("dad"), f.get("mom"), f.get("other")) if x is not None]
        fam_rows.append(dict(family="F%d" % (k + 1), kind=f["kind"], members=",".join(members), kid=kid if f["kind"] in ("trio", "duo") else NA,
                             dad=C.samples[f["dad"]] if f.get("dad") is not None else NA, mom=C.samples[f["mom"]] if f.get("mom") is not None else NA,
                             other=C.samples[f["other"]] if f.get("other") is not None else NA, kid_sex=sex_of[kid] or NA,
                             mendelian_error_rate=f["mie"], mendelian_error_rate_genome=f.get("mie_overall", np.nan), mendelian_error_rate_other_orientations=f["mie_other"], note=f["note"] or NA))
    write_tsv(os.path.join(pdir, "families.tsv"), ["family", "kind", "members", "kid", "dad", "mom", "other", "kid_sex", "mendelian_error_rate", "mendelian_error_rate_genome",
                                                   "mendelian_error_rate_other_orientations", "note"], fam_rows, nd=5)
    inferred = []
    for f in families:
        if f["kind"] in ("trio", "duo"):
            kid = C.samples[f["kid"]]
            dad = C.samples[f["dad"]] if f.get("dad") is not None else ""
            mom = C.samples[f["mom"]] if f.get("mom") is not None else ""
            inferred.append(dict(kid=kid, dad=dad, mom=mom, kid_sex=sex_of[kid], dad_sex=sex_of.get(dad, "") if dad else "", mom_sex=sex_of.get(mom, "") if mom else ""))
    write_pedigree(os.path.join(pdir, "pedigree.inferred.tsv"), inferred)
    # the check against the given pedigree
    check_rows = []
    flags = {}
    index = {s: i for i, s in enumerate(C.samples)}
    if pedigree:
        for r in given_rows:
            kid = r["kid"]
            if kid not in index:
                check_rows.append(dict(kid=kid, relative=NA, claimed="child of %s and %s" % (r["dad"] or "?", r["mom"] or "?"), found=NA, status="the child is not in the export", note=NA))
                continue
            for par, role in ((r["dad"], "father"), (r["mom"], "mother")):
                if not par:
                    continue
                if par not in index:
                    check_rows.append(dict(kid=kid, relative=par, claimed=role, found=NA, status="the %s is not in the export" % role, note=NA))
                    continue
                a, b = sorted((index[kid], index[par]))
                c = cls.get((a, b), "")
                k = kin[a, b]
                found = CLASS_NAME.get(c, "unrelated")
                if c == "po":
                    status = "agrees"
                    note = ""
                elif c in ("", "3rd", "2nd"):
                    status = "claimed %s unrelated" % role if c == "" else "claimed %s is a %s relative" % (role, found)
                    sc = sexes.get(par, {}).get("check", "")
                    dupe = any(cls.get((min(index[par], x), max(index[par], x)), "") == "dup" for x in range(S))
                    if sc and sc not in ("agrees", "not read", NA) and ("sample swap" in sc):
                        note = "a sample swap: the %s's sex check says %s" % (role, sc)
                    elif dupe:
                        note = "the %s has a duplicate in the cohort" % role
                    else:
                        note = "non-%s, or a swap elsewhere" % ("paternity" if role == "father" else "maternity")
                else:
                    status = "claimed %s is a %s" % (role, found)
                    note = ""
                check_rows.append(dict(kid=kid, relative=par, claimed=role, found=found, kinship=k, status=status, note=note or NA))
            # sex contradictions
            for s, gs in ((kid, r["kid_sex"]), (r["dad"], r["dad_sex"] or "M"), (r["mom"], r["mom_sex"] or "F")):
                if s and gs and s in sexes and sexes[s]["sex"] not in (NA, "") and sexes[s]["sex"] != gs:
                    check_rows.append(dict(kid=kid, relative=s, claimed="%s in the pedigree" % ("male" if gs == "M" else "female"), found=sexes[s]["complement"], status="sex contradiction", note=sexes[s]["check"]))
        claimed_pairs = {(min(index[r["kid"]], index[x]), max(index[r["kid"]], index[x])) for r in given_rows for x in (r["dad"], r["mom"]) if x and r["kid"] in index and x in index}
        for (a, b), c in sorted(cls.items()):
            if (a, b) in claimed_pairs or c in ("3rd",):
                continue
            # siblings of the same given parents are disclosed
            same_parents = False
            pa = [r for r in given_rows if r["kid"] == C.samples[a]]
            pb = [r for r in given_rows if r["kid"] == C.samples[b]]
            if pa and pb and (pa[0]["dad"] or pa[0]["mom"]) and (pa[0]["dad"], pa[0]["mom"]) == (pb[0]["dad"], pb[0]["mom"]):
                same_parents = True
            if c == "fs" and same_parents:
                check_rows.append(dict(kid=C.samples[a], relative=C.samples[b], claimed="siblings", found=CLASS_NAME[c], kinship=kin[a, b], status="agrees", note=NA))
                continue
            listed = disclosed.get(frozenset((C.samples[a], C.samples[b])))
            if listed:
                check_rows.append(dict(kid=C.samples[a], relative=C.samples[b], claimed=listed, found=CLASS_NAME[c], kinship=kin[a, b], status="agrees" if _agree(listed, c) else "listed as %s, found %s" % (listed, CLASS_NAME[c]), note="from the table's relative columns"))
                continue
            status = "undisclosed duplicate" if c == "dup" else "undisclosed %s pair" % CLASS_NAME[c]
            check_rows.append(dict(kid=C.samples[a], relative=C.samples[b], claimed="unrelated", found=CLASS_NAME[c], kinship=kin[a, b], status=status, note=NA))
            if c == "dup":
                flags[C.samples[a]] = {"add_flag": "duplicate"}
                flags[C.samples[b]] = {"add_flag": "duplicate"}
            elif c in ("po", "fs", "2nd"):
                flags.setdefault(C.samples[a], {})["add_flag"] = "unexpected_relative"
                flags.setdefault(C.samples[b], {})["add_flag"] = "unexpected_relative"
        write_tsv(os.path.join(pdir, "check.tsv"), ["kid", "relative", "claimed", "found", "kinship", "status", "note"], check_rows, nd=4)
    else:
        for (a, b), c in cls.items():
            if c == "dup":
                flags[C.samples[a]] = {"add_flag": "duplicate"}
                flags[C.samples[b]] = {"add_flag": "duplicate"}
    update_sample_qc(out, flags)
    # Hardy-Weinberg over the unrelated passing samples
    unrelated = passing.copy()
    for (a, b), c in sorted(cls.items(), key=lambda kv: -kin[kv[0]]):
        if unrelated[a] and unrelated[b]:
            unrelated[b] = False
    from .qc import recompute_hwe
    n_hwe = recompute_hwe(out, C, unrelated)
    n_class = {}
    for c in cls.values():
        n_class[CLASS_NAME[c]] = n_class.get(CLASS_NAME[c], 0) + 1
    kinds = {}
    for f in families:
        kinds[f["kind"]] = kinds.get(f["kind"], 0) + 1
    summary = dict(samples=S, probes=n_sel, pairs_related=len(rel_rows), classes=n_class, families=kinds, e_ibs0_unrelated=e_ibs0, ibs0_floor=floor, po_fs_cut=po_cut,
                   mie_floor=mie_floor, unrelated_passing=int(unrelated.sum()), hwe_failed_probes=n_hwe, check=_count(check_rows, "status") if pedigree else None, params=p, seconds=time.time() - t0)
    write_json(summary, os.path.join(pdir, "summary.json"))
    log_fn("%d related pairs (%s); families: %s; IBS0 floor %.4f, parent-offspring/full-sibling cut %.4f; %d unrelated passing samples, %d probes fail HWE among them" % (
        len(rel_rows), ", ".join("%s %d" % kv for kv in sorted(n_class.items())), ", ".join("%s %d" % kv for kv in sorted(kinds.items())) or "none", floor, po_cut, unrelated.sum(), n_hwe), "pedigree")
    if pedigree:
        log_fn("check against the given pedigree: %s" % ", ".join("%s %d" % kv for kv in sorted(summary["check"].items())), "pedigree")
    try:
        from .plots import fig_kinship
        fig_kinship(out, kin[iu], ibs0_rate[iu], [cls.get((int(a), int(b)), "") for a, b in zip(iu[0], iu[1])], p, po_cut, log_fn)
    except Exception as e:
        log_fn("figure skipped: %s" % e, "pedigree")
    mark_done(out, "pedigree", dict(seconds=time.time() - t0))
    return summary


def _agree(listed, c):
    return (listed == "siblings" and c in ("fs", "po")) or (listed == "second order" and c == "2nd") or (listed == "third order" and c in ("3rd", "2nd"))


def index_of(C, name):
    return C.sample_index(name)


def _count(rows, key):
    out = {}
    for r in rows:
        out[r[key]] = out.get(r[key], 0) + 1
    return out


def read_families(out):
    from .util import read_tsv
    p = os.path.join(out, "pedigree", "families.tsv")
    return read_tsv(p) if os.path.exists(p) else []
