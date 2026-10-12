"""arraykaryo mock: a simulated array export in the gtc2vcf style, with families, planted events and QC failures, and the
truth every test reads. No real data: every genotype is drawn from a population allele frequency, every intensity from a
copy-number model.

The model. Each probe has a B-allele frequency; founders carry two haplotypes per autosome (one per X in a man, one Y),
children inherit one homologue from each parent at every probe (Mendelian transmission; no linkage, which no step here uses).
A sample's signal at a probe follows from the copies present: in normal cells the constitutional homologues (A copies of
the B allele among C copies), in event cells C + dC copies with A + dB of B; in a cell fraction f the mean B-allele fraction
is ((1-f) A + f (A+dB)) / ((1-f) C + f (C+dC)) and the LRR the array's compressed level of the mean copy number (one copy
near -0.55, three near +0.35, the literature's PennCNV levels being -0.66 and +0.40). Noise: a per-sample LRR SD drawn in
0.10-0.25, a BAF SD at heterozygous calls in 0.03-0.06, a per-sample GC wave (amplitude about 0.03 on the probe's GC and
on its 1-Mb window), a per-probe cohort-wide offset (the array's structure, which the cohort reference removes), a per-batch
per-probe offset, genotype errors (0.2%), no-calls and missing values. Genotypes are called from the observed BAF as GenCall
would: below 0.15 AA, above 0.85 BB, between AB. A female's Y probes are no-calls with LRR near -4; PAR probes sit on the X
with two copies in a man (three in a 47,XXY, one in a 45,X, 2 - f under a mosaic loss of Y)."""
import os
import shutil
import subprocess
import time

import numpy as np
import pysam

from . import __version__
from .genome import load_genome
from .util import log, write_json, write_pedigree, write_tsv

DEFAULT_CONTIGS = "chr1,chr7,chr11,chr12,chr15,chr21,chrX,chrY"
LEVELS = np.array([-3.5, -0.55, 0.0, 0.35, 0.60, 0.80])      # LRR at 0, 1, 2, 3, 4, 5 copies: the array's compression
Y_ABSENT = -4.0
MIN_SAMPLES = 30


def level(c):
    """The LRR of a mean copy number (piecewise linear between the integer levels)."""
    c = np.clip(np.asarray(c, dtype=float), 0, 5)
    return np.interp(c, np.arange(6), LEVELS)


def level_y(c, rng):
    """The Y's LRR against the male cluster: 0.55 log2(copies) for a Y present in part of the cells, the absent level with wide noise below 5%."""
    c = np.asarray(c, dtype=float)
    out = np.where(c >= 0.05, 0.55 * np.log2(np.maximum(c, 0.05)), Y_ABSENT + rng.normal(0, 0.8, c.shape))
    return out


# ---------------------------------------------------------------------------------------------------------------- probes

def make_probes(G, chrom, density_kb, rng, planted):
    """Positions, GC, allele frequency and the per-probe planted failures of one chromosome."""
    L = G.length[chrom]
    cen = G.centromere[chrom]
    first = G.p_end[chrom] + 1 if chrom in ("chr13", "chr14", "chr15", "chr21", "chr22") else 1    # no probes on an acrocentric short arm
    spacing = rng.exponential(density_kb * 1000.0, size=int(L / (density_kb * 1000.0) * 1.3) + 10)
    pos = first + np.cumsum(spacing).astype(np.int64)
    pos = pos[pos <= L - 500]
    if cen[1]:
        pos = pos[(pos < cen[0] - 500_000) | (pos > cen[1] + 500_000)]              # the centromere and its flanks carry no probes
    n = pos.size
    # GC: a slow field along the chromosome plus a local term, so that the 1-Mb window GC carries information
    slow = np.cumsum(rng.normal(0, 1, n))
    slow = (slow - slow.mean()) / (slow.std() + 1e-9)
    gc = np.clip(0.41 + 0.06 * slow + 0.05 * rng.normal(0, 1, n), 0.18, 0.78).astype(np.float32)
    p = np.clip(rng.beta(0.6, 0.6, n), 0.0, 1.0)
    mono = rng.random(n) < 0.05
    p[mono] = 0.0
    kind = np.zeros(n, dtype=np.int8)               # 0 SNV, 1 intensity-only, 2 indel, 3 multiallelic
    kind[rng.random(n) < 0.03] = 1
    gentrain = np.clip(rng.beta(14, 1.5, n), 0, 1).astype(np.float32)
    csep = np.clip(rng.beta(6, 2, n), 0, 1).astype(np.float32)
    fail = {}                                       # probe index -> reason
    if planted and chrom not in ("chrX", "chrY") and n > 2000:
        cand = np.flatnonzero((kind == 0) & (p > 0.1) & (p < 0.9))
        rng.shuffle(cand)
        take = lambda k: [int(i) for i in cand[len(fail):len(fail) + k]]  # noqa: E731
        for i in take(6):
            fail[i] = "non_snv"
            kind[i] = 2 if rng.random() < 0.5 else 3
        for i in take(8):
            fail[i] = "call_rate"
        for i in take(8):
            fail[i] = "lrr_median"
        for i in take(8):
            fail[i] = "baf_shift"
        for i in take(8):
            fail[i] = "gentrain"
            gentrain[i] = rng.uniform(0.2, 0.55)
        for i in take(6):
            fail[i] = "cluster_sep"
            csep[i] = rng.uniform(0.05, 0.25)
        for i in take(6):
            fail[i] = "hwe"
            p[i] = 0.5
    bases = np.array(list("ACGT"))
    ref = bases[rng.integers(0, 4, n)].astype(object)
    alt = np.array([bases[(list("ACGT").index(r) + rng.integers(1, 4)) % 4] for r in ref]).astype(object)
    for i, why in fail.items():
        if why == "non_snv":
            if kind[i] == 2:
                ref[i] = ref[i] + "T"
            else:
                alt[i] = alt[i] + "," + bases[(list("ACGT").index(ref[i][0]) + 2) % 4]
    alt[kind == 1] = "."
    ids = np.array(["mock_%s_%d" % (chrom[3:], x) for x in pos], dtype=object)
    return dict(chrom=chrom, pos=pos, gc=gc, p=p, kind=kind, gentrain=gentrain, csep=csep, ref=ref.astype(object), alt=alt, id=ids, fail=fail, par=G.par_mask(chrom, pos))


# --------------------------------------------------------------------------------------------------------------- persons

class Person:
    def __init__(self, name, sex, kind="single"):
        self.name = name                  # the person; the export column may carry another name (the swap)
        self.sex = sex                    # "M" / "F": the sex the body has (XXY and 45,X carry their own complement)
        self.kind = kind
        self.complement = "XY" if sex == "M" else "XX"
        self.dad = self.mom = None
        self.haps = {}                    # chrom -> (H x P) int8, -1 where a homologue is absent at the probe
        self.events = []                  # dicts: chrom, start, end (probe indices), f, dB (array or None), dC, type, label, class
        self.sex_events = []
        self.lrr_sd = self.baf_sd = self.wave = self.wave2 = self.shift = 0.0
        self.batch = ""
        self.role = ""
        self.copy_of = None               # the duplicate pair: the person whose haplotypes are shared


def founder_haps(pr, chrom, person, rng, par):
    """Constitutional homologues of a founder at every probe of the chromosome."""
    n = pr["pos"].size
    p = pr["p"]
    draw = lambda: (rng.random(n) < p).astype(np.int8)  # noqa: E731
    comp = person.complement
    if chrom == "chrY":
        return y_haps(n, draw, comp, par)
    if chrom == "chrX":
        nx = comp.count("X")
        hs = [draw() for _ in range(max(nx, 2))]
        H = np.stack(hs[:max(nx, 2)]).astype(np.int8)
        if nx == 1:
            H[1, ~par] = -1                                           # one X: the second row is the Y's PAR copy, present in PAR only
            if "Y" not in comp:
                H[1, :] = -1                                          # 45,X: one copy everywhere, PAR included
        if nx == 2 and "Y" in comp:                                   # 47,XXY: three copies in the PAR
            H = np.vstack([H, np.where(par, draw(), -1).astype(np.int8)])
        return H
    return np.stack([draw(), draw()]).astype(np.int8)


def y_haps(n, draw, comp, par):
    """The Y: one homologue outside the PAR in a man (none in a woman); inside the PAR as many copies as the person has
    (X's plus the Y: two in XX and XY, three in XXY, one in 45,X), drawn fresh - the real export carries PAR probes on the Y
    with diploid calls in both sexes, and no step reads them as the Y."""
    n_par = comp.count("X") + comp.count("Y")
    rows = max(1, n_par)
    H = np.full((rows, n), -1, dtype=np.int8)
    if "Y" in comp:
        H[0] = draw()
    for r in range(n_par):
        H[r, par] = draw()[par]
    return H


def pick(H, rng):
    """One homologue per probe from a parent's valid homologues (Mendelian transmission; no linkage)."""
    n = H.shape[1]
    valid = H >= 0
    out = np.full(n, -1, dtype=np.int8)
    k = valid.sum(axis=0)
    r = rng.random(n)
    for h in range(H.shape[0]):
        # the h-th valid homologue is chosen with probability 1/k
        rank = np.cumsum(valid, axis=0)[h] - 1
        choose = valid[h] & (r < (rank + 1) / np.maximum(k, 1)) & (out < 0)
        out[choose] = H[h, choose]
    return out


def child_haps(pr, chrom, kid, dad, mom, rng, par, hidden_dad=None):
    """A child's constitutional homologues from the parents' (hidden_dad: the true father when the pedigree names another)."""
    n = pr["pos"].size
    D = (hidden_dad or dad).haps[chrom]
    M = mom.haps[chrom]
    if chrom == "chrY":
        draw = lambda: (rng.random(n) < pr["p"]).astype(np.int8)  # noqa: E731
        H = y_haps(n, draw, kid.complement, par)
        if kid.sex == "M":
            H[0, ~par] = D[0, ~par]                                       # the Y itself from the father; the PAR copies drawn fresh (no transmission modelled there)
        return H
    if chrom == "chrX":
        hm = pick(M, rng)
        if kid.sex == "M":
            h2 = np.full(n, -1, dtype=np.int8)
            h2[par] = pick(D, rng)[par]                                # the Y's PAR copy
            return np.stack([hm, h2])
        hd = pick(D, rng)
        return np.stack([hm, hd])
    return np.stack([pick(D, rng), pick(M, rng)])


# ---------------------------------------------------------------------------------------------------------------- events

def plant(person, pr, chrom, start_bp, end_bp, f, kind, label, cls, rng, hap=0):
    """A planted change on one homologue over [start_bp, end_bp] (bp, inclusive; None = chromosome end) in a cell fraction f."""
    pos = pr["pos"]
    a = int(np.searchsorted(pos, start_bp, side="left"))
    b = int(np.searchsorted(pos, end_bp if end_bp else pos[-1] + 1, side="right"))
    if b - a < 10:
        return None
    H = person.haps[chrom]
    present = (H[hap, a:b] >= 0).astype(np.int16)                      # the homologue exists at the probe (not in a man's X outside the PAR, row 1)
    hk = np.where(present > 0, H[hap, a:b], 0).astype(np.int16)
    if kind == "gain":
        dB, dC = hk, present
    elif kind == "loss":
        dB, dC = -hk, -present
    elif kind == "loh":                                                 # copy-neutral: the other homologue replaced by this one
        A = np.where(H[:, a:b] >= 0, H[:, a:b], 0).sum(axis=0)
        dB, dC = 2 * hk - A, np.zeros(b - a, dtype=np.int16)
    elif kind == "roh":                                                 # constitutional homozygosity: the homologues made identical
        other = 1 - hap
        H[other, a:b] = np.where(H[hap, a:b] >= 0, H[hap, a:b], H[other, a:b])
        dB, dC = None, 0
    else:
        raise ValueError(kind)
    ev = dict(sample=person.name, chrom=chrom, a=a, b=b, start=int(pos[a]), end=int(pos[b - 1]), f=float(f), dB=dB, dC=dC, type=kind, label=label, cls=cls, n_probes=int(b - a))
    person.events.append(ev)
    return ev


def iso_upd(person, chrom, hap):
    """A whole-chromosome isodisomy: both homologues the same one."""
    H = person.haps[chrom]
    H[1 - hap] = H[hap]


# ------------------------------------------------------------------------------------------------------------------ signal

def signal(person, pr, chrom, rng, batch_offset, probe_offset, gc_z, gc1mb_z, noise, cf):
    """GT (int8: 0/1/2, -1), LRR, BAF (unclipped; NaN where missing) of one person on one chromosome."""
    n = pr["pos"].size
    H = person.haps[chrom]
    valid = H >= 0
    A = np.where(valid, H, 0).sum(axis=0).astype(np.float64)
    C = valid.sum(axis=0).astype(np.float64)
    meanB, meanC = A.copy(), C.copy()
    for ev in person.events:
        if ev["chrom"] != chrom or ev["dB"] is None:
            continue
        a, b, f = ev["a"], ev["b"], min(ev["f"] * cf, 1.0)
        meanB[a:b] = (1 - f) * A[a:b] + f * (A[a:b] + ev["dB"])
        meanC[a:b] = (1 - f) * C[a:b] + f * (C[a:b] + ev["dC"])
    if chrom == "chrY":
        lrr_mean = np.where(pr["par"], level(meanC), level_y(meanC, rng))
    else:
        lrr_mean = level(meanC)
    baf_true = np.where(meanC > 0, meanB / np.maximum(meanC, 1e-9), np.nan)
    het = (baf_true > 0.02) & (baf_true < 0.98)
    sd = np.where(het, person.baf_sd * noise, 0.015 * noise)
    baf = baf_true + rng.normal(0, 1, n) * sd
    lrr = lrr_mean + person.shift + probe_offset + batch_offset + person.wave * gc_z + person.wave2 * gc1mb_z + rng.normal(0, person.lrr_sd * noise, n)
    # the planted probe failures that act on the signal
    for i, why in pr["fail"].items():
        if why == "lrr_median":
            lrr[i] += 0.5 if (i % 2) else -0.5
        elif why == "baf_shift" and het[i]:
            baf[i] += 0.12
    # genotype calls as GenCall would make them, with errors and no-calls
    gt = np.where(np.isnan(baf), -1, np.where(baf < 0.15, 0, np.where(baf > 0.85, 2, 1))).astype(np.int8)
    err = rng.random(n) < 0.002
    gt[err] = ((gt[err] + rng.integers(1, 3, err.sum())) % 3).astype(np.int8)
    nocall = rng.random(n) < person.nocall
    gt[nocall] = -1
    for i, why in pr["fail"].items():
        if why == "call_rate" and rng.random() < 0.2:
            gt[i] = -1
    io = pr["kind"] == 1
    gt[io] = -1
    baf[io] = np.nan
    if chrom == "chrY":
        absent = (meanC < 0.05) & ~pr["par"]
        gt[absent] = -1
        baf[absent] = np.nan
    miss = rng.random(n) < 0.001
    lrr[miss] = np.nan
    baf[miss] = np.nan
    return gt, lrr.astype(np.float32), baf.astype(np.float32)


# ------------------------------------------------------------------------------------------------------------------ cohort

def build(out, n_samples=60, seed=1, genome="chm13", contigs=DEFAULT_CONTIGS, full=False, density_kb=12.0, trios=4, duos=2, events=True, families=True,
          qc_failures=True, noise=1.0, cell_fraction=1.0, extra_events=None, lrr_sd=None, log_fn=log):
    """extra_events: a list of dicts (chrom, start, end (None: the chromosome end), f, kind (gain/loss/loh/roh), label, cls) each planted
    in its own singleton sample (the calibration study); lrr_sd: a fixed per-sample LRR noise instead of the 0.10-0.19 draw."""
    t0 = time.time()
    rng = np.random.default_rng(seed)
    G = load_genome(genome)
    chroms = G.chroms if full else [c for c in contigs.split(",") if c in G.length]
    if n_samples < MIN_SAMPLES:
        raise SystemExit("the mock needs at least %d samples for its planted roles (--samples)" % MIN_SAMPLES)
    os.makedirs(out, exist_ok=True)
    names = ["S%03d" % (i + 1) for i in range(n_samples)]
    # people and roles
    people = []
    fam_rows, families_truth = [], []
    i = 0
    if not families:
        trios = duos = 0
    for t in range(trios):
        dad, mom, kid = Person(names[i], "M", "dad"), Person(names[i + 1], "F", "mom"), Person(names[i + 2], rng.choice(["M", "F"]), "kid")
        kid.dad, kid.mom = dad, mom
        people += [dad, mom, kid]
        families_truth.append(dict(kid=kid.name, dad=dad.name, mom=mom.name, kind="trio", note=""))
        i += 3
    for d in range(duos):
        par = Person(names[i], rng.choice(["M", "F"]), "duo_parent")
        kid = Person(names[i + 1], rng.choice(["M", "F"]), "duo_kid")
        if par.sex == "M":
            kid.dad, kid.mom = par, Person("hidden_mom_%d" % d, "F")
        else:
            kid.mom, kid.dad = par, Person("hidden_dad_%d" % d, "M")
        people += [par, kid]
        families_truth.append(dict(kid=kid.name, dad=par.name if par.sex == "M" else "", mom=par.name if par.sex == "F" else "", kind="duo", note=""))
        i += 2
    singles = []
    while i < n_samples:
        p = Person(names[i], rng.choice(["M", "F"]))
        people.append(p)
        singles.append(p)
        i += 1
    by_name = {p.name: p for p in people}
    # per-sample noise and batch
    for k, p in enumerate(people):
        p.lrr_sd = float(rng.uniform(0.10, 0.19)) if lrr_sd is None else float(lrr_sd)
        p.baf_sd = float(rng.uniform(0.03, 0.05))
        p.wave = float(rng.normal(0, 1) * 0.010)
        p.wave2 = float(rng.normal(0, 1) * 0.008)
        p.shift = float(rng.normal(0, 0.02))
        p.nocall = float(rng.uniform(0.002, 0.008))
        p.batch = "plate1" if k % 2 == 0 else "plate2"
    # the roles
    roles = {}
    pool = [p for p in singles]
    rng.shuffle(pool)

    def take(sex=None):
        for q in pool:
            if sex is None or q.sex == sex:
                pool.remove(q)
                return q
        raise SystemExit("not enough singleton samples for the planted roles: raise --samples")

    truth_events, sex_truth, swap, nonpat, dup, noisy = [], [], None, None, None, []
    if events:
        roles["trisomy21"] = take(); roles["mosaic_gain_12"] = take(); roles["mosaic_loss_7q"] = take(); roles["mosaic_loh_11p"] = take()
        roles["del_1p36"] = take(); roles["roh_15mb"] = take()
        roles["xxy"] = take("M"); roles["x0"] = take("F"); roles["loy"] = take("M"); roles["mosaic_x_loss"] = take("F")
        if trios >= 1:
            roles["upd15_iso"] = people[2]                      # the first trio's child
        roles["xxy"].complement = "XXY"
        roles["x0"].complement = "X"
    if families and trios >= 2:
        woman = take("F")
        swap = (people[3].name, woman.name)                     # the second trio's father and an unrelated woman exchange columns
        roles["swap_dad"], roles["swap_woman"] = people[3], woman
    if families and trios >= 3:
        nonpat = dict(kid=people[8].name, claimed_dad=people[6].name)
        roles["nonpat_kid"] = people[8]
    if qc_failures:
        a, b = take(), take()
        b.copy_of = a                                            # the same person twice: the same haplotypes, the same sex, noise drawn again
        b.sex, b.complement = a.sex, a.complement
        a.role, b.role = "duplicate", "duplicate"
        dup = (a.name, b.name)
        n1, n2 = take(), take()
        n1.lrr_sd, n2.lrr_sd = 0.36, 0.45
        noisy = [dict(sample=n1.name, reason="lrr_sd_fail", value=0.36), dict(sample=n2.name, reason="lrr_sd_fail", value=0.45)]
        w1, w2 = take(), take()
        w1.wave, w1.wave2 = 0.06, 0.09
        w2.wave, w2.wave2 = -0.05, -0.11
        noisy += [dict(sample=w1.name, reason="wave", value=0.09), dict(sample=w2.name, reason="wave", value=-0.11)]
        b1 = take()
        b1.baf_sd = 0.09
        noisy.append(dict(sample=b1.name, reason="baf_sd", value=0.09))
        c1 = take()
        c1.nocall = 0.05
        noisy.append(dict(sample=c1.name, reason="call_rate", value=0.95))
        for q, r in ((n1, "lrr_sd_fail"), (n2, "lrr_sd_fail"), (w1, "wave"), (w2, "wave"), (b1, "baf_sd"), (c1, "call_rate")):
            q.role = (q.role + "," if q.role else "") + "qc_" + r
    for r, p in roles.items():
        p.role = (p.role + "," if p.role else "") + r
    # probes
    probes = {c: make_probes(G, c, density_kb, rng, qc_failures) for c in chroms}
    # haplotypes: founders first (including hidden parents), then children by generation
    hidden = {}
    for p in people:
        for q in (p.dad, p.mom):
            if q is not None and q.name not in by_name:
                hidden[q.name] = q
    hidden_dad = None
    if nonpat:
        hidden_dad = Person("hidden_true_dad", "M")
        hidden["hidden_true_dad"] = hidden_dad
    for c in chroms:
        pr = probes[c]
        par = pr["par"]
        for p in list(hidden.values()) + [q for q in people if q.dad is None and q.copy_of is None]:
            p.haps[c] = founder_haps(pr, c, p, rng, par)
        for p in people:
            if p.dad is not None:
                hd = hidden_dad if (nonpat and p.name == nonpat["kid"]) else None
                p.haps[c] = child_haps(pr, c, p, p.dad, p.mom, rng, par, hidden_dad=hd)
        for p in people:
            if p.copy_of is not None:
                p.haps[c] = p.copy_of.haps[c].copy()
    # planted events
    if events:
        cf = cell_fraction
        def ev(role, chrom, start, end, f, kind, label, cls, hap=0):
            if role not in roles or chrom not in probes:
                return
            e = plant(roles[role], probes[chrom], chrom, start, end, f, kind, label, cls, rng, hap)
            if e:
                truth_events.append(e)
        ev("trisomy21", "chr21", 1, None, 1.0, "gain", "trisomy 21, constitutional", "whole")
        ev("mosaic_gain_12", "chr12", 1, None, 0.30, "gain", "mosaic gain of chromosome 12 in 30% of cells", "whole")
        ev("mosaic_loss_7q", "chr7", G.p_end["chr7"] + 1, None, 0.25, "loss", "mosaic loss of 7q in 25% of cells", "q")
        ev("mosaic_loh_11p", "chr11", 1, G.p_end["chr11"], 0.40, "loh", "mosaic copy-neutral LOH of 11p in 40% of cells", "p")
        ev("del_1p36", "chr1", 2_000_000, 12_000_000, 1.0, "loss", "constitutional 10-Mb deletion of 1p36", "stretch")
        ev("roh_15mb", "chr1", 150_000_000, 165_000_000, 1.0, "roh", "a 15-Mb run of homozygosity on 1q", "stretch")
        if "upd15_iso" in roles and "chr15" in probes:
            p = roles["upd15_iso"]
            iso_upd(p, "chr15", 1)                                          # both copies the maternal homologue
            pr = probes["chr15"]
            truth_events.append(dict(sample=p.name, chrom="chr15", a=0, b=pr["pos"].size, start=int(pr["pos"][0]), end=int(pr["pos"][-1]), f=1.0, dB=None, dC=0, type="roh",
                                     label="whole-chromosome isodisomy of 15 (maternal) in a trio child", cls="whole", n_probes=int(pr["pos"].size)))
        if "chrY" in probes and "loy" in roles:
            p = roles["loy"]
            pr = probes["chrY"]
            e = dict(sample=p.name, chrom="chrY", a=0, b=pr["pos"].size, start=int(pr["pos"][0]), end=int(pr["pos"][-1]), f=0.5, dB=np.zeros(pr["pos"].size, dtype=np.int16),
                     dC=-(p.haps["chrY"][0] >= 0).astype(np.int16),
                     type="loss", label="mosaic loss of Y in 50% of cells", cls="whole", n_probes=int(pr["pos"].size))
            p.events.append(e)
            sex_truth.append(dict(sample=p.name, chrom="chrY", type="loss", f=0.5, label=e["label"]))
            if "chrX" in probes:                                             # the Y's PAR copy goes with it
                plant(p, probes["chrX"], "chrX", 1, None, 0.5, "loss", "the Y's PAR copy lost with the Y", "par", rng, hap=1)
                p.events[-1]["par_only"] = True
        if "chrX" in probes and "mosaic_x_loss" in roles:
            e = plant(roles["mosaic_x_loss"], probes["chrX"], "chrX", 1, None, 0.30, "loss", "mosaic loss of an X in 30% of cells of a woman", "whole", rng, hap=1)
            truth_events.remove(e) if e in truth_events else None
            sex_truth.append(dict(sample=roles["mosaic_x_loss"].name, chrom="chrX", type="loss", f=0.30, label=e["label"]))
        if "xxy" in roles:
            sex_truth.append(dict(sample=roles["xxy"].name, chrom="chrX", type="gain", f=1.0, label="47,XXY"))
        if "x0" in roles:
            sex_truth.append(dict(sample=roles["x0"].name, chrom="chrX", type="loss", f=1.0, label="45,X"))
    for k, x in enumerate(extra_events or []):
        if x["chrom"] not in probes:
            continue
        q = take()
        q.role = (q.role + "," if q.role else "") + "extra_%d" % k
        e = plant(q, probes[x["chrom"]], x["chrom"], x["start"], x.get("end"), x["f"], x["kind"], x.get("label", "planted"), x.get("cls", "stretch"), rng)
        if e:
            truth_events.append(e)
    # the export columns: the person's own name, except the swapped pair
    column_of = {p.name: p.name for p in people}
    if swap:
        column_of[swap[0]], column_of[swap[1]] = swap[1], swap[0]
    order = sorted(people, key=lambda p: names.index(column_of[p.name]))     # export order by column name
    # per-batch and per-probe offsets, then the VCF
    batches = sorted({p.batch for p in people})
    vcf_path = os.path.join(out, "mock.vcf.gz")
    h = header(G, order, column_of)
    fh = pysam.BGZFile(vcf_path, "wb")
    fh.write(h.encode())
    n_records = {}
    for c in chroms:                                                            # the wave is a genome-wide function of GC: standardise over all probes
        probes[c]["gc1mb"] = window_mean(probes[c]["pos"], probes[c]["gc"].astype(float), 500_000)
    all_gc = np.concatenate([probes[c]["gc"].astype(float) for c in chroms])
    all_gc1 = np.concatenate([probes[c]["gc1mb"] for c in chroms])
    gc_mean, gc_sd, gc1_mean, gc1_sd = all_gc.mean(), all_gc.std() + 1e-9, all_gc1.mean(), all_gc1.std() + 1e-9
    for c in chroms:
        pr = probes[c]
        n = pr["pos"].size
        probe_offset = rng.normal(0, 0.08, n)
        batch_offset = {b: rng.normal(0, 0.04, n) for b in batches}
        gc = pr["gc"].astype(float)
        gc_z = (gc - gc_mean) / gc_sd
        gc1_z = (pr["gc1mb"] - gc1_mean) / gc1_sd
        cols = []
        for p in order:
            gt, lrr, baf = signal(p, pr, c, rng, batch_offset[p.batch], probe_offset, gc_z, gc1_z, noise, cell_fraction)
            cols.append((gt, lrr, baf))
        if qc_failures:
            for i, why in pr["fail"].items():
                if why == "hwe":
                    for gt, lrr, baf in cols:
                        if gt[i] >= 0:
                            gt[i] = 1
                            baf[i] = 0.5 + rng.normal(0, 0.04)
        write_records(fh, pr, cols)
        n_records[c] = int(n)
        log_fn("%s: %d probes" % (c, n), "mock")
    fh.close()
    pysam.tabix_index(vcf_path, preset="vcf", force=True)
    bcf_path = os.path.join(out, "mock.bcf")
    bcftools = shutil.which("bcftools")
    if bcftools:
        subprocess.run([bcftools, "view", "-Ob", "-o", bcf_path, vcf_path], check=True)
        subprocess.run([bcftools, "index", "-c", "-f", bcf_path], check=True)
    else:
        log_fn("bcftools not on PATH: mock.bcf not written (mock.vcf.gz is indexed)", "mock")
    # the sample sheet (what the lab believes, by column name) and the pedigree
    sheet = []
    for p in people:
        col = column_of[p.name]
        believed = by_name[col]                                                  # the person the lab believes sits in that column
        sheet.append(dict(sample=col, reported_sex=believed.sex, batch=p.batch, family=family_of(believed, families_truth), note=believed.role))
    sheet.sort(key=lambda r: names.index(r["sample"]))
    write_tsv(os.path.join(out, "sample_sheet.tsv"), ["sample", "reported_sex", "batch", "family", "note"], sheet)
    ped = []
    for f in families_truth:
        kid, dad, mom = by_name[f["kid"]], by_name.get(f["dad"]), by_name.get(f["mom"])
        ped.append(dict(kid=f["kid"], dad=f["dad"], mom=f["mom"], kid_sex=kid.sex, dad_sex=dad.sex if dad else "", mom_sex=mom.sex if mom else ""))
    write_pedigree(os.path.join(out, "pedigree.tsv"), ped)
    # the truth, by export column: what the data under each name carries
    samples = []
    for p in people:
        col = column_of[p.name]
        comp = p.complement
        samples.append(dict(sample=col, person=p.name, sex=p.sex, complement=comp, batch=p.batch, role=p.role, lrr_sd=p.lrr_sd, baf_sd=p.baf_sd, wave=p.wave, wave2=p.wave2,
                            nocall=p.nocall, reported_sex=by_name[col].sex, family=family_of(p, families_truth), duplicate_of=column_of.get(p.copy_of.name) if p.copy_of else ""))
    samples.sort(key=lambda r: names.index(r["sample"]))
    ev_out = []
    for e in truth_events:
        if e.get("par_only"):
            continue
        ev_out.append(dict(sample=column_of[e["sample"]], chrom=e["chrom"], start=e["start"], end=e["end"], n_probes=e["n_probes"], type=e["type"], f=e["f"] * (cell_fraction if e["f"] < 1 else 1),
                           **{"class": e["cls"]}, label=e["label"], bands=G.band_range(e["chrom"], e["start"], e["end"])))
    for e in sex_truth:
        e["sample"] = column_of[e["sample"]]
    affected = {e["sample"] for e in ev_out} | {e["sample"] for e in sex_truth} | {column_of[x["sample"]] for x in noisy}
    if swap:
        affected |= set(swap)
    if dup:
        affected |= {column_of[dup[0]], column_of[dup[1]]}
    nulls = [r["sample"] for r in samples if r["sample"] not in affected and r["complement"] in ("XX", "XY")]
    probe_fail = {}
    for c in chroms:
        for i, why in probes[c]["fail"].items():
            probe_fail[str(probes[c]["id"][i])] = why
    truth = dict(version=__version__, seed=seed, genome=G.name, contigs=chroms, density_kb=density_kb, n_samples=n_samples, samples=samples, families=families_truth,
                 events=ev_out, sex_events=sex_truth, swap=list(swap) if swap else [], non_paternity=nonpat or {}, duplicate=[column_of[dup[0]], column_of[dup[1]]] if dup else [],
                 qc_fail=[dict(sample=column_of[x["sample"]], reason=x["reason"], value=x["value"]) for x in noisy], probe_fail=probe_fail, null_samples=nulls,
                 levels=dict(one_copy=float(LEVELS[1]), three_copies=float(LEVELS[3]), y_absent=Y_ABSENT), noise=noise, cell_fraction=cell_fraction,
                 records=n_records, export=dict(vcf="mock.vcf.gz", bcf="mock.bcf" if bcftools else ""), seconds=time.time() - t0)
    write_json(truth, os.path.join(out, "truth.json"))
    with open(os.path.join(out, "MOCK_DATA.txt"), "w") as fhm:
        fhm.write("MOCK DATA: simulated by arraykaryo mock %s (seed %d); no real sample anywhere in this directory.\n" % (__version__, seed))
    log_fn("%d samples, %d probes on %d chromosomes, %d events, %d families -> %s (%.1f s)" % (n_samples, sum(n_records.values()), len(chroms), len(ev_out), len(families_truth), out, time.time() - t0), "mock")
    return truth


def family_of(p, fams):
    for k, f in enumerate(fams):
        if p.name in (f["kid"], f["dad"], f["mom"]):
            return "F%d" % (k + 1)
    return ""


def window_mean(pos, v, half):
    """The mean of v over the probes within +-half bp of each probe (a cumulative-sum window)."""
    cs = np.concatenate([[0.0], np.cumsum(v)])
    a = np.searchsorted(pos, pos - half, side="left")
    b = np.searchsorted(pos, pos + half, side="right")
    return (cs[b] - cs[a]) / np.maximum(b - a, 1)


def header(G, order, column_of):
    lines = ["##fileformat=VCFv4.2", "##source=arraykaryo mock %s (simulated data; no real sample)" % __version__]
    for c in G.chroms:
        lines.append("##contig=<ID=%s,length=%d>" % (c, G.length[c]))
    lines += ['##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description="A allele">', '##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description="B allele">',
              '##INFO=<ID=INTENSITY_ONLY,Number=0,Type=Flag,Description="Locus with intensity information only and no genotyping information">',
              '##INFO=<ID=GC,Number=1,Type=Float,Description="GC ratio content around the variant">',
              '##INFO=<ID=GenTrain_Score,Number=1,Type=Float,Description="The SNP cluster quality from the GenTrain clustering algorithm">',
              '##INFO=<ID=Cluster_Sep,Number=1,Type=Float,Description="The cluster separation measurement for the SNP that ranges between 0 and 1">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', '##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R Ratio">',
              '##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">']
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(column_of[p.name] for p in order))
    return "\n".join(lines) + "\n"


GT_STR = np.array(["0/0", "0/1", "1/1", "./."])


def fmt4(x):
    """Floats as '%.4f' strings, in bulk: integer arithmetic on the scaled values (np.char.mod is a Python loop per element)."""
    x = np.asarray(x, dtype=np.float64)
    nan = np.isnan(x)
    v = np.rint(np.where(nan, 0.0, x) * 10000.0).astype(np.int64)
    neg = v < 0
    v = np.abs(v)
    whole = (v // 10000).astype(str)
    frac = np.char.zfill((v % 10000).astype(str), 4)
    s = np.char.add(np.char.add(whole, "."), frac)
    s = np.where(neg, np.char.add("-", s), s)
    return np.where(nan, "nan", s)


def write_records(fh, pr, cols, block=1000):
    n = pr["pos"].size
    gt = np.stack([c[0] for c in cols], axis=1)
    lrr = np.stack([c[1] for c in cols], axis=1)
    baf = np.stack([c[2] for c in cols], axis=1)
    for i0 in range(0, n, block):
        i1 = min(i0 + block, n)
        gts = GT_STR[np.where(gt[i0:i1] < 0, 3, gt[i0:i1])]
        cells = np.char.add(np.char.add(np.char.add(np.char.add(gts, ":"), fmt4(lrr[i0:i1])), ":"), fmt4(baf[i0:i1]))
        rows = cells.tolist()
        out = []
        for k, i in enumerate(range(i0, i1)):
            info = "ALLELE_A=%d;ALLELE_B=%d;GC=%.4f;GenTrain_Score=%.4f;Cluster_Sep=%.4f" % (i % 2, 1 - i % 2, pr["gc"][i], pr["gentrain"][i], pr["csep"][i])
            if pr["kind"][i] == 1:
                info += ";INTENSITY_ONLY"
            out.append("%s\t%d\t%s\t%s\t%s\t.\t.\t%s\tGT:LRR:BAF\t%s\n" % (pr["chrom"], pr["pos"][i], pr["id"][i], pr["ref"][i], pr["alt"][i], info, "\t".join(rows[k])))
        fh.write("".join(out).encode())
