"""Assemblies an array export may be aligned to, recognised by their contig lengths, and what each carries: chromosome
lengths, the end of every short arm (the centromere, the p/q boundary of the acen bands), the pseudoautosomal regions on
the X and the Y, the acrocentric short arms (no probes: an arm with no probes is not a loss), the cytogenetic bands, and
the regions where array signal is known to mislead (the MHC, the 8p23.1 and 17q21.31 inversions, the pericentromeric and
subtelomeric flanks), resolved per assembly from the band table.

GRCh37 (hg19, b37, with or without the chr prefix), GRCh38 (hg38) and T2T-CHM13v2.0 are vendored: data/cytoband.<asm>.tsv
holds each assembly's bands (UCSC hg19/hg38 cytoBand; the T2T consortium's chm13v2.0 cytobands - scripts/vendor_genome_tables.py
records the sources and checksums), from which the lengths, the centromeres and the acrocentric arms follow. The
pseudoautosomal regions are typed below from UCSC's hg19/hg38 `par` tables and the T2T consortium's chm13v2.0_PAR.bed,
1-based inclusive, and were checked empirically on the public example: a male's heterozygous X calls are frequent inside
PAR1 and near absent beyond it (docs/METHODS.md). Any other assembly is given as a table in the same format by
--genome-table (optional `##par` lines for its pseudoautosomal regions)."""
import os

import numpy as np

# chr1, chrX, chrY lengths: enough to tell the three apart
ASSEMBLIES = {
    "grch37": dict(chr1=249250621, chrX=155270560, chrY=59373566),
    "grch38": dict(chr1=248956422, chrX=156040895, chrY=57227415),
    "chm13": dict(chr1=248387328, chrX=154259566, chrY=62460029),
}
ALIASES = {"grch37": "grch37", "hg19": "grch37", "b37": "grch37", "grch38": "grch38", "hg38": "grch38", "chm13": "chm13", "t2t": "chm13",
           "chm13v2": "chm13", "chm13v2.0": "chm13", "t2t-chm13v2.0": "chm13", "t2t-chm13": "chm13"}
LABEL = {"grch37": "GRCh37 (hg19)", "grch38": "GRCh38 (hg38)", "chm13": "T2T-CHM13v2.0"}

# pseudoautosomal regions, 1-based inclusive, on the X and on the Y (the Y copy of each is written at the Y's own coordinates)
PAR = {
    "grch37": {"chrX": ((60001, 2699520), (154931044, 155260560)), "chrY": ((10001, 2649520), (59034050, 59363566))},
    "grch38": {"chrX": ((10001, 2781479), (155701383, 156030895)), "chrY": ((10001, 2781479), (56887903, 57217415))},
    "chm13": {"chrX": ((1, 2394410), (153925835, 154259566)), "chrY": ((1, 2458320), (62122810, 62460029))},
}

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]
CHROMS = AUTOSOMES + ["chrX", "chrY"]
ACROCENTRIC = ("chr13", "chr14", "chr15", "chr21", "chr22")

# regions where array signal is known to mislead, by band so that they resolve on any assembly with a band table
ARTEFACT_BANDS = (
    ("MHC", "chr6", "p22.1", "p21.32", "the major histocompatibility complex: dense structural variation; probe clusters shift with the haplotype"),
    ("8p23.1 inversion", "chr8", "p23.1", "p23.1", "a common inversion polymorphism flanked by segmental duplications; copy-number-variable flanks"),
    ("17q21.31 inversion", "chr17", "q21.31", "q21.31", "the H1/H2 inversion polymorphism (MAPT) with copy-number-variable flanks"),
    ("X-transposed region", "chrX", "q21.31", "q21.32", "99% identical to a stretch of Yp: in a man its probes read the X and the Y copies (two copies, heterozygous-looking calls)"),
)
PERICENTROMERIC_BP = 2_000_000     # each side of the centromere
SUBTELOMERIC_BP = 1_000_000        # at each chromosome end


def normalize_contig(name):
    """'1', 'chr1', 'X', 'chrx', 'MT' -> 'chr1', 'chrX', 'chrM'; a contig that is none of chr1-22, X, Y, M is returned as given."""
    c = str(name)
    if not c.startswith("chr"):
        c = "chr" + c
    tail = c[3:]
    if tail.upper() in ("X", "Y", "M", "MT"):
        tail = "M" if tail.upper() in ("M", "MT") else tail.upper()
        c = "chr" + tail
    return c


def guess_assembly(contig_lengths):
    """The assembly whose chr1, chrX and chrY lengths all match the header's; '' when none does (an unknown assembly, or a header
    without lengths). contig_lengths: {name: length}, any naming."""
    lengths = {normalize_contig(k): v for k, v in contig_lengths.items() if v}
    for name, ref in ASSEMBLIES.items():
        if all(lengths.get(c) == L for c, L in ref.items()):
            return name
    return ""


def is_par(assembly, contig, pos):
    """Whether a position lies in a pseudoautosomal region of the assembly, on the X or the Y (False for an unknown assembly)."""
    c = normalize_contig(contig)
    for a, b in PAR.get(assembly, {}).get(c, ()):
        if a <= pos <= b:
            return True
    return False


def data_path(name):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", name)


def read_band_table(path):
    """A cytoband table (chrom, start 0-based, end exclusive, band, stain; '#' comments and a header allowed) -> ({chrom: [(start, end,
    band, stain)]}, {chrom: [(a, b)]} from any '##par<TAB>chrom<TAB>a<TAB>b' lines, 1-based inclusive)."""
    bands, par = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##par"):
                _, c, a, b = line.rstrip("\n").split("\t")[:4]
                par.setdefault(normalize_contig(c), []).append((int(a), int(b)))
                continue
            if line.startswith("#") or line.startswith("chrom\t") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            c, s, e, b = parts[:4]
            stain = parts[4] if len(parts) > 4 else ""
            bands.setdefault(normalize_contig(c), []).append((int(s), int(e), b, stain))
    for c in bands:
        bands[c].sort()
    return bands, par


class Genome:
    """One assembly: what every step asks of it. Positions 1-based inclusive at the interface; the band table is 0-based
    half-open inside, as the sources are."""

    def __init__(self, name, bands, par, label=""):
        self.name = name
        self.label = label or name
        self.bands = bands
        order = {c: i for i, c in enumerate(CHROMS)}
        self.chroms = sorted(bands, key=lambda c: (order.get(c, 99), c))
        self.length = {c: bands[c][-1][1] for c in self.chroms}
        self.p_end, self.centromere = {}, {}
        for c in self.chroms:
            acen = [r for r in bands[c] if r[3] == "acen"]
            if acen:
                p = [r for r in acen if r[2].startswith("p")]
                self.p_end[c] = p[-1][1] if p else acen[0][1]              # the last base of the short arm (1-based) = the p acen band's 0-based end
                self.centromere[c] = (acen[0][0] + 1, acen[-1][1])
            else:                                                            # no acen band (a custom table): the middle, so that every event is a stretch
                self.p_end[c] = 0
                self.centromere[c] = (0, 0)
        self.par = {normalize_contig(c): tuple(v) for c, v in par.items()}

    @property
    def autosomes(self):
        return [c for c in self.chroms if c not in ("chrX", "chrY", "chrM")]

    def normalize(self, contig):
        """A contig name as the genome names it, or None for one it does not carry (alts, decoys, chrM)."""
        c = normalize_contig(contig)
        return c if c in self.length else None

    def is_par(self, chrom, pos):
        return any(a <= pos <= b for a, b in self.par.get(chrom, ()))

    def par_mask(self, chrom, positions):
        """A boolean array over 1-based positions: inside a pseudoautosomal region of this chromosome."""
        pos = np.asarray(positions)
        m = np.zeros(pos.shape, dtype=bool)
        for a, b in self.par.get(chrom, ()):
            m |= (pos >= a) & (pos <= b)
        return m

    def arm(self, chrom, pos):
        """'p' or 'q' for a 1-based position ('' when the chromosome has no centromere)."""
        if not self.p_end.get(chrom):
            return ""
        return "p" if pos <= self.p_end[chrom] else "q"

    def arm_bounds(self, chrom, arm):
        """The 1-based inclusive span of an arm."""
        if arm == "p":
            return 1, self.p_end[chrom]
        return self.p_end[chrom] + 1, self.length[chrom]

    def is_acrocentric_p(self, chrom, arm):
        return arm == "p" and chrom in ACROCENTRIC

    def band(self, chrom, pos):
        """The cytogenetic band holding a 1-based position ('q22.1'); '' outside the table."""
        rows = self.bands.get(chrom)
        if not rows:
            return ""
        starts = np.fromiter((r[0] for r in rows), dtype=np.int64, count=len(rows))
        i = int(np.searchsorted(starts, pos - 1, side="right")) - 1
        if i < 0 or not (rows[i][0] < pos <= rows[i][1]):
            return ""
        return rows[i][2]

    def band_range(self, chrom, start, end):
        """The bands spanned by [start, end] (1-based inclusive) in ISCN form: 'q22.1q31.1', or one band name; '' without a table."""
        a, b = self.band(chrom, start), self.band(chrom, max(end, start))
        if not a or not b:
            return ""
        return a if a == b else a + b

    def band_span(self, chrom, band_from, band_to):
        """The 1-based inclusive span from the start of band_from to the end of band_to; None when a band is missing."""
        rows = self.bands.get(chrom, [])
        s = [r for r in rows if r[2] == band_from]
        e = [r for r in rows if r[2] == band_to]
        if not s or not e:
            return None
        return s[0][0] + 1, e[-1][1]

    def artefacts(self):
        """The regions where an event is tagged (never removed): [(name, chrom, start, end, why)], 1-based inclusive."""
        out = []
        for name, c, b0, b1, why in ARTEFACT_BANDS:
            if c in self.bands:
                sp = self.band_span(c, b0, b1)
                if sp:
                    out.append((name, c, sp[0], sp[1], why))
        for c in self.chroms:
            if c == "chrM":
                continue
            a, b = self.centromere[c]
            if b:
                out.append(("pericentromeric", c, max(1, a - PERICENTROMERIC_BP), min(self.length[c], b + PERICENTROMERIC_BP),
                            "the centromere and %d kb each side: satellite repeats, few probes, mis-mapped signal" % (PERICENTROMERIC_BP // 1000)))
            out.append(("subtelomeric", c, 1, min(self.length[c], SUBTELOMERIC_BP), "the first %d kb: probe-sparse, GC-rich" % (SUBTELOMERIC_BP // 1000)))
            out.append(("subtelomeric", c, max(1, self.length[c] - SUBTELOMERIC_BP + 1), self.length[c], "the last %d kb: probe-sparse, GC-rich" % (SUBTELOMERIC_BP // 1000)))
            if c in ACROCENTRIC and self.p_end.get(c):
                out.append(("acrocentric p arm", c, 1, self.p_end[c], "the short arm of an acrocentric chromosome (ribosomal DNA, no probes): an arm with no probes is not a loss"))
        return out

    def chrom_order(self):
        return {c: i for i, c in enumerate(self.chroms)}


def load_genome(name_or_path):
    """A Genome by name (grch37/hg19/b37, grch38/hg38, chm13/t2t) from the vendored tables, or from a custom band table path."""
    key = ALIASES.get(str(name_or_path).lower())
    if key:
        bands, _ = read_band_table(data_path("cytoband.%s.tsv" % key))
        return Genome(key, bands, PAR[key], LABEL[key])
    if os.path.exists(str(name_or_path)):
        bands, par = read_band_table(name_or_path)
        return Genome(os.path.splitext(os.path.basename(name_or_path))[0], bands, par, "custom table %s" % os.path.basename(name_or_path))
    raise ValueError("unknown genome %r: grch37 (hg19), grch38 (hg38) and chm13 are built in; pass a band table path for another" % (name_or_path,))


def genome_for(contig_lengths, override=""):
    """The genome an export is on: `override` (a name or a table path) when given, else the assembly recognised from the header's
    contig lengths. Raises with a clear message when neither identifies one."""
    if override:
        return load_genome(override)
    asm = guess_assembly(contig_lengths)
    if not asm:
        raise ValueError("the export's contig lengths match none of GRCh37, GRCh38 and T2T-CHM13v2.0 (chr1 %s, chrX %s, chrY %s); name the assembly "
                         "or a band table with --genome" % tuple(contig_lengths.get(c, contig_lengths.get(c[3:], "absent")) for c in ("chr1", "chrX", "chrY")))
    return load_genome(asm)
