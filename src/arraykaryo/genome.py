"""Assemblies an array export may be aligned to, recognised by their contig lengths; the pseudoautosomal regions per assembly.

An export names its contigs in its header (##contig=<ID=chr1,length=...>). The lengths identify the assembly: a GRCh37 export
(hg19, b37 with or without the chr prefix), a GRCh38 one, or a T2T-CHM13v2.0 one. Nothing here is fetched: the tables are
constants. VERIFY: the pseudoautosomal coordinates below are the published ones to the author's knowledge and must be
checked against the assembly's own annotation before they gate anything (the GRCh38 and GRCh37 ones are the UCSC/Ensembl
PAR definitions; the CHM13 ones are the T2T consortium's). A cheap empirical check: a male's heterozygous calls on the X
are common inside a PAR and near absent beyond it."""

# chr1, chrX, chrY lengths: enough to tell the three apart
ASSEMBLIES = {
    "grch37": dict(chr1=249250621, chrX=155270560, chrY=59373566),
    "grch38": dict(chr1=248956422, chrX=156040895, chrY=57227415),
    "chm13": dict(chr1=248387328, chrX=154259566, chrY=62460029),
}

# pseudoautosomal regions on the X, 1-based inclusive (VERIFY, see above)
PAR = {
    "grch37": (("chrX", 60001, 2699520), ("chrX", 154931044, 155260560)),
    "grch38": (("chrX", 10001, 2781479), ("chrX", 155701383, 156030895)),
    "chm13": (("chrX", 1, 2394410), ("chrX", 153925835, 154259566)),
}

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]


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
    """Whether a position lies in a pseudoautosomal region of the assembly (False for an unknown assembly)."""
    for c, a, b in PAR.get(assembly, ()):
        if normalize_contig(contig) == c and a <= pos <= b:
            return True
    return False
