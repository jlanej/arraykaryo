"""Vendor the per-assembly genome tables arraykaryo ships (src/arraykaryo/data/): the cytogenetic bands (from which the
chromosome lengths, the centromeres and the acrocentric short arms follow) for GRCh37, GRCh38 and T2T-CHM13v2.0.

Sources (public; the files are read as downloaded, no edit):
  GRCh37  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz   (UCSC hg19 cytoBand)
  GRCh38  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz   (UCSC hg38 cytoBand)
  CHM13   https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed
          (T2T consortium, CHM13v2.0 cytobands)
The pseudoautosomal regions are typed into genome.py from UCSC's hg19/hg38 `par` tables and the T2T consortium's
chm13v2.0_PAR.bed (analysis_set), and verified empirically on the public example (docs/METHODS.md).

  python scripts/vendor_genome_tables.py --grch37 cytoBand.hg19.txt.gz --grch38 cytoBand.hg38.txt.gz --chm13 chm13v2.0_cytobands_allchrs.bed[.gz]
"""
import argparse
import gzip
import hashlib
import os

CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
SOURCE = {
    "grch37": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
    "grch38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
    "chm13": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed",
}


def read(path):
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            c, s, e, b, stain = line.rstrip("\n").split("\t")[:5]
            if c in CHROMS:
                rows.append((c, int(s), int(e), b, stain))
    rows.sort(key=lambda r: (CHROMS.index(r[0]), r[1]))
    return rows


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grch37"), ap.add_argument("--grch38"), ap.add_argument("--chm13")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "arraykaryo", "data"))
    a = ap.parse_args()
    for asm in ("grch37", "grch38", "chm13"):
        path = getattr(a, asm)
        if not path:
            continue
        rows = read(path)
        out = os.path.join(a.out, "cytoband.%s.tsv" % asm)
        with open(out, "w") as fh:
            fh.write("# cytogenetic bands of %s: %s (sha256 of the file read: %s); the 24 chromosomes; start 0-based, end exclusive\n" % (asm, SOURCE[asm], sha256(path)))
            fh.write("chrom\tstart\tend\tband\tstain\n")
            for r in rows:
                fh.write("%s\t%d\t%d\t%s\t%s\n" % r)
        # a check: every chromosome's bands tile it from 0 without a gap, and the acen bands come in a p/q pair
        for c in CHROMS:
            b = [r for r in rows if r[0] == c]
            assert b[0][1] == 0 and all(b[i][2] == b[i + 1][1] for i in range(len(b) - 1)), c
            acen = [r for r in b if r[4] == "acen"]
            assert len(acen) == 2 and acen[0][3].startswith("p") and acen[1][3].startswith("q"), (c, acen)
        print("%s: %d bands -> %s (chr1 %d, chrX %d, chrY %d)" % (asm, len(rows), out, *[max(r[2] for r in rows if r[0] == c) for c in ("chr1", "chrX", "chrY")]))


if __name__ == "__main__":
    main()
