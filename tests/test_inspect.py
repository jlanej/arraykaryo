"""arraykaryo inspect on a tiny synthetic export written here (no real data): the header's samples, contigs and assembly, the
FORMAT gate, the conventions counted over the records (no-call, missing LRR/BAF, an intensity-only probe, BAF beyond 1)."""
import os
import shutil
import subprocess

import pysam
import pytest

from arraykaryo.cli import main
from arraykaryo.genome import guess_assembly, is_par, normalize_contig
from arraykaryo.inspect import inspect

CHM13 = dict(chr1=248387328, chr21=45090682, chrX=154259566, chrY=62460029)


def write_tiny(path, with_lrr=True, bcf=False):
    """Three samples, a dozen probes over chr1, chr21, chrX and chrY, in the gtc2vcf style: diploid GT everywhere, LRR and BAF per call,
    nan where missing, ./. for a no-call, one intensity-only probe without an ALT, a BAF slightly above 1."""
    h = pysam.VariantHeader()
    for c, L in CHM13.items():
        h.contigs.add(c, length=L)
    h.add_line('##INFO=<ID=INTENSITY_ONLY,Number=0,Type=Flag,Description="Locus with intensity information only and no genotyping information">')
    h.add_line('##INFO=<ID=GC,Number=1,Type=Float,Description="GC ratio content around the variant">')
    h.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if with_lrr:
        h.add_line('##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R Ratio">')
    h.add_line('##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">')
    for s in ("S1", "S2", "S3"):
        h.add_sample(s)
    out = pysam.VariantFile(path, "wb" if bcf else "wz", header=h)
    nan = float("nan")
    rows = [("chr1", 1000, "T", "C", [((0, 0), -0.1, 0.01), ((0, 1), 0.05, 0.49), ((1, 1), 0.2, 1.002)]),
            ("chr1", 2000, "A", "G", [((0, 1), 0.0, 0.52), ((None, None), nan, nan), ((0, 0), -0.3, 0.0)]),
            ("chr1", 3000, "G", None, [((None, None), 0.1, nan), ((None, None), 0.2, nan), ((None, None), -0.05, nan)]),   # intensity-only
            ("chr21", 5000, "C", "T", [((0, 0), 0.0, 0.0), ((0, 0), 0.1, 0.01), ((0, 1), 0.0, 0.5)]),
            ("chrX", 50_000_000, "A", "C", [((0, 0), -0.2, 1.0), ((0, 1), 0.3, 0.5), ((0, 0), 0.1, 0.99)]),
            ("chrY", 12_000_000, "T", "G", [((0, 0), 0.1, 0.0), ((None, None), -4.1, 0.6), ((None, None), -3.8, 0.3)])]
    for c, pos, ref, alt, vals in rows:
        rec = out.new_record(contig=c, start=pos - 1, stop=pos, alleles=(ref, alt if alt else "."), id="probe_%s_%d" % (c, pos))   # ALT "." as gtc2vcf writes an intensity-only probe
        rec.info["GC"] = 0.41
        if alt is None:
            rec.info["INTENSITY_ONLY"] = True
        for s, (gt, lrr, baf) in zip(("S1", "S2", "S3"), vals):
            rec.samples[s]["GT"] = gt
            if with_lrr:
                rec.samples[s]["LRR"] = lrr
            rec.samples[s]["BAF"] = baf
        out.write(rec)
    out.close()
    if bcf:
        subprocess.run([shutil.which("bcftools"), "index", "-c", "-f", path], check=True)
    else:
        pysam.tabix_index(path, preset="vcf", force=True)
    return path


def test_genome_tables():
    assert guess_assembly(CHM13) == "chm13" and guess_assembly(dict(chr1=248956422, chrX=156040895, chrY=57227415)) == "grch38"
    assert guess_assembly({"1": 249250621, "X": 155270560, "Y": 59373566}) == "grch37" and guess_assembly(dict(chr1=1)) == ""
    assert normalize_contig("1") == "chr1" and normalize_contig("x") == "chrX" and normalize_contig("MT") == "chrM" and normalize_contig("chrY") == "chrY"
    assert is_par("chm13", "chrX", 1_000_000) and not is_par("chm13", "chrX", 50_000_000) and is_par("grch38", "X", 155_800_000) and not is_par("", "chrX", 1)


@pytest.mark.parametrize("bcf", [False, True])
def test_inspect_reads_the_conventions(tmp_path, bcf):
    if bcf and not shutil.which("bcftools"):
        pytest.skip("bcftools not on PATH: the BCF is indexed with it")
    p = write_tiny(str(tmp_path / ("tiny.bcf" if bcf else "tiny.vcf.gz")), bcf=bcf)
    d = inspect(p)
    assert d["samples"] == 3 and d["assembly"] == "chm13" and d["chrom_prefixed"] is True and d["required_missing"] == []
    assert set(("GT", "LRR", "BAF")) <= set(d["formats"]) and "INTENSITY_ONLY" in d["info_of_interest"] and "GC" in d["info_of_interest"]
    s = d["sampled"]
    assert s["records"] == 6 and s["samples_read"] == 3 and s["by_contig"] == {"chr1": 3, "chr21": 1, "chrX": 1, "chrY": 1}
    assert s["gt_nocall"] == 6 and s["baf_missing"] == 4 and s["lrr_missing"] == 1 and s["intensity_only"] == 1 and s["alt_missing"] == 1
    assert s["baf_max"] > 1.0 and s["het"] == 4
    if d["records_per_contig"]:                                     # bcftools on PATH: the index agrees
        assert sum(d["records_per_contig"].values()) == 6
    assert main(["inspect", p, "--json", str(tmp_path / "i.json")]) == 0 and os.path.getsize(str(tmp_path / "i.json")) > 100


def test_inspect_refuses_an_export_without_lrr(tmp_path, capsys):
    p = write_tiny(str(tmp_path / "nolrr.vcf.gz"), with_lrr=False)
    assert main(["inspect", p]) == 2
    assert "REFUSED" in capsys.readouterr().out


def test_unimplemented_commands_say_where_they_are_specified(tmp_path):
    with pytest.raises(SystemExit) as e:
        main(["report", "--out", str(tmp_path)])
    assert "PROMPT.md" in str(e.value)
