"""The mock writes what the tests need (section 4.7) and the extraction round-trips it exactly (section 4.1)."""
import math
import os
import shutil

import numpy as np
import pysam
import pytest

from arraykaryo.cli import main
from arraykaryo.genome import load_genome
from arraykaryo.util import read_pedigree, read_sample_sheet


def test_truth_lists_every_planted_role(truth, mock_dir):
    names = {s["sample"] for s in truth["samples"]}
    assert len(names) == truth["n_samples"]
    types = {(e["type"], e["class"]) for e in truth["events"]}
    assert {("gain", "whole"), ("loss", "q"), ("loh", "p"), ("loss", "stretch"), ("roh", "stretch"), ("roh", "whole")} <= types
    labels = {e["label"] for e in truth["sex_events"]}
    assert any("47,XXY" in l for l in labels) and any("45,X" in l for l in labels) and any("loss of Y" in l for l in labels) and any("loss of an X" in l for l in labels)
    assert len(truth["swap"]) == 2 and truth["non_paternity"]["kid"] and len(truth["duplicate"]) == 2
    assert {q["reason"] for q in truth["qc_fail"]} == {"lrr_sd_fail", "wave", "baf_sd", "call_rate"} and len(truth["qc_fail"]) == 6
    assert len(truth["null_samples"]) >= 20 and set(truth["null_samples"]) <= names
    assert all(e["sample"] in names for e in truth["events"] + truth["sex_events"])
    assert set(truth["probe_fail"].values()) == {"non_snv", "call_rate", "lrr_median", "baf_shift", "gentrain", "cluster_sep", "hwe"}
    assert os.path.exists(os.path.join(mock_dir, "MOCK_DATA.txt"))


def test_sample_sheet_and_pedigree_carry_the_lab_view(truth, mock_dir):
    sheet, cols = read_sample_sheet(os.path.join(mock_dir, "sample_sheet.tsv"))
    assert {"sample", "reported_sex", "batch"} <= set(cols) and len(sheet) == truth["n_samples"]
    a, b = truth["swap"]
    by = {s["sample"]: s for s in truth["samples"]}
    assert sheet[a]["reported_sex"] != by[a]["sex"] and sheet[b]["reported_sex"] != by[b]["sex"]        # the swap: the sheet names the other person
    assert sum(1 for s in truth["samples"] if sheet[s["sample"]]["reported_sex"] == s["sex"]) == truth["n_samples"] - 2
    rows, sex_of, _ = read_pedigree(os.path.join(mock_dir, "pedigree.tsv"))
    assert len(rows) == len(truth["families"]) and sum(1 for r in rows if r["dad"] and r["mom"]) == sum(1 for f in truth["families"] if f["kind"] == "trio")
    assert sex_of[truth["non_paternity"]["claimed_dad"]] == "M"


def test_export_conventions(mock_dir):
    """Diploid GT everywhere, nan for missing, BAF unclipped, intensity-only probes with ALT '.', PAR probes, a female's Y as no-calls."""
    from arraykaryo.inspect import inspect
    d = inspect(os.path.join(mock_dir, "mock.vcf.gz"), records=3000, samples=None)
    s = d["sampled"]
    assert d["assembly"] == "chm13" and d["required_missing"] == [] and s["baf_max"] > 1.0 and s["baf_min"] < 0.0 and s["intensity_only"] > 0 and s["alt_missing"] == s["intensity_only"]
    assert s["lrr_missing"] > 0 and s["gt_nocall"] > 0 and s["samples_read"] == d["samples"]


@pytest.mark.parametrize("route", ["bcftools", "cyvcf2"])
def test_extract_round_trips_the_export_exactly(mock_dir, tmp_path, route):
    if route == "cyvcf2":
        pytest.importorskip("cyvcf2")
    from arraykaryo.cache import Cache
    from arraykaryo.extract import extract
    cache_dir = str(tmp_path / ("cache_" + route))
    extract(os.path.join(mock_dir, "mock.bcf"), cache_dir, threads=1, route=route, chroms=["chr21", "chrX", "chrY"], log_fn=lambda *a, **k: None)
    C = Cache(cache_dir)
    assert C.chroms == ["chr21", "chrX", "chrY"] and C.manifest["route"] == route and C.manifest["assembly"] == "chm13"
    vf = pysam.VariantFile(os.path.join(mock_dir, "mock.vcf.gz"))
    names = C.samples[:3] + C.samples[-2:]
    vf.subset_samples(names)
    idx = [C.sample_index(s) for s in names]
    n = 0
    for chrom in C.chroms:
        gt, lrr, baf, pr = C.gt(chrom), C.lrr(chrom), C.baf(chrom), C.probes(chrom)
        assert gt.shape == (C.n_samples, pr.n) == lrr.shape == baf.shape and gt.dtype == np.int8 and lrr.dtype == np.float32 and baf.dtype == np.float16
        for i, rec in enumerate(vf.fetch(chrom)):
            assert rec.pos == pr.pos[i] and pr.ids[i] == rec.id
            for j, s in zip(idx, names):
                v = rec.samples[s]
                g = v["GT"]
                exp = -1 if (g is None or any(a is None for a in g)) else int(sum(min(a, 1) for a in g))
                assert gt[j, i] == exp
                L, B = v["LRR"], v["BAF"]
                assert (math.isnan(L) and np.isnan(lrr[j, i])) or np.float32(L) == lrr[j, i]
                assert (math.isnan(B) and np.isnan(baf[j, i])) or np.float16(B) == baf[j, i]
                n += 1
    assert n > 10000
    px = C.probes("chrX")
    assert px.par.sum() > 0 and px.par[px.pos <= 2394410].all() and not px.par[(px.pos > 2394410) & (px.pos < 153925835)].any()
    p21 = C.probes("chr21")
    assert p21.intensity_only.sum() > 0 and (~p21.snv[p21.intensity_only]).all() and np.isfinite(p21.gc).all() and np.isfinite(p21.gentrain).all()


def test_extract_is_resumable_and_idempotent(mock_dir, tmp_path, capsys):
    from arraykaryo.extract import extract
    cache_dir = str(tmp_path / "cache_resume")
    msgs = []
    extract(os.path.join(mock_dir, "mock.bcf"), cache_dir, threads=1, route="bcftools", chroms=["chr21", "chrY"], log_fn=lambda m, s="": msgs.append(m))
    os.remove(os.path.join(cache_dir, "done", "chrY.json"))
    msgs.clear()
    extract(os.path.join(mock_dir, "mock.bcf"), cache_dir, threads=1, route="bcftools", chroms=["chr21", "chrY"], log_fn=lambda m, s="": msgs.append(m))
    assert any("1 to do" in m for m in msgs) and any(m.startswith("chrY:") for m in msgs) and not any(m.startswith("chr21:") for m in msgs)
    assert main(["extract", "--bcf", os.path.join(mock_dir, "mock.bcf"), "--out", str(tmp_path / "o"), "--cache", cache_dir, "--chroms", "chr21,chrY"]) == 0


def test_genome_tables_agree_with_the_assemblies():
    for name, c1 in (("grch37", 249250621), ("grch38", 248956422), ("chm13", 248387328)):
        G = load_genome(name)
        assert G.length["chr1"] == c1 and len(G.chroms) == 24 and all(G.p_end[c] > 0 for c in G.chroms)
        assert G.arm("chr1", 1000) == "p" and G.arm("chr1", G.length["chr1"]) == "q" and G.band("chr1", 1000).startswith("p36")
        assert G.is_par("chrX", G.par["chrX"][0][1]) and not G.is_par("chrX", G.par["chrX"][0][1] + 1) and G.is_par("chrY", G.par["chrY"][0][0])
        names = {a[0] for a in G.artefacts()}
        assert {"MHC", "8p23.1 inversion", "17q21.31 inversion", "pericentromeric", "subtelomeric", "acrocentric p arm"} <= names
    with pytest.raises(ValueError):
        load_genome("hg17")


def test_mock_on_grch38(tmp_path):
    from arraykaryo.mock import build
    from arraykaryo.inspect import inspect
    d = str(tmp_path / "m38")
    t = build(d, n_samples=30, seed=3, genome="grch38", contigs="chr21,chrX,chrY", trios=0, duos=0, log_fn=lambda *a, **k: None)
    assert t["genome"] == "grch38" and inspect(os.path.join(d, "mock.vcf.gz"), records=10)["assembly"] == "grch38"
