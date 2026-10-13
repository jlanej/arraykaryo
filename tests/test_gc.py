"""The GC of a probe's neighbourhood from a reference FASTA, and the QC step's behaviour when an export carries no INFO/GC."""
import os
import subprocess

import numpy as np
import pytest

from arraykaryo.util import read_json, read_tsv
from conftest import has_bcftools


def write_fasta(path, contigs, seed=3):
    """A random FASTA (no real sequence) with the given contig lengths, indexed."""
    import pysam
    rng = np.random.default_rng(seed)
    with open(path, "w") as fh:
        for name, n in contigs.items():
            seq = rng.choice(np.frombuffer(b"ACGTN", dtype=np.uint8), size=n, p=[0.3, 0.2, 0.2, 0.29, 0.01]).tobytes().decode()
            fh.write(">%s\n" % name)
            for i in range(0, n, 60):
                fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(path)
    return path


def test_gc_from_fasta_matches_a_direct_count(tmp_path):
    from arraykaryo.gc import gc_from_fasta, HALF
    fa = write_fasta(str(tmp_path / "ref.fa"), {"chr21": 30_000, "22": 5_000})
    import pysam
    seqs = {c: pysam.FastaFile(fa).fetch(c) for c in ("chr21", "22")}
    pos = np.array([1, 150, 201, 7_777, 29_900, 30_000, 30_001, 45_000])
    got = gc_from_fasta(fa, "chr21", pos, chunk=10_000)                        # several chunks, windows across chunk edges
    for p, g in zip(pos, got):
        if p > 30_000:
            assert np.isnan(g)
            continue
        w = seqs["chr21"][max(0, p - 1 - HALF):p + HALF].upper()
        n = sum(w.count(b) for b in "ACGT")
        assert abs(g - (w.count("G") + w.count("C")) / n) < 1e-6
    got22 = gc_from_fasta(fa, "chr22", np.array([2_500]))                      # the FASTA names it '22'
    w = seqs["22"][2_500 - 1 - HALF:2_500 + HALF]
    assert abs(got22[0] - (w.count("G") + w.count("C")) / sum(w.count(b) for b in "ACGT")) < 1e-6
    assert np.isnan(gc_from_fasta(fa, "chr1", np.array([100]))).all()           # a contig the FASTA lacks


@pytest.mark.skipif(not has_bcftools(), reason="bcftools strips INFO/GC from the mock")
def test_qc_without_info_gc_skips_the_wave_and_says_so_then_a_fasta_restores_it(tmp_path, mock_dir):
    from arraykaryo.extract import extract
    from arraykaryo.qc import run_qc, Correction
    from arraykaryo.cache import Cache
    nogc = str(tmp_path / "nogc.bcf")
    subprocess.run(["bcftools", "annotate", "-x", "INFO/GC", "-Ob", "-o", nogc, os.path.join(mock_dir, "mock.bcf")], check=True)
    subprocess.run(["bcftools", "index", nogc], check=True)
    out = str(tmp_path / "out")
    lines = []
    extract(nogc, os.path.join(out, "cache"), threads=2, route="bcftools", log_fn=lambda m, *a, **k: lines.append(m))
    C = Cache(os.path.join(out, "cache"))
    assert C.manifest["gc_source"] == "none" and not np.isfinite(C.probes("chr1").gc).any()
    assert any("no INFO/GC" in l for l in lines)
    run_qc(out, cache_dir=os.path.join(out, "cache"), sample_sheet=os.path.join(mock_dir, "sample_sheet.tsv"), threads=2, log_fn=lambda m, *a, **k: lines.append(m))
    s = read_json(os.path.join(out, "qc", "summary.json"))
    assert s["gc_source"] == "none" and s["gc"]["available"] is False
    assert any("NO GC" in l for l in lines)
    corr = Correction(out)
    assert not corr.gc_available and not corr.coef.any()
    L = corr.apply(C, "chr21", sample_idx=np.array([0, 1]))                    # the reference and offset still apply, nothing is NaN that was not
    assert L.shape == (2, C.n_probes["chr21"]) and np.isfinite(L).sum() > 0.9 * L.size
    q = read_tsv(os.path.join(out, "qc", "samples.tsv"))
    assert all(float(r["wf"]) >= 0 for r in q if r["wf"] != "NA")              # unsigned without GC
    # a FASTA for chr1's first 3 Mb (random sequence, the contig named '1'): the probes there get a GC, the others none, and qc fits the wave on them
    fa = write_fasta(str(tmp_path / "ref.fa"), {"1": 3_000_000})
    run_qc(out, cache_dir=os.path.join(out, "cache"), sample_sheet=os.path.join(mock_dir, "sample_sheet.tsv"), threads=2, fasta=fa, log_fn=lambda m, *a, **k: lines.append(m))
    C = Cache(os.path.join(out, "cache"))
    gc = C.probes("chr1").gc
    inside = C.probes("chr1").pos <= 3_000_000
    assert inside.sum() > 50
    assert np.isfinite(gc[inside]).all() and not np.isfinite(gc[~inside]).any() and not np.isfinite(C.probes("chr21").gc).any()
    assert C.manifest["gc_source"] == "INFO/GC + fasta" and C.manifest["gc_fasta"] == "ref.fa"
    tsv = [r for r in C.probe_table() if r["chrom"] == "chr1"]
    assert sum(r["gc"] != "NA" for r in tsv) == int(inside.sum())
    s = read_json(os.path.join(out, "qc", "summary.json"))
    assert s["gc"]["available"] is True and s["gc_source"] == "INFO/GC + fasta"
    assert Correction(out).gc_available
