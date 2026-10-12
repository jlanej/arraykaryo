"""QC flags the planted failures and nothing else; the GC-wave correction removes the planted wave; the sex calls recover
every planted complement and mosaic within 0.1 of the fraction, and name the swap (sections 4.2, 4.3)."""
import os

import numpy as np
import pytest

from arraykaryo.util import read_json, read_tsv


@pytest.fixture(scope="session")
def qc_out(out_dir, mock_dir):
    from arraykaryo.qc import run_qc
    run_qc(out_dir, sample_sheet=os.path.join(mock_dir, "sample_sheet.tsv"), threads=2, log_fn=lambda *a, **k: None)
    return out_dir


@pytest.fixture(scope="session")
def sex_out(qc_out):
    from arraykaryo.sex import run_sex
    run_sex(qc_out, log_fn=lambda *a, **k: None)
    return qc_out


def test_sample_qc_flags_the_planted_failures_and_nothing_else(qc_out, truth):
    rows = {r["sample"]: r for r in read_tsv(os.path.join(qc_out, "qc", "samples.tsv"))}
    assert len(rows) == truth["n_samples"]
    planted = {q["sample"]: q["reason"] for q in truth["qc_fail"]}
    for s, why in planted.items():
        assert why in rows[s]["flags"].split(","), (s, why, rows[s]["flags"])
        assert rows[s]["pass"] == "0"
        assert rows[s]["hard_fail"] == ("1" if why == "lrr_sd_fail" else "0")
    others = [s for s, r in rows.items() if r["flags"] != "NA" and s not in planted]
    assert others == [], others
    by = {s["sample"]: s for s in truth["samples"]}
    for s, r in rows.items():                                       # the measured noise is the planted noise
        if s not in planted:
            assert abs(float(r["lrr_sd"]) - by[s]["lrr_sd"]) < 0.02 and abs(float(r["baf_sd"]) - by[s]["baf_sd"]) < 0.012
    summ = read_json(os.path.join(qc_out, "qc", "summary.json"))
    assert summ["passing"] == truth["n_samples"] - len(planted) and summ["hard_fail"] == 2 and summ["use_batch"] is True


def test_gc_wave_correction_removes_the_wave(qc_out, truth):
    """On the null samples (no event: an event on a GC-rich chromosome is a real correlation, not a wave) the correlation of the
    corrected LRR with the 1-Mb GC falls under 0.05; on the planted strong waves it falls by an order of magnitude."""
    rows = {r["sample"]: r for r in read_tsv(os.path.join(qc_out, "qc", "samples.tsv"))}
    nulls = [rows[s] for s in truth["null_samples"]]
    after = np.array([abs(float(r["gc_corr_after"])) for r in nulls])
    before = np.array([abs(float(r["gc_corr_before"])) for r in nulls])
    assert after.max() < 0.05 and np.median(before) > np.median(after)
    by = {s["sample"]: s for s in truth["samples"]}
    strong = [r for r in rows.values() if "wave" in r["flags"]]
    assert all(abs(by[r["sample"]]["wave2"]) >= 0.05 for r in strong)       # only the planted strong waves are flagged
    assert all(abs(float(r["gc_corr_after"])) < 0.1 * abs(float(r["gc_corr_before"])) for r in strong)


def test_probe_qc_flags_the_planted_probes(qc_out, truth):
    rows = {r["id"]: r for r in read_tsv(os.path.join(qc_out, "qc", "probes.tsv"))}
    planted = truth["probe_fail"]
    hit = {}
    for pid, why in planted.items():
        r = rows[pid]
        got = r["reason"].split(",")
        hit[why] = hit.get(why, 0) + int(why in got)
        if why in got:
            assert r["keep_gt"] == "0"
            if why == "lrr_median":
                assert r["keep"] == "0"
    n = {}
    for why in planted.values():
        n[why] = n.get(why, 0) + 1
    for why in n:
        assert hit[why] >= 0.9 * n[why], (why, hit[why], n[why])
    # beyond the planted ones, a reason fires rarely (intensity-only and the top-1% noisy probes are by construction)
    by_reason = {}
    for pid, r in rows.items():
        for why in r["reason"].split(","):
            if why != "NA" and planted.get(pid) != why:
                by_reason[why] = by_reason.get(why, 0) + 1
    total = len(rows)
    for why, k in by_reason.items():
        if why not in ("intensity_only", "lrr_sd", "non_snv", "duplicate_position"):
            assert k < 0.01 * total, (why, k)
    assert all(r["par"] == "1" for pid, r in rows.items() if r["chrom"] == "chrX" and int(r["pos"]) <= 2394410)


def test_sex_recovers_every_planted_complement_and_mosaic(sex_out, truth):
    rows = {r["sample"]: r for r in read_tsv(os.path.join(sex_out, "sex", "sex.tsv"))}
    by = {s["sample"]: s for s in truth["samples"]}
    for e in truth["sex_events"]:
        r = rows[e["sample"]]
        if e["f"] >= 0.9:
            assert r["complement"] == by[e["sample"]]["complement"], (e, r)
            assert r["karyotype"] in ("47,XXY", "45,X")
        else:
            assert r["mosaic"] != "NA" and abs(float(r["f"]) - e["f"]) <= 0.1, (e, r["mosaic"], r["f"])
            assert ("loss of Y" in r["mosaic"]) == (e["chrom"] == "chrY")
            if e["chrom"] == "chrX":
                assert abs(float(r["f_baf"]) - e["f"]) <= 0.1 and r["karyotype"].startswith("mos 45,X[")
    planted = {e["sample"] for e in truth["sex_events"]}
    for s, r in rows.items():
        if s not in planted:
            assert r["complement"] == by[s]["complement"] and r["mosaic"] == "NA", (s, r)
            assert r["sex"] == by[s]["sex"]


def test_sex_check_names_the_swap(sex_out, truth):
    rows = {r["sample"]: r for r in read_tsv(os.path.join(sex_out, "sex", "sex.tsv"))}
    a, b = truth["swap"]
    assert "sample swap" in rows[a]["check"] and "sample swap" in rows[b]["check"]
    xxy = next(e["sample"] for e in truth["sex_events"] if "47,XXY" in e["label"])
    x0 = next(e["sample"] for e in truth["sex_events"] if "45,X" in e["label"])
    assert rows[xxy]["check"].startswith("XXY in a reported male") and rows[x0]["check"].startswith("X in a reported female")
    agree = [s for s, r in rows.items() if r["check"] == "agrees"]
    assert len(agree) == truth["n_samples"] - 4
    qc = {r["sample"]: r for r in read_tsv(os.path.join(sex_out, "qc", "samples.tsv"))}
    assert "sex_mismatch" in qc[a]["flags"] and "sex_mismatch" in qc[b]["flags"] and "sex_mismatch" not in qc[xxy]["flags"]
    assert os.path.exists(os.path.join(sex_out, "sex", "figures", "x_vs_y.png")) and os.path.exists(os.path.join(sex_out, "sex", "figures", "tracks_%s.png" % xxy))


def test_qc_and_sex_run_from_the_cli_and_are_idempotent(sex_out, capsys):
    from arraykaryo.cli import main
    assert main(["qc", "--out", sex_out]) == 0 and main(["sex", "--out", sex_out]) == 0
    err = capsys.readouterr().err
    assert "qc is done" in err and "sex is done" in err
