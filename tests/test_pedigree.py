"""The pedigree recovers every planted family and flags the non-paternity, the swap and the duplicate, with no spurious
first-degree pair; the kinship matches PLINK 2's KING-robust when plink2 is on PATH (section 4.4)."""
import os
import shutil
import subprocess

import numpy as np
import pytest

from arraykaryo.util import read_json, read_tsv


@pytest.fixture(scope="session")
def ped_out(out_dir, mock_dir):
    from arraykaryo.qc import run_qc
    from arraykaryo.sex import run_sex
    from arraykaryo.pedigree import run_pedigree
    from arraykaryo.util import is_done
    quiet = lambda *a, **k: None  # noqa: E731
    if not is_done(out_dir, "qc"):
        run_qc(out_dir, sample_sheet=os.path.join(mock_dir, "sample_sheet.tsv"), threads=2, log_fn=quiet)
    if not is_done(out_dir, "sex"):
        run_sex(out_dir, log_fn=quiet)
    run_pedigree(out_dir, pedigree=os.path.join(mock_dir, "pedigree.tsv"), sample_sheet=os.path.join(mock_dir, "sample_sheet.tsv"), log_fn=quiet)
    return out_dir


def test_families_recovered_and_no_spurious_first_degree_pair(ped_out, truth):
    fams = read_tsv(os.path.join(ped_out, "pedigree", "families.tsv"))
    swap = dict(zip(truth["swap"], truth["swap"][::-1]))
    nonpat = truth["non_paternity"]
    trios = {(f["kid"], f["dad"], f["mom"]) for f in fams if f["kind"] == "trio"}
    duos = {(f["kid"], f["dad"], f["mom"]) for f in fams if f["kind"] == "duo"}
    for f in truth["families"]:
        dad = swap.get(f["dad"], f["dad"])                               # the true father's data sits under the swapped name
        if f["kind"] == "trio":
            if f["kid"] == nonpat["kid"]:
                assert (f["kid"], "NA", f["mom"]) in duos, (f, duos)      # the mother remains as a duo
            else:
                assert (f["kid"], dad, f["mom"]) in trios, (f, trios)
        else:
            assert (f["kid"], f["dad"] or "NA", f["mom"] or "NA") in duos, (f, duos)
    assert any("Mendelian errors concentrated on chr15" in f["note"] for f in fams if f["kind"] == "trio")       # the isodisomy child
    rel = read_tsv(os.path.join(ped_out, "pedigree", "relatedness.tsv"))
    first = [(r["sample_a"], r["sample_b"], r["class"]) for r in rel if r["class"] in ("parent-offspring", "full sibling", "duplicate or monozygotic")]
    expected_pairs = set()
    for f in truth["families"]:
        dad = swap.get(f["dad"], f["dad"])
        for par in (dad if f["kid"] != nonpat["kid"] else "", f["mom"]):
            if par:
                expected_pairs.add(frozenset((f["kid"], par)))
    expected_pairs.add(frozenset(truth["duplicate"]))
    got = {frozenset((a, b)) for a, b, _ in first}
    assert got == expected_pairs, (got ^ expected_pairs)
    assert all(c == "duplicate or monozygotic" for a, b, c in first if frozenset((a, b)) == frozenset(truth["duplicate"]))
    for r in rel:
        assert float(r["ibs0"]) < 0.03 if r["class"] != "duplicate or monozygotic" else float(r["kinship"]) > 0.45


def test_check_against_the_given_pedigree(ped_out, truth):
    chk = read_tsv(os.path.join(ped_out, "pedigree", "check.tsv"))
    a, b = truth["swap"]
    nonpat = truth["non_paternity"]
    rows = {(r["kid"], r["relative"], r["status"]): r for r in chk}
    swap_row = next(r for r in chk if r["relative"] == a and r["status"] == "claimed father unrelated")
    assert "sample swap" in swap_row["note"]
    np_row = next(r for r in chk if r["kid"] == nonpat["kid"] and r["relative"] == nonpat["claimed_dad"])
    assert np_row["status"] == "claimed father unrelated" and "non-paternity" in np_row["note"]
    assert any(r["status"] == "undisclosed parent-offspring pair" and {r["kid"], r["relative"]} & {b} for r in chk)
    assert any(r["status"] == "undisclosed duplicate" and {r["kid"], r["relative"]} == set(truth["duplicate"]) for r in chk)
    assert any(r["status"] == "sex contradiction" and r["relative"] == a for r in chk)
    agrees = [r for r in chk if r["status"] == "agrees"]
    n_claimed = sum((1 if f["dad"] else 0) + (1 if f["mom"] else 0) for f in truth["families"])
    assert len(agrees) == n_claimed - 2                                   # the swapped father and the non-father
    qc = {r["sample"]: r for r in read_tsv(os.path.join(ped_out, "qc", "samples.tsv"))}
    assert all("duplicate" in qc[s]["flags"] for s in truth["duplicate"]) and "unexpected_relative" in qc[b]["flags"]
    summ = read_json(os.path.join(ped_out, "pedigree", "summary.json"))
    assert summ["families"]["trio"] == 3 and summ["hwe_failed_probes"] >= 30 and os.path.exists(os.path.join(ped_out, "pedigree", "figures", "kinship_vs_ibs0.png"))
    inferred = read_tsv(os.path.join(ped_out, "pedigree", "pedigree.inferred.tsv"))
    assert len(inferred) == 6 and all(r["kid_sex"] in ("1", "2") for r in inferred)


def test_thinning_does_not_move_the_kinship(ped_out):
    """Halving the thinning distance leaves the related pairs' kinship within 0.02."""
    from arraykaryo.cache import Cache
    from arraykaryo.pedigree import PARAMS, kinship_matrices, king_robust, select_probes
    from arraykaryo.qc import probe_masks
    C = Cache(os.path.join(ped_out, "cache"))
    masks = probe_masks(ped_out, C)
    ks = []
    for thin in (50_000, 25_000):
        p = dict(PARAMS, thin_bp=thin)
        hh, ibs0, both, het_i, e, n = kinship_matrices(C, select_probes(C, masks, p))
        ks.append(king_robust(hh, ibs0, het_i))
    rel = ks[0] >= 0.0442
    assert rel.any() and np.nanmax(np.abs(ks[0][rel] - ks[1][rel])) < 0.02


@pytest.mark.skipif(shutil.which("plink2") is None, reason="plink2 not on PATH")
def test_kinship_matches_plink2(ped_out, mock_dir, tmp_path):
    from arraykaryo.cache import Cache
    from arraykaryo.pedigree import PARAMS, select_probes
    from arraykaryo.qc import probe_masks
    C = Cache(os.path.join(ped_out, "cache"))
    sel = select_probes(C, probe_masks(ped_out, C), PARAMS)
    ids = [i for c, idx in sel for i in C.probes(c).ids[idx].tolist()]
    (tmp_path / "probes.txt").write_text("\n".join(ids) + "\n")
    r = subprocess.run(["plink2", "--vcf", os.path.join(mock_dir, "mock.vcf.gz"), "--autosome", "--extract", str(tmp_path / "probes.txt"), "--make-king-table", "--king-table-filter", "0.02",
                        "--out", str(tmp_path / "king"), "--silent"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    z = np.load(os.path.join(ped_out, "pedigree", "kinship.npz"))
    names = list(z["samples"])
    rows = [l.lstrip("#").rstrip("\n").split("\t") for l in open(tmp_path / "king.kin0")]
    ix = {h: i for i, h in enumerate(rows[0])}
    n = 0
    for r in rows[1:]:
        a, b = names.index(r[ix["IID1"]]), names.index(r[ix["IID2"]])
        assert abs(float(r[ix["KINSHIP"]]) - z["kinship"][a, b]) < 0.002 and abs(float(r[ix["IBS0"]]) - z["ibs0"][a, b]) < 1e-4
        n += 1
    assert n >= 9
