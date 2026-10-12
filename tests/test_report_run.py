"""The report and the digest contain what the tables contain; the CLI runs end to end on the mock (sections 4.6, 4.8)."""
import os
import re
import subprocess
import sys

import pytest

from arraykaryo.util import read_tsv


@pytest.fixture(scope="session")
def report_out(out_dir, mock_dir):
    from arraykaryo.events import run_events
    from arraykaryo.pedigree import run_pedigree
    from arraykaryo.qc import run_qc
    from arraykaryo.report import run_report
    from arraykaryo.sex import run_sex
    from arraykaryo.util import is_done
    quiet = lambda *a, **k: None  # noqa: E731
    if not is_done(out_dir, "qc"):
        run_qc(out_dir, sample_sheet=os.path.join(mock_dir, "sample_sheet.tsv"), threads=2, log_fn=quiet)
    if not is_done(out_dir, "sex"):
        run_sex(out_dir, log_fn=quiet)
    if not is_done(out_dir, "pedigree"):
        run_pedigree(out_dir, pedigree=os.path.join(mock_dir, "pedigree.tsv"), log_fn=quiet)
    if not is_done(out_dir, "events"):
        run_events(out_dir, threads=2, figures="flagged", log_fn=quiet)
    run_report(out_dir, log_fn=quiet)
    return out_dir


def test_report_contains_what_the_tables_contain(report_out, truth):
    html = open(os.path.join(report_out, "report", "index.html")).read()
    assert html.startswith("<!doctype html>") and "MOCK DATA" in html and "data:image/png;base64," in html
    ev = read_tsv(os.path.join(report_out, "events", "events.all.tsv"))
    for e in ev:
        assert e["sample"] in html and e["bands"] in html if e["bands"] != "NA" else True
    for s, r in {r["sample"]: r for r in read_tsv(os.path.join(report_out, "qc", "samples.tsv"))}.items():
        assert s in html
        if r["flags"] != "NA":
            assert os.path.exists(os.path.join(report_out, "report", "samples", s + ".html"))
    chk = read_tsv(os.path.join(report_out, "pedigree", "check.tsv"))
    for r in chk:
        if r["status"] != "agrees":
            assert r["status"] in html
    for s in {e["sample"] for e in ev}:
        page = open(os.path.join(report_out, "report", "samples", s + ".html")).read()
        assert s in page and "data:image/png;base64," in page and "Events" in page
    assert len(re.findall(r"<table", html)) >= 6 and 'id="events"' in html and 'id="samples"' in html
    guide = open(os.path.join(report_out, "report", "guide.html")).read()
    assert "Pattern cards" in guide and "f = 4d/(1-2d)" in guide and "kinship" in guide and "data:image/png;base64," in guide


def test_digest_ranks_the_largest_events(report_out, truth):
    rows = read_tsv(os.path.join(report_out, "report", "digest.tsv"))
    assert rows and rows[0]["span"] == "whole"
    tiers = [{"whole": 0, "p": 1, "q": 1, "stretch": 2}[r["span"]] for r in rows]
    assert tiers == sorted(tiers)
    html = open(os.path.join(report_out, "report", "digest.html")).read()
    for r in rows:
        assert r["sample"] in html and r["event"].split(" in ")[0] in html
    whole_const = [r for r in rows if r["span"] == "whole" and float(r["f"]) >= 0.9]
    assert whole_const and rows.index(whole_const[0]) < min(rows.index(r) for r in rows if r["span"] == "whole" and float(r["f"]) < 0.9)
    left = read_tsv(os.path.join(report_out, "report", "digest_left_out.tsv"))
    assert all(r["left_out"] for r in left)


@pytest.mark.slow
def test_cli_runs_end_to_end(mock_dir, tmp_path):
    out = str(tmp_path / "run")
    exe = os.path.join(os.path.dirname(sys.executable), "arraykaryo")
    cmd = [exe, "run", "--bcf", os.path.join(mock_dir, "mock.bcf"), "--out", out, "--threads", "2", "--sample-sheet", os.path.join(mock_dir, "sample_sheet.tsv"), "--pedigree", os.path.join(mock_dir, "pedigree.tsv")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    assert "all done" in r.stderr and "report/index.html" in r.stderr
    for p in ("cache/manifest.json", "qc/samples.tsv", "sex/sex.tsv", "pedigree/families.tsv", "pedigree/check.tsv", "events/events.all.tsv", "events/karyotypes.tsv", "report/index.html", "report/guide.html", "report/digest.html"):
        assert os.path.exists(os.path.join(out, p)), p
    r2 = subprocess.run(cmd, capture_output=True, text=True)                      # idempotent: every step is skipped
    assert r2.returncode == 0 and r2.stderr.count("is done") >= 5
