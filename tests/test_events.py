"""Every planted event is found with its type, class and cell fraction within 0.08, and the ROH; null samples receive no
event in 19 of 20; the karyotype strings match the truth; figures exist for samples with an event (section 4.5)."""
import os

import pytest

from arraykaryo.util import read_json, read_tsv

TYPE_MAP = {"gain": "gain", "loss": "loss", "loh": "loh", "roh": "roh"}


@pytest.fixture(scope="session")
def events_out(out_dir, mock_dir):
    from arraykaryo.events import run_events
    from arraykaryo.pedigree import run_pedigree
    from arraykaryo.qc import run_qc
    from arraykaryo.sex import run_sex
    from arraykaryo.util import is_done
    quiet = lambda *a, **k: None  # noqa: E731
    if not is_done(out_dir, "qc"):
        run_qc(out_dir, sample_sheet=os.path.join(mock_dir, "sample_sheet.tsv"), threads=2, log_fn=quiet)
    if not is_done(out_dir, "sex"):
        run_sex(out_dir, log_fn=quiet)
    if not is_done(out_dir, "pedigree"):
        run_pedigree(out_dir, pedigree=os.path.join(mock_dir, "pedigree.tsv"), log_fn=quiet)
    run_events(out_dir, threads=2, figures="flagged", log_fn=quiet)
    return out_dir


def _overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1) / float(min(a1 - a0 + 1, b1 - b0 + 1))


def test_every_planted_event_is_found(events_out, truth):
    ev = read_tsv(os.path.join(events_out, "events", "events.all.tsv"))
    for t in truth["events"]:
        cands = [e for e in ev if e["sample"] == t["sample"] and e["chrom"] == t["chrom"] and _overlap(int(e["start"]), int(e["end"]), t["start"], t["end"]) >= 0.8]
        assert cands, ("not found", t)
        e = cands[0]
        assert e["type"] == TYPE_MAP[t["type"]], (t, e)
        assert e["span"] == t["class"], (t["label"], e["span"])
        if t["f"] < 0.9:
            assert abs(float(e["f"]) - t["f"]) <= 0.08, (t["label"], e["f"])
            assert e["constitutional"] == "0" and abs(float(e["f_baf"]) - t["f"]) <= 0.08
        else:
            assert e["constitutional"] == "1"
        assert e["bands"] == t["bands"] or e["bands"][:3] == t["bands"][:3]
    for t in truth["sex_events"]:
        cands = [e for e in ev if e["sample"] == t["sample"] and e["chrom"] == t["chrom"] and e["source"] == "sex"]
        assert cands and cands[0]["type"] == t["type"] and abs(float(cands[0]["f"]) - t["f"]) <= 0.1, (t, cands)


def test_null_samples_receive_no_event(events_out, truth):
    ev = read_tsv(os.path.join(events_out, "events", "events.all.tsv"))
    nulls = set(truth["null_samples"])
    with_event = {e["sample"] for e in ev if e["sample"] in nulls}
    assert len(with_event) <= len(nulls) // 20, sorted(with_event)
    planted = {e["sample"] for e in truth["events"]} | {e["sample"] for e in truth["sex_events"]}
    extra = [e for e in ev if e["sample"] in planted and not any(e["chrom"] == t["chrom"] for t in truth["events"] + truth["sex_events"] if t["sample"] == e["sample"])]
    assert len(extra) <= 1, extra                                        # at most one unplanted call on a planted sample


def test_roh_and_f_roh(events_out, truth):
    roh = read_tsv(os.path.join(events_out, "events", "roh.all.tsv"))
    kar = {r["sample"]: r for r in read_tsv(os.path.join(events_out, "events", "karyotypes.tsv"))}
    planted = [t for t in truth["events"] if t["type"] == "roh"]
    for t in planted:
        runs = [r for r in roh if r["sample"] == t["sample"] and r["chrom"] == t["chrom"] and _overlap(int(r["start"]), int(r["end"]), t["start"], t["end"]) >= 0.9]
        assert runs and runs[0]["is_event"] == "1", t["label"]
        assert float(kar[t["sample"]]["f_roh"]) > 0.01
    whole = [t for t in planted if t["class"] == "whole"]
    assert whole and float(kar[whole[0]["sample"]]["f_roh"]) > 0.05
    for s in truth["null_samples"]:
        assert float(kar[s]["f_roh"]) < 0.01


def test_karyotype_strings(events_out, truth):
    kar = {r["sample"]: r["karyotype"] for r in read_tsv(os.path.join(events_out, "events", "karyotypes.tsv"))}
    by = {t["sample"]: t for t in truth["events"]}
    sx = {t["sample"]: t for t in truth["sex_events"]}
    for s, t in by.items():
        k = kar[s]
        n = t["chrom"][3:]
        if t["type"] == "gain" and t["class"] == "whole" and t["f"] >= 0.9:
            assert k.startswith("47,X") and ",+%s" % n in k, k
        elif t["type"] == "gain" and t["class"] == "whole":
            assert k.startswith("mos 47,X") and "+%s[" % n in k, k
        elif t["type"] == "loss" and t["class"] in ("p", "q"):
            assert "del(%s)(" % n in k and "[0." in k, k
        elif t["type"] == "loh":
            assert "loh(%s)(" % n in k and "[0." in k, k
        elif t["type"] == "loss" and t["class"] == "stretch":
            assert "del(%s)(" % n in k and "mos" not in k, k
        elif t["type"] == "roh" and t["class"] == "whole":
            assert "upd(%s)(iso)?" % n in k, k
        elif t["type"] == "roh":
            assert "roh(%s)(" % n in k, k
    for s, t in sx.items():
        k = kar[s]
        if "47,XXY" in t["label"]:
            assert k == "47,XXY", k
        elif "45,X" in t["label"] and t["f"] >= 0.9:
            assert k == "45,X", k
        elif t["chrom"] == "chrY":
            assert k.startswith("mos 45,X[0.") and k.endswith("/46,XY"), k
        else:
            assert k.startswith("mos 45,X[0.") and k.endswith("/46,XX"), k
    for s in truth["null_samples"]:
        assert kar[s] in ("46,XX", "46,XY"), (s, kar[s])


def test_figures_tables_and_summary(events_out, truth):
    summ = read_json(os.path.join(events_out, "events", "summary.json"))
    assert summ["events"] >= len(truth["events"]) + len(truth["sex_events"]) and "the median LRR of 1 constitutional trisomy" in summ["calibration"]["three_copies_source"]
    assert abs(summ["calibration"]["one_copy_lrr"] - truth["levels"]["one_copy"]) < 0.06 and abs(summ["calibration"]["three_copies_lrr"] - truth["levels"]["three_copies"]) < 0.06
    for t in truth["events"][:3]:
        d = os.path.join(events_out, "events", "figures", t["sample"])
        assert os.path.exists(os.path.join(d, "genome.png")) and os.path.exists(os.path.join(d, "genome.txt")) and os.path.exists(os.path.join(d, "chr_%s.png" % t["chrom"]))
    qc = {r["sample"]: r for r in read_tsv(os.path.join(events_out, "qc", "samples.tsv"))}
    assert all("n_events" in r for r in qc.values()) and not any("event_count" in r["flags"] for r in qc.values())
