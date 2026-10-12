"""The detection curve: mock cohorts at the array's own probe density with events planted over a grid of cell fractions and
sizes, run through the pipeline, and tabulated - detection, type, class and cell-fraction error per cell of the grid, with
the false-positive count on the null samples. Writes docs/calibration/calibration.tsv, false_positives.tsv, calibration.md
and calibration.png. No real data: everything is simulated.

  python scripts/calibrate.py --out docs/calibration [--density-kb 1.3] [--replicates 2] [--threads 8] [--work /tmp/calib]
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from arraykaryo.mock import build  # noqa: E402
from arraykaryo.util import read_tsv, write_tsv  # noqa: E402

FRACTIONS = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0]
SIZES_MB = [2, 5, 10, 25, "arm", "whole"]
KINDS = ["gain", "loss", "loh"]
CONTIGS = "chr1,chr7,chr11,chr12,chr15,chr21"


def grid(G):
    """Planted events: every (kind, f, size) on a chromosome chosen to fit, one per sample."""
    out = []
    spots = {"chr1": (20_000_000, 115_000_000), "chr7": (70_000_000, 150_000_000), "chr11": (60_000_000, 130_000_000), "chr12": (40_000_000, 125_000_000), "chr15": (25_000_000, 95_000_000)}
    k = 0
    for kind in KINDS:
        for f in FRACTIONS:
            if kind == "loh" and f == 1.0:
                continue
            for size in SIZES_MB:
                chrom = list(spots)[k % len(spots)]
                k += 1
                if size == "whole":
                    out.append(dict(chrom="chr21" if kind != "loh" else "chr15", start=1, end=None, f=f, kind=kind, cls="whole", size=size))
                elif size == "arm":
                    out.append(dict(chrom=chrom, start=G.p_end[chrom] + 1, end=None, f=f, kind=kind, cls="q", size=size))
                else:
                    a = spots[chrom][0] + (k * 7_000_003) % (spots[chrom][1] - spots[chrom][0] - size * 1_000_000)
                    out.append(dict(chrom=chrom, start=a, end=a + size * 1_000_000, f=f, kind=kind, cls="stretch", size=size))
    return out


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1) / float(min(a1 - a0 + 1, b1 - b0 + 1))


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit("failed: %s\n%s" % (" ".join(cmd), r.stderr[-2000:]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="docs/calibration")
    ap.add_argument("--work", default="calibration_work")
    ap.add_argument("--density-kb", type=float, default=1.3)
    ap.add_argument("--lrr-sd", type=float, default=0.15)
    ap.add_argument("--replicates", type=int, default=2)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--keep", action="store_true", help="keep the work directory")
    a = ap.parse_args()
    from arraykaryo.genome import load_genome
    G = load_genome("chm13")
    planted = grid(G)
    rows, fps = [], []
    t0 = time.time()
    exe = os.path.join(os.path.dirname(sys.executable), "arraykaryo")
    for rep in range(a.replicates):
        wd = os.path.join(a.work, "rep%d" % rep)
        md, od = os.path.join(wd, "mock"), os.path.join(wd, "out")
        shutil.rmtree(wd, ignore_errors=True)
        n = len(planted) + 40
        truth = build(md, n_samples=n, seed=100 + rep, contigs=CONTIGS, density_kb=a.density_kb, trios=0, duos=0, events=False, families=False, qc_failures=False,
                      extra_events=planted, lrr_sd=a.lrr_sd, log_fn=lambda *x, **k: None)
        run([exe, "extract", "--bcf", os.path.join(md, "mock.bcf"), "--out", od, "--threads", str(a.threads)])
        run([exe, "qc", "--out", od, "--sample-sheet", os.path.join(md, "sample_sheet.tsv"), "--threads", str(a.threads)])
        run([exe, "sex", "--out", od])
        run([exe, "pedigree", "--out", od])
        run([exe, "events", "--out", od, "--threads", str(a.threads), "--figures", "none"])
        ev = read_tsv(os.path.join(od, "events", "events.all.tsv"))
        by_sample = {}
        for e in ev:
            by_sample.setdefault(e["sample"], []).append(e)
        truth_by = {}
        for t in truth["events"]:
            truth_by[t["sample"]] = t
        spec = {}
        for t in truth["events"]:
            for x in planted:
                if x["chrom"] == t["chrom"] and x["kind"] == t["type"] and abs(x["f"] - t["f"]) < 1e-6 and abs(x["start"] - t["start"]) < 3_000_000:
                    spec[t["sample"]] = x
        for s, t in truth_by.items():
            x = spec.get(s, {})
            hits = [e for e in by_sample.get(s, []) if e["chrom"] == t["chrom"] and overlap(int(e["start"]), int(e["end"]), t["start"], t["end"]) >= 0.5]
            e = hits[0] if hits else None
            f_called = float(e["f"]) if e and e["f"] != "NA" else np.nan
            rows.append(dict(replicate=rep, kind=t["type"], f=t["f"], size=x.get("size", "?"), chrom=t["chrom"], n_probes=t["n_probes"], detected=int(bool(e)), type_called=e["type"] if e else "NA",
                             span_called=e["span"] if e else "NA", f_called=f_called, f_lrr=float(e["f_lrr"]) if e and e["f_lrr"] != "NA" else np.nan, f_baf=float(e["f_baf"]) if e and e["f_baf"] != "NA" else np.nan,
                             type_ok=int(bool(e) and e["type"] == t["type"]), class_ok=int(bool(e) and e["span"] == t["class"]), f_err=abs(f_called - t["f"]) if e else np.nan,
                             other_calls=len(by_sample.get(s, [])) - len(hits)))
        for s in truth["null_samples"]:
            for e in by_sample.get(s, []):
                fps.append(dict(replicate=rep, sample=s, chrom=e["chrom"], start=e["start"], end=e["end"], type=e["type"], span=e["span"], f=e["f"], lrr=e["lrr"], d=e["d"], note=e["note"]))
        n_null = len(truth["null_samples"])
        fps.append(dict(replicate=rep, sample="(null samples: %d)" % n_null, chrom="", start="", end="", type="", span="", f="", lrr="", d="", note=""))
        print("replicate %d: %d planted, %d detected, %d false-positive events on %d null samples (%.0f s)" % (rep, len(truth_by), sum(r["detected"] for r in rows if r["replicate"] == rep),
                                                                                                                 sum(1 for r in fps if r["replicate"] == rep and r["chrom"]), n_null, time.time() - t0))
        if not a.keep:
            shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(a.out, exist_ok=True)
    cols = ["replicate", "kind", "f", "size", "chrom", "n_probes", "detected", "type_called", "span_called", "f_called", "f_lrr", "f_baf", "type_ok", "class_ok", "f_err", "other_calls"]
    write_tsv(os.path.join(a.out, "calibration.tsv"), cols, rows, nd=4)
    write_tsv(os.path.join(a.out, "false_positives.tsv"), ["replicate", "sample", "chrom", "start", "end", "type", "span", "f", "lrr", "d", "note"], fps)
    # the tables
    md = ["# Calibration: the detection curve on the mock\n",
          "Mock cohorts at a probe spacing of %.1f kb (the example array's), per-sample LRR SD %.2f, BAF SD drawn in 0.03-0.05, %d replicate(s); every (kind, cell fraction, size) planted "
          "in its own sample on chr1, 7, 11, 12, 15 (stretches and arms) or chr21/chr15 (whole chromosomes); 40 null samples per replicate. Detection = an event of any type overlapping "
          "the planted span by half or more; the cell-fraction error is |called - planted| on the primary estimate.\n" % (a.density_kb, a.lrr_sd, a.replicates)]
    for kind in KINDS:
        md.append("\n## %s\n\nDetection rate (type right / class right), and the median cell-fraction error, by cell fraction (rows) and size (columns, Mb; arm = a whole q arm; whole = a whole chromosome).\n" % kind)
        md.append("| f | " + " | ".join(str(s) for s in SIZES_MB) + " |")
        md.append("|---|" + "---|" * len(SIZES_MB))
        for f in FRACTIONS:
            if kind == "loh" and f == 1.0:
                continue
            cells = []
            for size in SIZES_MB:
                rr = [r for r in rows if r["kind"] == kind and abs(r["f"] - f) < 1e-6 and str(r["size"]) == str(size)]
                if not rr:
                    cells.append("-")
                    continue
                det = np.mean([r["detected"] for r in rr])
                tok = np.mean([r["type_ok"] for r in rr])
                cok = np.mean([r["class_ok"] for r in rr])
                errs = [r["f_err"] for r in rr if np.isfinite(r["f_err"])]
                cells.append("%.0f%% (%.0f%%/%.0f%%) err %s" % (100 * det, 100 * tok, 100 * cok, ("%.2f" % np.median(errs)) if errs else "-"))
            md.append("| %.2f | " % f + " | ".join(cells) + " |")
    n_fp = sum(1 for x in fps if x["chrom"])
    n_null_samples = sum(int(x["sample"].split(": ")[1].rstrip(")")) for x in fps if x["sample"].startswith("(null"))
    fp_samples = len({(x["replicate"], x["sample"]) for x in fps if x["chrom"]})
    md.append("\n## False positives\n\n%d events on %d of %d null samples (the brief asks for no event in 19 of 20): see false_positives.tsv.\n" % (n_fp, fp_samples, n_null_samples))
    with open(os.path.join(a.out, "calibration.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    # the figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, kind in zip(axes, KINDS):
        fr = [f for f in FRACTIONS if not (kind == "loh" and f == 1.0)]
        M = np.full((len(fr), len(SIZES_MB)), np.nan)
        for i, f in enumerate(fr):
            for j, size in enumerate(SIZES_MB):
                rr = [r for r in rows if r["kind"] == kind and abs(r["f"] - f) < 1e-6 and str(r["size"]) == str(size)]
                if rr:
                    M[i, j] = np.mean([r["detected"] * r["type_ok"] for r in rr])
        im = ax.imshow(M, vmin=0, vmax=1, cmap="viridis", aspect="auto", origin="lower")
        ax.set_xticks(range(len(SIZES_MB))), ax.set_xticklabels([str(s) for s in SIZES_MB])
        ax.set_yticks(range(len(fr))), ax.set_yticklabels(["%.2f" % f for f in fr])
        ax.set_xlabel("size (Mb; arm; whole chromosome)"), ax.set_ylabel("cell fraction")
        ax.set_title("%s: detected with the right type" % kind)
        for i in range(len(fr)):
            for j in range(len(SIZES_MB)):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, "%.0f" % (100 * M[i, j]), ha="center", va="center", fontsize=7, color="white" if M[i, j] < 0.6 else "black")
    fig.colorbar(im, ax=axes, fraction=0.02, label="detection rate")
    fig.suptitle("arraykaryo detection on the mock at %.1f-kb probe spacing, LRR SD %.2f (%d replicates; %d false-positive events on %d null samples)" % (a.density_kb, a.lrr_sd, a.replicates, n_fp, n_null_samples), fontsize=9)
    fig.savefig(os.path.join(a.out, "calibration.png"), dpi=150, bbox_inches="tight")
    print("-> %s (%.0f s)" % (a.out, time.time() - t0))


if __name__ == "__main__":
    main()
