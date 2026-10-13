#!/usr/bin/env python3
"""Build docs/index.html, the showcase of arraykaryo on the public 1000 Genomes example, from a finished run (GitHub Pages
serves docs/). Every number and every drawing on the page comes from the run's tables; the dense figures are copied into
docs/assets. No path is written anywhere: the run directory is an argument.

    python scripts/build_site.py RUN_DIR docs
"""
import csv
import html
import json
import os
import shutil
import sys
from collections import Counter, defaultdict
from string import Template

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from arraykaryo import __version__                       # noqa: E402
from arraykaryo.genome import load_genome                # noqa: E402

GITHUB = "https://github.com/jlanej/arraykaryo"
COL = dict(gain="#D55E00", loss="#0072B2", loh="#CC79A7", roh="#E69F00", XX="#E69F00", XY="#0072B2", XXY="#CC79A7", X="#009E73")
STAIN = dict(gneg="#f4f1ea", gpos25="#cfc9bd", gpos50="#a8a195", gpos75="#7d7669", gpos100="#4a443b", acen="#c0392b", gvar="#d9d3c7", stalk="#bdb5a8")


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fmt_int(n):
    return "{:,}".format(int(n))


def esc(s):
    return html.escape(str(s), quote=True)


def pct(f):
    return "every cell" if f >= 0.995 else "%d%% of cells" % round(100 * f)


# ------------------------------------------------------------------------------------------------------------ the data

def load(run):
    d = dict(run=run)
    d["manifest"] = json.load(open(os.path.join(run, "cache", "manifest.json")))
    d["qc"] = json.load(open(os.path.join(run, "qc", "summary.json")))
    d["sexs"] = json.load(open(os.path.join(run, "sex", "summary.json")))
    d["ped"] = json.load(open(os.path.join(run, "pedigree", "summary.json")))
    d["ev"] = json.load(open(os.path.join(run, "events", "summary.json")))
    d["sex"] = tsv(os.path.join(run, "sex", "sex.tsv"))
    d["events"] = tsv(os.path.join(run, "events", "events.all.tsv"))
    d["kar"] = {r["sample"]: r for r in tsv(os.path.join(run, "events", "karyotypes.tsv"))}
    d["digest"] = tsv(os.path.join(run, "report", "digest.tsv"))
    d["G"] = load_genome(d["manifest"]["assembly"])
    return d


def tagged(e):
    return e["artefact"] != "NA"


def event_words(e):
    kind = dict(gain="gain", loss="loss", loh="copy-neutral LOH", roh="run of homozygosity")[e["type"]]
    where = e["chrom"][3:] + ("" if e["span"] == "whole" else e["bands"].replace("chr", ""))
    f = num(e["f"])
    cells = "" if e["type"] == "roh" else (", " + pct(f) if f == f else "")
    return "%s: %s %s (%.1f Mb%s)" % (e["sample"], kind, where, num(e["size_bp"]) / 1e6, cells)


# --------------------------------------------------------------------------------------------------------- the drawings

def karyogram_svg(d):
    """Every untagged event of the run on the assembly's ideograms: copy changes above each chromosome, runs of homozygosity
    below; the sex step's findings on X and Y. Hover gives the sample."""
    G = d["G"]
    chroms = [c for c in G.chroms if c != "chrM"]
    W, x0, right = 1000, 54, 12
    scale = (W - x0 - right) / max(G.length[c] for c in chroms)
    above_lanes, below_lanes, lane_h = 9, 4, 2.6
    row_h = above_lanes * lane_h + 12 + below_lanes * lane_h + 12
    H = 8 + row_h * len(chroms)
    ev = [e for e in d["events"] if not tagged(e) and e["chrom"] in G.length and e["source"] != "sex"]   # the sex step's whole-X and whole-Y findings are on the sex plot
    by = defaultdict(list)
    for e in ev:
        by[e["chrom"]].append(e)
    out = ['<svg class="karyo" viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="every untagged event on the chromosome ideograms">' % (W, H)]
    out.append("<defs>")
    for i, c in enumerate(chroms):
        y = 8 + i * row_h + above_lanes * lane_h + 1
        out.append('<clipPath id="cp%s"><rect x="%.1f" y="%.1f" width="%.1f" height="10" rx="5"/></clipPath>' % (c, x0, y, G.length[c] * scale))
    out.append("</defs>")
    for i, c in enumerate(chroms):
        ytop = 8 + i * row_h
        y = ytop + above_lanes * lane_h + 1
        out.append('<text x="%d" y="%.1f" class="klabel">%s</text>' % (x0 - 8, y + 9, c[3:]))
        out.append('<g clip-path="url(#cp%s)">' % c)
        for (s, e, band, stain) in G.bands[c]:
            if stain == "acen":
                out.append('<rect x="%.2f" y="%.1f" width="%.2f" height="6" fill="%s"/>' % (x0 + s * scale, y + 2, max((e - s) * scale, 0.6), STAIN["acen"]))
            else:
                out.append('<rect x="%.2f" y="%.1f" width="%.2f" height="10" fill="%s"><title>%s%s</title></rect>' % (x0 + s * scale, y, max((e - s) * scale, 0.5), STAIN.get(stain, "#eee"), c[3:], band))
        out.append("</g>")
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="10" rx="5" fill="none" stroke="#8a8378" stroke-width="0.6"/>' % (x0, y, G.length[c] * scale))
        # lanes
        for side, kinds, nl in (("above", ("gain", "loss", "loh"), above_lanes), ("below", ("roh",), below_lanes)):
            rows = sorted([e for e in by[c] if e["type"] in kinds], key=lambda e: (-num(e["size_bp"]), num(e["start"])))   # the longest nearest the ideogram
            ends = []
            for e in rows:
                s, t = num(e["start"]), num(e["end"])
                lane = next((k for k, taken in enumerate(ends) if all(t + 300_000 < a or s > b + 300_000 for a, b in taken)), None)
                if lane is None:
                    ends.append([(s, t)])
                    lane = len(ends) - 1
                else:
                    ends[lane].append((s, t))
                k = lane % nl
                yy = (y - 2 - (k + 1) * lane_h) if side == "above" else (y + 12 + 1 + k * lane_h)
                x1, w = x0 + s * scale, max((t - s) * scale, 1.4)
                out.append('<rect class="ev ev-%s" x="%.2f" y="%.2f" width="%.2f" height="%.1f" fill="%s" fill-opacity="0.9" data-t="%s"><title>%s</title></rect>' % (
                    e["type"], x1, yy, w, lane_h - 0.7, COL[e["type"]], esc(event_words(e)), esc(event_words(e))))
    out.append("</svg>")
    return "\n".join(out), len(ev)


def sex_svg(d):
    """X level against Y level of every sample, from sex.tsv, coloured by the complement; mosaics ringed; the odd ones named."""
    rows = [r for r in d["sex"] if num(r["x_lrr"]) == num(r["x_lrr"]) and num(r["y_lrr"]) == num(r["y_lrr"])]
    W, H, ml, mr, mt, mb = 640, 440, 54, 16, 14, 44
    xs = [num(r["x_lrr"]) for r in rows]
    ys = [num(r["y_lrr"]) for r in rows]
    xmin, xmax = min(xs) - 0.04, max(xs) + 0.04
    ymin, ymax = min(ys) - 0.25, max(ys) + 0.35
    sx = lambda v: ml + (v - xmin) / (xmax - xmin) * (W - ml - mr)
    sy = lambda v: mt + (ymax - v) / (ymax - ymin) * (H - mt - mb)
    o = ['<svg class="sexplot" viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="X level against Y level of every sample">' % (W, H)]
    o.append('<rect x="%d" y="%d" width="%d" height="%d" class="plotbg"/>' % (ml, mt, W - ml - mr, H - mt - mb))
    xm, ym = d["sexs"]["x_modes"], d["sexs"]["y_modes"]
    for v, lab in ((xm["one_copy"], "one X"), (xm["two_copies"], "two X's")):
        o.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" class="mode"/><text x="%.1f" y="%d" class="modelab" text-anchor="middle">%s</text>' % (sx(v), sx(v), mt, H - mb, sx(v), H - mb + 14, lab))
    for v, lab in ((ym["present"], "Y present"), (ym["absent"], "Y absent")):
        o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" class="mode"/><text x="%d" y="%.1f" class="modelab">%s</text>' % (ml, W - mr, sy(v), sy(v), ml + 6, sy(v) - 4, lab))
    for t in (-0.6, -0.4, -0.2, 0.0, 0.2):
        if xmin < t < xmax:
            o.append('<text x="%.1f" y="%d" class="tick" text-anchor="middle">%g</text>' % (sx(t), H - mb + 26, t))
    for t in (0, -1, -2, -3, -4, -5):
        if ymin < t < ymax:
            o.append('<text x="%d" y="%.1f" class="tick" text-anchor="end">%g</text>' % (ml - 6, sy(t) + 4, t))
    o.append('<text x="%d" y="%d" class="axis" text-anchor="middle">X level: median corrected LRR outside the PARs</text>' % ((ml + W - mr) / 2, H - 4))
    o.append('<text transform="translate(14 %d) rotate(-90)" class="axis" text-anchor="middle">Y level: median LRR over the Y-specific probes</text>' % ((mt + H - mb) / 2))
    plain, odd = [], []
    for r in rows:
        comp = r["complement"]
        mos = r["mosaic"] if r["mosaic"] != "NA" else ""
        ring = mos.startswith("mosaic")
        c = COL.get(comp, "#888")
        tip = "%s: %s%s; reported %s, %s" % (r["sample"], r["karyotype"], (" - " + mos) if mos else "", r["reported_sex"] if r["reported_sex"] != "NA" else "?", r["check"] if r["check"] != "NA" else "not checked")
        el = '<circle class="%s" cx="%.1f" cy="%.1f" r="%s" fill="%s" fill-opacity="0.85" data-t="%s"><title>%s</title></circle>' % (
            "odd" if (ring or comp not in ("XX", "XY")) else "pt", sx(num(r["x_lrr"])), sy(num(r["y_lrr"])), "4.2" if (ring or comp not in ("XX", "XY")) else "2.6", c, esc(tip), esc(tip))
        (odd if (ring or comp not in ("XX", "XY") or r["check"] not in ("agrees", "NA")) else plain).append(el)
    o.extend(plain)
    o.extend(odd)
    named = sorted([r for r in rows if r["check"] not in ("agrees", "NA")], key=lambda r: (num(r["y_lrr"]), num(r["x_lrr"])))
    for i, r in enumerate(named):                                               # the disagreements, named, alternating above and below
        x, y = sx(num(r["x_lrr"])), sy(num(r["y_lrr"]))
        dy = -11 if i % 2 == 0 else 20
        o.append('<text x="%.1f" y="%.1f" class="ptlab">%s %s</text>' % (x + 7, y + dy, r["sample"], r["karyotype"]))
    lx, ly = ml + 14, sy(-1.2)                                                  # the legend in the empty lower-left quarter
    for k, lab in (("XY", "46,XY"), ("XX", "46,XX"), ("XXY", "47,XXY"), ("X", "45,X")):
        o.append('<circle cx="%d" cy="%d" r="4" fill="%s"/><text x="%d" y="%d" class="tick">%s</text>' % (lx, ly, COL[k], lx + 9, ly + 4, lab))
        ly += 16
    o.append('<circle class="odd" cx="%d" cy="%d" r="4" fill="none"/><text x="%d" y="%d" class="tick">X or Y lost or gained in part of the cells</text>' % (lx, ly, lx + 9, ly + 4))
    o.append("</svg>")
    return "\n".join(o)


# ----------------------------------------------------------------------------------------------------------- the page

def numbers(d):
    n = {}
    qc, sx, pd, ev = d["qc"], d["sexs"], d["ped"], d["ev"]
    n["samples"] = qc["samples"]
    n["probes"] = qc["probes"]
    n["called"] = ev["called"]
    n["hard"] = qc["hard_fail"]
    n["flagged"] = qc["flagged"]
    n["events"] = ev["events"]
    n["samples_with_events"] = ev["samples_with_events"]
    n["roh_runs"] = ev["roh_runs"]
    auto = [e for e in d["events"] if e["chrom"] not in ("chrX", "chrY")]
    un = [e for e in auto if not tagged(e)]
    n["auto"] = len(auto)
    n["auto_tagged"] = len(auto) - len(un)
    n["untagged"] = len(un)
    n["untagged_samples"] = len({e["sample"] for e in un})
    n["untagged_baf"] = sum(1 for e in un if e["type"] in ("gain", "loss", "loh") and e["f_baf"] != "NA")
    n["untagged_roh"] = sum(1 for e in un if e["type"] == "roh")
    n["untagged_lrr_only"] = sum(1 for e in un if e["type"] in ("gain", "loss") and e["f_baf"] == "NA")
    n["whole"] = [e for e in auto if e["span"] == "whole"]
    checks = sx["checks"]
    n["sex_agree"] = checks.get("agrees", 0)
    n["sex_checked"] = sum(v for k, v in checks.items() if k != "NA")
    n["sex_disagree"] = [(r["sample"], r["karyotype"], r["reported_sex"], r["check"]) for r in d["sex"] if r["check"] not in ("agrees", "NA")]
    mos = [r for r in d["sex"] if r["mosaic"].startswith("mosaic")]
    n["loy"] = [num(r["f"]) for r in mos if "loss of Y" in r["mosaic"]]
    n["xloss"] = [num(r["f"]) for r in mos if "loss of an X" in r["mosaic"]]
    n["xgain"] = [num(r["f"]) for r in mos if "gain of an X" in r["mosaic"]]
    n["xxy"] = sum(1 for r in d["sex"] if r["complement"] == "XXY")
    n["x0"] = sum(1 for r in d["sex"] if r["complement"] == "X")
    n["comp"] = sx["complements"]
    n["x_modes"], n["y_modes"] = sx["x_modes"], sx["y_modes"]
    n["trios"], n["duos"] = pd["families"].get("trio", 0), pd["families"].get("duo", 0)
    n["classes"] = pd["classes"]
    n["pairs"] = pd["pairs_related"]
    n["ped_probes"] = pd["probes"]
    n["check_agree"] = pd["check"].get("agrees", 0)
    n["undisclosed"] = pd["check"].get("undisclosed full sibling pair", 0) + pd["check"].get("undisclosed second degree pair", 0)
    n["froh_high"] = sum(1 for r in d["kar"].values() if r["f_roh"] != "NA" and num(r["f_roh"]) > 0.0625)
    n["digest_n"] = len(d["digest"])
    n["one_copy"] = ev["calibration"]["one_copy_lrr"]
    n["assembly"] = d["manifest"].get("genome_label", d["manifest"]["assembly"])
    return n


def digest_table(d):
    rows = []
    for r in d["digest"]:
        f = num(r["f"])
        cells = "" if r["type"] == "roh" else pct(f) if f == f else ""
        rows.append("<tr><td>%s</td><td><code>%s</code></td><td>%s</td><td class=\"num\">%s</td><td>%s</td><td class=\"k\">%s</td></tr>" % (
            r["rank"], esc(r["sample"]), esc(r["event"]), r["size_mb"], cells, esc(r["karyotype"])))
    return "\n".join(rows)


def fig(src, caption, alt=""):
    return '<figure><a href="%s" class="zoom"><img src="%s" alt="%s" loading="lazy"></a><figcaption>%s</figcaption></figure>' % (src, src, esc(alt or caption), caption)


ASSETS = [("events/figures/HG01103/chr_chr2.png", "HG01103_chr2.png"), ("events/figures/HG01103/genome.png", "HG01103_genome.png"),
          ("events/figures/HG01103/chr_chr14.png", "HG01103_chr14.png"), ("events/figures/NA12748/chr_chr13.png", "NA12748_chr13.png"),
          ("events/figures/NA12748/chr_chr21.png", "NA12748_chr21.png"), ("events/figures/HG02322/genome.png", "HG02322_genome.png"),
          ("events/figures/NA19208/chr_chr9.png", "NA19208_chr9.png"), ("events/figures/NA20509/chr_chr21.png", "NA20509_chr21.png"),
          ("events/figures/HG02372/chr_chrX.png", "HG02372_chrX.png"), ("events/figures/NA06984/chr_chr18.png", "NA06984_chr18.png"),
          ("pedigree/figures/kinship_vs_ibs0.png", "kinship_vs_ibs0.png")]


def build(run, docs):
    d = load(run)
    n = numbers(d)
    os.makedirs(os.path.join(docs, "assets"), exist_ok=True)
    for src, dst in ASSETS:
        shutil.copyfile(os.path.join(run, src), os.path.join(docs, "assets", dst))
    karyo, n_drawn = karyogram_svg(d)
    whole = sorted(n["whole"], key=lambda e: -num(e["f"]))
    whole_words = "; ".join("%s of chromosome %s in %s (%s)" % (dict(gain="a gain", loss="a loss", loh="copy-neutral LOH", roh="a whole chromosome without heterozygosity")[e["type"]], e["chrom"][3:], e["sample"], "every cell" if e["type"] == "roh" else pct(num(e["f"]))) for e in whole)
    rng = lambda v: "%d-%d%%" % (round(100 * min(v)), round(100 * max(v))) if v else "-"
    t = Template(PAGE)
    page = t.substitute(
        version=__version__, github=GITHUB, assembly=esc(n["assembly"]),
        samples=fmt_int(n["samples"]), probes=fmt_int(n["probes"]), called=fmt_int(n["called"]), hard=n["hard"], flagged=fmt_int(n["flagged"]),
        sex_agree=fmt_int(n["sex_agree"]), sex_checked=fmt_int(n["sex_checked"]), trios=n["trios"], duos=n["duos"], pairs=n["pairs"], ped_probes=fmt_int(n["ped_probes"]),
        po=n["classes"].get("parent-offspring", 0), fs=n["classes"].get("full sibling", 0), second=n["classes"].get("second degree", 0), third=n["classes"].get("third degree", 0),
        check_agree=n["check_agree"], undisclosed=n["undisclosed"],
        events=fmt_int(n["events"]), samples_with_events=fmt_int(n["samples_with_events"]), auto=fmt_int(n["auto"]), auto_tagged=fmt_int(n["auto_tagged"]),
        untagged=n["untagged"], untagged_samples=n["untagged_samples"], untagged_baf=n["untagged_baf"], untagged_roh=n["untagged_roh"], untagged_lrr_only=n["untagged_lrr_only"],
        roh_runs=fmt_int(n["roh_runs"]), froh_high=n["froh_high"], n_drawn=fmt_int(n_drawn), whole_words=whole_words, n_whole=len(whole),
        xx=n["comp"].get("XX", 0) + n["comp"].get("XX mosaic", 0), xy=n["comp"].get("XY", 0) + n["comp"].get("XY mosaic", 0), xxy=n["xxy"], x0=n["x0"],
        n_loy=len(n["loy"]), loy_range=rng(n["loy"]), n_xloss=len(n["xloss"]), xloss_range=rng(n["xloss"]), n_xgain=len(n["xgain"]), xgain_range=rng(n["xgain"]),
        x_one="%.2f" % n["x_modes"]["one_copy"], x_two="%.2f" % n["x_modes"]["two_copies"], y_present="%.2f" % n["y_modes"]["present"], y_absent="%.2f" % n["y_modes"]["absent"],
        n_disagree=len(n["sex_disagree"]), disagree="; ".join("%s is %s in a reported %s" % (s, k, {"M": "man", "F": "woman"}.get(r, r)) for s, k, r, c in n["sex_disagree"]),
        digest_rows=digest_table(d), digest_n=n["digest_n"], karyogram=karyo, sexplot=sex_svg(d), one_copy="%.2f" % n["one_copy"],
    )
    with open(os.path.join(docs, "index.html"), "w") as fh:
        fh.write(page)
    print("docs/index.html: %d bytes; %d assets; %d events drawn" % (len(page.encode()), len(ASSETS), n_drawn))


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>arraykaryo: karyotypes, mosaics and pedigrees from a genotyping array</title>
<meta name="description" content="arraykaryo read 2,141 public 1000 Genomes samples from one Illumina array export: sex-chromosome complements, relatedness and families, mosaic gains, losses and copy-neutral LOH with cell fractions, karyotype strings and a report. This page is built from that run.">
<style>
:root{--paper:#fbfaf7;--ink:#1b1b1f;--muted:#5f5a52;--line:#e4dfd4;--card:#ffffff;--accent:#D55E00;--blue:#0072B2;--pink:#CC79A7;--gold:#E69F00;--green:#009E73;--hero1:#0f1620;--hero2:#15314a;--hero3:#0b4f6c;--code:#f1ede4;--max:1120px}
@media (prefers-color-scheme:dark){:root{--paper:#121418;--ink:#ecebe6;--muted:#a9a49a;--line:#2a2e36;--card:#1a1e25;--code:#20242c}}
*{box-sizing:border-box}html{scroll-behavior:smooth}
body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;-webkit-font-smoothing:antialiased}
a{color:var(--blue)}code{font:0.92em ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.wrap{max-width:var(--max);margin:0 auto;padding:0 20px}
nav.top{position:sticky;top:0;z-index:20;background:rgba(15,22,32,.92);backdrop-filter:blur(8px);color:#fff;border-bottom:1px solid rgba(255,255,255,.08)}
nav.top .wrap{display:flex;gap:18px;align-items:center;height:48px;overflow-x:auto;white-space:nowrap}
nav.top a{color:#dfe6ee;text-decoration:none;font-size:14px}nav.top a:hover{color:#fff}nav.top .brand{font-weight:700;letter-spacing:.02em;color:#fff;margin-right:8px}nav.top .gh{margin-left:auto}
header.hero{background:radial-gradient(1200px 600px at 10% -10%,#1e4d6b 0%,transparent 60%),linear-gradient(135deg,var(--hero1),var(--hero2) 55%,var(--hero3));color:#fff;padding:72px 0 56px}
.hero h1{font-size:clamp(40px,7vw,78px);line-height:1;margin:0 0 14px;letter-spacing:-.02em}
.hero h1 span{color:#ffb469}
.hero p.lead{font-size:clamp(18px,2.3vw,23px);max-width:46em;margin:0 0 10px;color:#e8eef4}
.hero p.sub{color:#b9c6d2;margin:0 0 28px;max-width:52em}
.cta{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:38px}
.cta a{display:inline-block;padding:10px 16px;border-radius:8px;text-decoration:none;font-weight:600;font-size:15px}
.cta a.primary{background:#ffb469;color:#1b1b1f}.cta a.ghost{border:1px solid rgba(255,255,255,.35);color:#fff}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:14px}
.tile{background:rgba(255,255,255,.07);border:1px solid rgba(255,255,255,.12);border-radius:12px;padding:16px 18px}
.tile b{display:block;font-size:30px;line-height:1.1;font-weight:700;letter-spacing:-.01em}
.tile span{color:#c9d4de;font-size:14px}
section{padding:54px 0 10px}section h2{font-size:clamp(26px,3.4vw,36px);letter-spacing:-.015em;margin:0 0 6px}
section p.kicker{color:var(--accent);font-weight:700;text-transform:uppercase;letter-spacing:.12em;font-size:12px;margin:0 0 6px}
section > .wrap > p{max-width:62em}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:22px;margin-top:18px}
figure{margin:0;background:var(--card);border:1px solid var(--line);border-radius:12px;padding:10px;box-shadow:0 1px 2px rgba(0,0,0,.04)}
figure img{width:100%;height:auto;display:block;border-radius:6px}
figcaption{font-size:14px;color:var(--muted);padding:10px 4px 2px;line-height:1.45}
figure.wide{grid-column:1/-1}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px 18px}
pre{background:var(--code);border:1px solid var(--line);border-radius:10px;padding:14px 16px;overflow:auto;font-size:14px;line-height:1.5}
.karyo{width:100%;height:auto;display:block}
.karyo .klabel{font-size:11px;fill:var(--muted);text-anchor:end;font-weight:600}
.karyo .ev{cursor:pointer}.karyo.hide-gain .ev-gain,.karyo.hide-loss .ev-loss,.karyo.hide-loh .ev-loh,.karyo.hide-roh .ev-roh{display:none}
.legend{display:flex;gap:16px;flex-wrap:wrap;font-size:14px;color:var(--muted);margin:10px 0 4px}
.legend label{display:inline-flex;align-items:center;gap:6px;cursor:pointer}.legend i{display:inline-block;width:14px;height:10px;border-radius:2px}
.sexplot{width:100%;max-width:760px;height:auto;display:block;margin:0 auto}
.sexplot .plotbg{fill:var(--card);stroke:var(--line)}.sexplot .mode{stroke:var(--muted);stroke-dasharray:2 3;stroke-width:.8}
.sexplot .modelab,.sexplot .tick,.sexplot .ptlab{font-size:11px;fill:var(--muted)}.sexplot .axis{font-size:12px;fill:var(--ink)}.sexplot .ptlab{fill:var(--ink);font-weight:600}
.sexplot circle{cursor:pointer}.sexplot circle.odd{stroke:var(--ink);stroke-width:1.2}
table.digest{width:100%;border-collapse:collapse;font-size:14px;background:var(--card);border:1px solid var(--line);border-radius:12px;overflow:hidden}
table.digest th{text-align:left;background:var(--code);padding:9px 10px;font-weight:600;position:sticky;top:48px}
table.digest td{padding:7px 10px;border-top:1px solid var(--line);vertical-align:top}table.digest td:nth-child(5){white-space:nowrap}table.digest td.num{text-align:right;white-space:nowrap}table.digest td.k{font-family:ui-monospace,Menlo,monospace;font-size:12.5px;color:var(--muted)}
.tablewrap{overflow:auto;max-height:620px;border-radius:12px}
#tip{position:fixed;pointer-events:none;background:#1b1b1f;color:#fff;font-size:13px;padding:6px 9px;border-radius:6px;max-width:360px;display:none;z-index:50;box-shadow:0 4px 14px rgba(0,0,0,.3)}
#lightbox{position:fixed;inset:0;background:rgba(10,12,16,.92);display:none;align-items:center;justify-content:center;z-index:60;padding:24px;cursor:zoom-out}
#lightbox img{max-width:100%;max-height:100%;border-radius:6px;background:#fff}
footer{border-top:1px solid var(--line);margin-top:48px;padding:28px 0 48px;color:var(--muted);font-size:14px}
ul.tight{margin:8px 0 0 0;padding-left:20px}ul.tight li{margin:4px 0}
.callout{border-left:4px solid var(--accent);background:var(--card);padding:12px 16px;border-radius:0 10px 10px 0;margin:16px 0}
</style>
</head>
<body>
<nav class="top"><div class="wrap">
<a class="brand" href="#top">arraykaryo</a>
<a href="#cohort">Cohort</a><a href="#validated">Validated</a><a href="#culture">Cell lines</a><a href="#sex">Sex chromosomes</a><a href="#pedigree">Pedigree</a><a href="#digest">Digest</a><a href="#limits">Limits</a><a href="#run">Run it</a>
<a class="gh" href="$github">GitHub &rarr;</a>
</div></nav>

<header class="hero" id="top"><div class="wrap">
<h1>Karyotypes, mosaics and pedigrees<br>from a <span>genotyping array</span></h1>
<p class="lead">arraykaryo reads a cohort's array export once and writes what a cytogeneticist would want to know about every sample: its sex-chromosome complement, who it is related to, and every near-chromosomal gain, loss, copy-neutral LOH or run of homozygosity, with the share of cells that carries it.</p>
<p class="sub">Everything on this page comes from one run of <code>arraykaryo run</code> on a public 1000 Genomes Project cohort: $samples samples, $probes probes, $assembly coordinates. Nothing was curated by hand; the numbers and drawings are read from the run's own tables.</p>
<div class="cta"><a class="primary" href="#cohort">See what it found</a><a class="ghost" href="$github">Code and docs</a><a class="ghost" href="$github/blob/main/docs/METHODS.md">Methods</a></div>
<div class="tiles">
<div class="tile"><b>$samples</b><span>samples read from one BCF; $called pass QC, $hard fail hard</span></div>
<div class="tile"><b>$sex_agree of $sex_checked</b><span>reported sexes agree with the genotypes; the $n_disagree that do not are explained</span></div>
<div class="tile"><b>$trios trios</b><span>and $duos duos assembled from the genotypes alone; the published pedigree holds</span></div>
<div class="tile"><b>$untagged events</b><span>of 2 Mb or more that survive the artefact map, on $untagged_samples samples, with cell fractions</span></div>
</div>
</div></header>

<section id="cohort"><div class="wrap">
<p class="kicker">The whole cohort at once</p>
<h2>Every event, on the chromosomes</h2>
<p>Each row is a chromosome of $assembly drawn to scale, bands shaded as a karyogram. Above it, every gain (orange), loss (blue) and copy-neutral LOH (pink) called in the cohort that survives the artefact map; below it, the runs of homozygosity of 10 Mb or more (gold). $n_drawn events on the picture, the longest nearest the chromosome; hover any one for the sample and its cell fraction. The X's and Y's lost or gained in part of the cells, the sex step's findings, are on the sex plot below rather than here. Whole chromosomes stand out as bars across the row: $whole_words.</p>
<div class="legend">
<label><input type="checkbox" checked data-k="gain"><i style="background:#D55E00"></i> gain</label>
<label><input type="checkbox" checked data-k="loss"><i style="background:#0072B2"></i> loss</label>
<label><input type="checkbox" checked data-k="loh"><i style="background:#CC79A7"></i> copy-neutral LOH</label>
<label><input type="checkbox" checked data-k="roh"><i style="background:#E69F00"></i> run of homozygosity (10 Mb or more)</label>
</div>
<div class="card" style="padding:8px 10px">$karyogram</div>
<p style="margin-top:14px">The run called $events events on $samples_with_events samples in all. Most are not on the picture: $auto_tagged of the $auto autosomal calls carry a tag - <em>recurrent</em> (the same reading in 2% or more of the cohort over the same megabase: a probe cluster the array reads oddly, or a common copy-number polymorphism), <em>probe-sparse</em>, or inside a listed artefact region (the MHC, the 8p23.1 and 17q21.31 inversions, centromere flanks, chromosome ends) - and almost all of those are read from the LRR alone, with no heterozygous BAF to split. They stay in the event table with their tags and off the digest. Of the $untagged that remain, $untagged_baf have the BAF split behind them, $untagged_lrr_only are copy changes in every cell read from the LRR alone, and $untagged_roh are runs of homozygosity of 10 Mb or more.</p>
</div></section>

<section id="validated"><div class="wrap">
<p class="kicker">Checked against sequencing</p>
<h2>The same mosaics, read from the array</h2>
<p>Two of these cell lines were read earlier from whole-genome sequencing by <a href="https://github.com/jlanej/triokaryo">triokaryo</a>. The array finds the same events at the same places. The cell fractions differ because the array's DNA and the sequenced DNA are different cultures of a lymphoblastoid line, and a clone's share changes with passage.</p>
<div class="grid2">
<figure class="wide"><a class="zoom" href="assets/HG01103_chr2.png"><img src="assets/HG01103_chr2.png" alt="HG01103 chromosome 2: LRR, BAF and heterozygous ratio with the loss of 2q22.1-q31.2 shaded" loading="lazy"></a><figcaption><b>HG01103, chromosome 2.</b> A loss of 2q22.1-q31.2 in 47% of cells and a further loss at q33.1 in 55%. The LRR (top) drops by a quarter; the heterozygous BAF (middle, green) splits into two bands exactly where a one-copy loss in 47% of the cells puts them (dashed); the heterozygous ratio (bottom) sags where calls near the bands are lost. Sequencing of another culture read the 2q loss at 67%.</figcaption></figure>
<figure><a class="zoom" href="assets/HG01103_chr14.png"><img src="assets/HG01103_chr14.png" alt="HG01103 chromosome 14 with the loss of 14q13.2-q21.1 shaded" loading="lazy"></a><figcaption><b>HG01103, chromosome 14.</b> A 6-Mb loss of 14q13.2-q21.1 in 55% of cells (sequencing: 64%). Her parents HG01101 and HG01102 carry none of these events: the losses arose in her line.</figcaption></figure>
<figure><a class="zoom" href="assets/HG01103_genome.png"><img src="assets/HG01103_genome.png" alt="HG01103 genome-wide tracks" loading="lazy"></a><figcaption><b>HG01103, genome-wide.</b> The same three losses in the context of every chromosome; the karyotype string above the tracks is written by the run: <code>mos 46,XX,del(2)(q22.1q31.2)[0.47]/46,XX,del(2)(q33.1)[0.55]/46,XX,del(14)(q13.2q21.1)[0.55]/46,XX</code>.</figcaption></figure>
<figure><a class="zoom" href="assets/NA12748_chr13.png"><img src="assets/NA12748_chr13.png" alt="NA12748 chromosome 13 with the gain of 13q31.2-q34 shaded" loading="lazy"></a><figcaption><b>NA12748, chromosome 13.</b> A gain of 13q31.2-q34 (26 Mb) in 19% of cells: a small rise of the LRR and a narrow BAF split, both significant over 25,000 probes. Sequencing of another culture read 69%.</figcaption></figure>
<figure><a class="zoom" href="assets/NA12748_chr21.png"><img src="assets/NA12748_chr21.png" alt="NA12748 chromosome 21 with the terminal loss of 21q22.3 shaded" loading="lazy"></a><figcaption><b>NA12748, chromosome 21.</b> A terminal loss of 21q22.3 (4 Mb) in 20% of cells, where sequencing read 59%. His partner NA12749 carries neither event.</figcaption></figure>
</div>
</div></section>

<section id="culture"><div class="wrap">
<p class="kicker">What a culture does to a cell line</p>
<h2>Whole chromosomes gained, and a chromosome without heterozygosity</h2>
<p>Lymphoblastoid lines acquire changes in culture, and the 1000 Genomes lines are no exception: $n_whole whole autosomes are changed in this cohort. arraykaryo lists them as findings with their cell fraction rather than treating them as errors; a digest for hand-off opens with them.</p>
<div class="grid2">
<figure><a class="zoom" href="assets/HG02322_genome.png"><img src="assets/HG02322_genome.png" alt="HG02322 genome-wide tracks with a gain of chromosome 12 in 84% of cells" loading="lazy"></a><figcaption><b>HG02322: a gain of chromosome 12 in 84% of cells.</b> The commonest acquired change of a lymphoblastoid line; the LRR of the whole chromosome rises and the heterozygous BAF splits into the bands of three copies in most cells.</figcaption></figure>
<figure><a class="zoom" href="assets/NA19208_chr9.png"><img src="assets/NA19208_chr9.png" alt="NA19208 chromosome 9 gained in 72% of cells" loading="lazy"></a><figcaption><b>NA19208: a gain of chromosome 9 in 72% of cells.</b> The BAF bands sit where f = 0.72 predicts for a third copy (dashed).</figcaption></figure>
<figure><a class="zoom" href="assets/NA20509_chr21.png"><img src="assets/NA20509_chr21.png" alt="NA20509 chromosome 21 at two copies without any heterozygous call" loading="lazy"></a><figcaption><b>NA20509: a whole chromosome 21 without a heterozygous call, at two copies.</b> An isodisomic uniparental disomy, or the two copies of a consanguineous union; the sample's genome-wide F<sub>ROH</sub> says which.</figcaption></figure>
<figure><a class="zoom" href="assets/HG02372_chrX.png"><img src="assets/HG02372_chrX.png" alt="HG02372 X chromosome at two copies with heterozygosity in a reported man" loading="lazy"></a><figcaption><b>HG02372: 47,XXY.</b> A reported man whose X sits at two copies with full heterozygosity while the Y is present. The sex step calls the complement; the report words the disagreement with the sample sheet rather than hiding it.</figcaption></figure>
<figure class="wide"><a class="zoom" href="assets/NA06984_chr18.png"><img src="assets/NA06984_chr18.png" alt="NA06984 chromosome 18 with a 10.9-Mb deletion of 18q22.2-q23 in every cell" loading="lazy"></a><figcaption><b>NA06984: a 10.9-Mb deletion of 18q22.2-q23 in every cell.</b> One copy: the LRR at the one-copy level of this platform ($one_copy) and no heterozygous call across the segment. A constitutional deletion of 18q at the distal end, the kind of event a karyotype would report.</figcaption></figure>
</div>
</div></section>

<section id="sex"><div class="wrap">
<p class="kicker">Sex chromosomes</p>
<h2>Every sample's X against its Y</h2>
<p>The X level of a sample (its median corrected LRR outside the pseudoautosomal regions) against its Y level (over the probes that are Y-specific in this cohort), from the run's <code>sex.tsv</code>. The cohort's own modes set the scale: one X at $x_one, two at $x_two; the Y present at $y_present and absent at $y_absent. Everything between the modes is a mosaic, confirmed on the X by the BAF split at heterozygous calls.</p>
<div class="card" style="padding:12px 8px">$sexplot</div>
<p style="margin-top:14px">Complements: 46,XX $xx, 46,XY $xy, 47,XXY $xxy, 45,X $x0. In culture, $n_loy men lost the Y from $loy_range of their cells, $n_xloss women lost an X from $xloss_range and $n_xgain gained one in $xgain_range. Against the reported sex, $sex_agree of $sex_checked agree. The ones that do not: $disagree.</p>
</div></section>

<section id="pedigree"><div class="wrap">
<p class="kicker">Relatedness</p>
<h2>Who is related to whom, from the genotypes alone</h2>
<div class="grid2">
<figure><a class="zoom" href="assets/kinship_vs_ibs0.png"><img src="assets/kinship_vs_ibs0.png" alt="KING-robust kinship against IBS0 for every pair of samples" loading="lazy"></a><figcaption>KING-robust kinship against IBS0 (the share of probes at which a pair is opposite homozygous) for every pair, from $ped_probes autosomal probes. Parent-offspring pairs have no IBS0; full siblings share the same kinship with some.</figcaption></figure>
<div class="card">
<p style="margin-top:0"><b>$pairs related pairs</b>: $po parent-offspring, $fs full-sibling, $second second-degree, $third third-degree. From them <b>$trios trios</b> and <b>$duos duos</b>, each trio oriented by its Mendelian-error rate against the cohort's own floor.</p>
<p>Against the published 1000 Genomes pedigree, $check_agree claimed relationships agree and no claimed parent is unrelated. The genotypes also find $undisclosed related pairs the table does not list, the cryptic relatedness of the collection, and one pair listed as third order that is in fact full siblings.</p>
<p style="margin-bottom:0">The trios the brief named are recovered: NA12878 with NA12891 and NA12892, NA19240 with NA19239 and NA19238, HG00733 with HG00731 and HG00732, HG01103 with HG01101 and HG01102. The "HG00702 trio" is a three-generation family, which is what the table says once read closely and what the genotypes find.</p>
</div>
</div>
</div></section>

<section id="digest"><div class="wrap">
<p class="kicker">For hand-off</p>
<h2>The digest: $digest_n findings, ranked as a cytogeneticist ranks</h2>
<p>A whole chromosome before an arm before a stretch; an event in every cell before one in part of the cells; a change of copy number before copy-neutral LOH before a run of homozygosity. Stretches under 10 Mb, mosaics under 10% of cells, runs under 20 Mb and the tagged readings are left out and counted with their reason. The X's and Y's lost or gained in culture are kept off it by default. Each row links to the sample's figures in the full report.</p>
<div class="tablewrap"><table class="digest">
<thead><tr><th>#</th><th>sample</th><th>event</th><th>Mb</th><th>cells</th><th>karyotype</th></tr></thead>
<tbody>$digest_rows</tbody>
</table></div>
</div></section>

<section id="limits"><div class="wrap">
<p class="kicker">What it can and cannot see</p>
<h2>Calibrated, and honest about the floor</h2>
<div class="grid2">
<figure><a class="zoom" href="calibration/calibration.png"><img src="calibration/calibration.png" alt="the smallest cell fraction recovered by event size and type on the simulated cohort" loading="lazy"></a><figcaption>On a simulated cohort at this array's probe spacing: a whole chromosome or arm is found from 10% of cells (copy-neutral LOH from 5%), a 5-Mb stretch from 15%, a 2-Mb stretch from 15-20%, with no false positive on 80 clean samples. <a href="calibration/">The calibration study.</a></figcaption></figure>
<div class="card">
<p style="margin-top:0"><b>Written down, not hidden.</b> This platform compresses one copy to an LRR of $one_copy rather than the textbook -0.55, so the cohort's own scale is used and fractions read from the LRR alone run low. The baseline BAF split of a sample sets what the BAF can confirm. Centromere flanks, the 9p11-q21 heterochromatin of CHM13 and a few common copy-number polymorphisms produce the commonest raw readings; the cohort's own artefact map tags them.</p>
<p>A low-level copy-neutral LOH pattern on 4q, at 8-14% of cells in 45 lines, sits at the floor of the LOH rule in BAF-noisy samples and is treated as the platform's until one is confirmed. The 0.06 BAF-SD threshold flags 12% of this cohort, which is read with the flag beside every call.</p>
<p style="margin-bottom:0">Every event is a reading of LRR, BAF and genotypes, to be confirmed by a clinical method before it means anything about a person. These are research cell lines.</p>
</div>
</div>
</div></section>

<section id="run"><div class="wrap">
<p class="kicker">Run it</p>
<h2>One command on your own export</h2>
<p>An Illumina array cohort exported to VCF/BCF with FORMAT GT, LRR and BAF (as <code>bcftools +gtc2vcf</code> writes), on GRCh37, GRCh38 or T2T-CHM13. A sample sheet with the reported sex and a pedigree are optional and are checked, never trusted over the genotypes.</p>
<pre>pip install git+$github
arraykaryo run --bcf cohort.bcf --out out --threads 8 --sample-sheet samples.tsv --pedigree trios.tsv
#   a cohort of thousands: add --pages digest;  an export without INFO/GC: add --fasta ref.fa

# or the container, which carries bcftools and everything else
docker run --rm -v "$$PWD:/work" ghcr.io/jlanej/arraykaryo:main arraykaryo run --bcf /work/cohort.bcf --out /work/out --threads 8</pre>
<p>Out come <code>qc/</code>, <code>sex/</code>, <code>pedigree/</code>, <code>events/</code> (every event with its evidence, every run of homozygosity, a karyotype string per sample, figures) and <code>report/</code>: a self-contained HTML report with per-sample pages, a guide to every column and rule, and the digest. The methods, with their thresholds and sources, are in <a href="$github/blob/main/docs/METHODS.md">METHODS.md</a>; every output column in <a href="$github/blob/main/docs/OUTPUT.md">OUTPUT.md</a>; the example's full numbers in <a href="$github/blob/main/docs/example/README.md">the example README</a>.</p>
</div></section>

<footer><div class="wrap">
<p>arraykaryo $version, MIT. The cohort is the 1000 Genomes Project's lymphoblastoid cell lines genotyped on an Illumina array and exported with <code>bcftools +gtc2vcf</code> on T2T-CHM13; the reported sexes and the known pedigree are the IGSR's public tables. All of it is public data, and nothing here is a clinical finding about a person.</p>
<p>Built by <code>scripts/build_site.py</code> from the run's tables: the karyogram and the sex plot are drawn from <code>events.all.tsv</code> and <code>sex.tsv</code> as SVG, the dense figures are the run's own PNGs.</p>
</div></footer>

<div id="tip"></div>
<div id="lightbox" aria-hidden="true"><img alt=""></div>
<script>
(function(){
  var tip=document.getElementById('tip');
  function show(e){var t=e.target.getAttribute&&e.target.getAttribute('data-t');if(!t){tip.style.display='none';return;}tip.textContent=t;tip.style.display='block';move(e);}
  function move(e){var x=e.clientX+14,y=e.clientY+14;if(x+370>window.innerWidth)x=e.clientX-370;tip.style.left=x+'px';tip.style.top=y+'px';}
  document.querySelectorAll('svg.karyo, svg.sexplot').forEach(function(s){s.addEventListener('mouseover',show);s.addEventListener('mousemove',move);s.addEventListener('mouseout',function(){tip.style.display='none';});});
  var k=document.querySelector('svg.karyo');
  document.querySelectorAll('.legend input').forEach(function(c){c.addEventListener('change',function(){k.classList.toggle('hide-'+c.getAttribute('data-k'),!c.checked);});});
  var lb=document.getElementById('lightbox'),im=lb.querySelector('img');
  document.querySelectorAll('a.zoom').forEach(function(a){a.addEventListener('click',function(e){e.preventDefault();im.src=a.getAttribute('href');im.alt=a.querySelector('img').alt;lb.style.display='flex';});});
  lb.addEventListener('click',function(){lb.style.display='none';im.src='';});
  document.addEventListener('keydown',function(e){if(e.key==='Escape'){lb.style.display='none';im.src='';}});
})();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    build(sys.argv[1], sys.argv[2])
