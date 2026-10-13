"""arraykaryo report: one self-contained HTML page per cohort (figures embedded as base64, readable from a file:// URL),
per-sample pages for flagged samples and samples with events (--pages all for every sample), the guide that explains every
column and colour with pattern cards of what each event kind looks like on an array, and the digest of the largest events
for hand-off, ranked as triokaryo's digest ranks: a whole chromosome first, then an arm, then a stretch by size; an event in
every cell before one in part of the cells; a change of copy number before a copy-neutral LOH before a run of homozygosity."""
import base64
import html
import io
import os
import time

import numpy as np

from . import __version__
from .cache import Cache
from .events import read_events, read_karyotypes, read_rohs
from .pedigree import read_families
from .qc import THRESHOLDS, read_sample_qc
from .sex import read_sex
from .util import NA, finite, log, mark_done, num, read_json, read_tsv, write_tsv

CSS = ("html{background:#fff;color-scheme:light}body{font-family:Helvetica,Arial,sans-serif;max-width:1280px;margin:1.5em auto;padding:0 1em;color:#222;background:#fff;line-height:1.35}"
       "table{border-collapse:collapse;font-size:12px;margin:0.5em 0}th,td{border:1px solid #ddd;padding:2px 6px;text-align:left;vertical-align:top}th{background:#f3f3f3;cursor:pointer}"
       "img{max-width:100%}figure{margin:1em 0}figcaption{font-size:12px;color:#444}h2{margin-top:1.6em;border-bottom:1px solid #ddd;padding-bottom:0.2em}h3{margin-top:1.2em}"
       ".mock{background:#fff3cd;padding:0.5em;border:1px solid #e0c060;margin:0.5em 0}.tiles{display:flex;flex-wrap:wrap;gap:10px;margin:1em 0}"
       ".tile{border:1px solid #ddd;border-radius:6px;padding:0.5em 0.9em;min-width:120px}.tile b{font-size:1.5em;display:block}.tile span{font-size:11px;color:#555}"
       ".muted{color:#777;font-size:12px}.flag{color:#b03a2e}.ok{color:#2a7a2a}.box{background:#f6f6f6;border:1px solid #e3e3e3;padding:0.6em 1em;margin:1em 0}"
       ".ev{border:1px solid #ddd;border-radius:6px;padding:0.6em 1em;margin:1.2em 0}.facts td:first-child{color:#666;white-space:nowrap}.tag{display:inline-block;background:#eee;border-radius:9px;padding:0 8px;font-size:11px;margin-left:6px}"
       "input#filter{width:40%;padding:3px}.nav a{margin-right:1em}")
JS = ("""<script>
function sortTable(th){var t=th.closest('table'),i=Array.from(th.parentNode.children).indexOf(th),rows=Array.from(t.tBodies[0].rows),asc=!th.asc;th.asc=asc;
rows.sort(function(a,b){var x=a.cells[i].textContent.trim(),y=b.cells[i].textContent.trim();var nx=parseFloat(x),ny=parseFloat(y);
if(!isNaN(nx)&&!isNaN(ny))return asc?nx-ny:ny-nx;if(x==='NA')return 1;if(y==='NA')return -1;return asc?x.localeCompare(y):y.localeCompare(x);});
rows.forEach(function(r){t.tBodies[0].appendChild(r);});}
function filterTable(inp,id){var q=inp.value.toLowerCase(),t=document.getElementById(id);Array.from(t.tBodies[0].rows).forEach(function(r){r.style.display=r.textContent.toLowerCase().indexOf(q)>=0?'':'none';});}
document.addEventListener('DOMContentLoaded',function(){document.querySelectorAll('table.sortable th').forEach(function(th){th.addEventListener('click',function(){sortTable(th);});});});
</script>""")
TYPE_WORD = {"gain": "gain", "loss": "loss", "loh": "copy-neutral loss of heterozygosity", "roh": "run of homozygosity"}
TIER = {"whole": 0, "p": 1, "q": 1, "stretch": 2}
TYPE_ORDER = {"gain": 0, "loss": 0, "loh": 1, "roh": 2}


def esc(x):
    return html.escape(str(x))


def embed_png(path, max_width=1400):
    """A PNG on disk as a data URI, downscaled to at most max_width pixels wide."""
    if not path or not os.path.exists(path):
        return ""
    try:
        from PIL import Image
        im = Image.open(path)
        if im.width > max_width:
            im = im.convert("RGB")
            im.thumbnail((max_width, int(im.height * max_width / im.width)))
            buf = io.BytesIO()
            im.save(buf, format="PNG", optimize=True)
            data = buf.getvalue()
        else:
            data = open(path, "rb").read()
    except Exception:
        data = open(path, "rb").read()
    return "data:image/png;base64," + base64.b64encode(data).decode()


def figure(path, caption, alt="", max_width=1400):
    uri = embed_png(path, max_width)
    if not uri:
        return '<p class="muted">figure not available: %s</p>' % esc(os.path.basename(path or ""))
    return '<figure><img src="%s" alt="%s"><figcaption>%s</figcaption></figure>' % (uri, esc(alt or caption[:60]), caption)


def table(cols, rows, id_=None, sortable=True, link=None, headers=None, max_rows=None):
    """rows: dicts; link: (column, fn(row) -> href or ''); headers: display names."""
    h = headers or cols
    w = ['<table%s%s><thead><tr>' % ((' id="%s"' % id_) if id_ else "", ' class="sortable"' if sortable else "") + "".join("<th>%s</th>" % esc(x) for x in h) + "</tr></thead><tbody>"]
    shown = rows if max_rows is None else rows[:max_rows]
    for r in shown:
        cells = []
        for c in cols:
            v = r.get(c, NA)
            if link and c == link[0]:
                href = link[1](r)
                cells.append('<td><a href="%s">%s</a></td>' % (href, esc(v)) if href else "<td>%s</td>" % esc(v))
            else:
                cells.append("<td>%s</td>" % esc(v))
        w.append("<tr>" + "".join(cells) + "</tr>")
    w.append("</tbody></table>")
    if max_rows is not None and len(rows) > max_rows:
        w.append('<p class="muted">%d of %d rows shown; the full table is in the TSV.</p>' % (max_rows, len(rows)))
    return "\n".join(w)


def what(e):
    """The event in words."""
    t = e["type"]
    c = e["chrom"][3:]
    f = num(e.get("f", ""))
    cells = "" if (not finite(f) or f >= 0.9) else " in %.0f%% of cells" % (100 * f)
    size = "%.1f-Mb " % (int(e["size_bp"]) / 1e6) if e["span"] == "stretch" else ""
    bands = e.get("bands", "")
    where = {"whole": "chromosome " + c, "p": "the short arm of chromosome %s (%sp)" % (c, c), "q": "the long arm of chromosome %s (%sq)" % (c, c)}.get(e["span"], "%s%s" % (c, bands if bands and bands != NA else ""))
    if e.get("source") == "sex":
        return e.get("note", "").replace("constitutional ", "") if e.get("note", NA) != NA else "%s of %s" % (t, where)
    if t == "roh" and e["span"] == "whole":
        return "chromosome %s without heterozygosity (an isodisomic uniparental disomy, or consanguinity)" % c
    if t == "roh":
        return "a %srun of homozygosity on %s" % (size, where)
    if t == "loss" and finite(f) and f >= 0.9 and e["span"] == "stretch":
        return "a %sdeletion of %s in every cell (one copy)" % (size, where)
    head = {"gain": "gain of %s", "loss": "loss of %s", "loh": "copy-neutral loss of heterozygosity of %s"}[t] % where
    if e["span"] in ("whole",) and t == "gain" and finite(f) and f >= 0.9:
        head = "an extra copy of chromosome %s (trisomy %s)" % (c, c)
    return ("mosaic " if cells else "") + size + head + cells


def rank_key(e):
    f = num(e.get("f", ""))
    return (TIER.get(e["span"], 2), 0 if (finite(f) and f >= 0.9) else 1, TYPE_ORDER.get(e["type"], 3), -int(e["size_bp"]), -(f if finite(f) else 0.0))


def mock_note(C):
    d = os.path.dirname(os.path.abspath(C.manifest.get("export", "")))
    p = os.path.join(d, "MOCK_DATA.txt")
    return open(p).read().strip() if os.path.exists(p) else ""


class Report:
    def __init__(self, out, cache_dir=None):
        self.out = out
        self.cache_dir = cache_dir or os.path.join(out, "cache")
        self.C = Cache(self.cache_dir)
        self.G = self.C.genome
        self.qc = read_sample_qc(out) if os.path.exists(os.path.join(out, "qc", "samples.tsv")) else {}
        self.qc_summary = read_json(os.path.join(out, "qc", "summary.json")) if os.path.exists(os.path.join(out, "qc", "summary.json")) else {}
        self.sex = read_sex(out)
        self.sex_summary = read_json(os.path.join(out, "sex", "summary.json")) if os.path.exists(os.path.join(out, "sex", "summary.json")) else {}
        self.families = read_families(out)
        p = os.path.join(out, "pedigree", "check.tsv")
        self.check = read_tsv(p) if os.path.exists(p) else None
        p = os.path.join(out, "pedigree", "relatedness.tsv")
        self.relatedness = read_tsv(p) if os.path.exists(p) else []
        self.ped_summary = read_json(os.path.join(out, "pedigree", "summary.json")) if os.path.exists(os.path.join(out, "pedigree", "summary.json")) else {}
        self.events = read_events(out)
        self.rohs = read_rohs(out)
        self.kar = read_karyotypes(out)
        self.ev_summary = read_json(os.path.join(out, "events", "summary.json")) if os.path.exists(os.path.join(out, "events", "summary.json")) else {}
        self.mock = mock_note(self.C)
        self.pages = set()

    def sample_href(self, s, from_sub=False):
        return ("%s.html" if from_sub else "samples/%s.html") % s if s in self.pages else ""

    # ----------------------------------------------------------------------------------------------- the cohort page

    def tiles(self):
        S = self.C.n_samples
        n_pass = sum(1 for r in self.qc.values() if r.get("pass") == "1")
        n_hard = sum(1 for r in self.qc.values() if r.get("hard_fail") == "1")
        flags = {}
        for r in self.qc.values():
            for f in r.get("flags", NA).split(","):
                if f and f != NA:
                    flags[f] = flags.get(f, 0) + 1
        fam = {}
        for f in self.families:
            fam[f["kind"]] = fam.get(f["kind"], 0) + 1
        by = {}
        for e in self.events:
            k = TYPE_WORD[e["type"]] + (", whole chromosome" if e["span"] == "whole" else ", arm" if e["span"] in ("p", "q") else ", stretch")
            by[k] = by.get(k, 0) + 1
        t = [("samples", S, ""), ("passing QC", n_pass, "%d flagged, %d not to be read" % (S - n_pass, n_hard)),
             ("families", sum(v for k, v in fam.items() if k in ("trio", "duo")), ", ".join("%s %d" % kv for kv in sorted(fam.items())) or "none"),
             ("samples with events", len({e["sample"] for e in self.events}), "%d events" % len(self.events)),
             ("not plainly XX or XY", sum(1 for r in self.sex.values() if r["mosaic"] != NA or r["complement"] not in ("XX", "XY")), "%d sex checks disagree" % sum(1 for r in self.sex.values() if r["check"] not in ("agrees", NA, "not read")))]
        w = ['<div class="tiles">'] + ['<div class="tile"><b>%s</b>%s<span>%s</span></div>' % (esc(v), esc(k), esc(sub)) for k, v, sub in t] + ["</div>"]
        if flags:
            w.append('<p class="muted">Flags: %s.</p>' % ", ".join("%s %d" % kv for kv in sorted(flags.items())))
        if by:
            w.append('<p class="muted">Events: %s.</p>' % ", ".join("%s %d" % kv for kv in sorted(by.items())))
        return "\n".join(w)

    def section_qc(self):
        w = ["<h2 id=\"qc\">Array quality</h2>"]
        s = self.qc_summary
        if not s:
            return "\n".join(w + ['<p class="muted">The qc step has not run.</p>'])
        w.append("<p>%d probes, %d kept for LRR and %d for genotypes (%s). %d of %d samples pass; a flag is shown beside every call of a flagged sample; a hard failure (LRR SD above %g) is not read. "
                 "Corrections: the GC wave per sample (the correlation of LRR with the 1-Mb GC falls from a median of %.3f to %.3f), then the cohort reference per probe%s.</p>" % (
                     s["probes"], s["probes_keep"], s["probes_keep_gt"], ", ".join("%s %d" % kv for kv in sorted(s["reasons"].items()) if kv[0] != NA), s["passing"], s["samples"],
                     s["thresholds"]["sample_lrr_sd_fail"]["value"], abs(s["cohort"]["gc_corr_before"]["median"]), abs(s["cohort"]["gc_corr_after"]["median"]),
                     " (per batch: %s)" % ", ".join(s["batches"]) if s.get("use_batch") else ""))
        w.append(figure(os.path.join(self.out, "qc", "figures", "qc_distributions.png"), "Sample QC metrics with the thresholds that flag (dotted) or fail (dashed) a sample; the black line is the cohort median. Every metric is defined in the guide."))
        w.append(figure(os.path.join(self.out, "qc", "figures", "gc_wave.png"), "The GC wave before and after correction for a typical sample and for the one with the strongest wave."))
        w.append("<h3>Thresholds</h3>")
        w.append(table(["name", "value", "source"], [dict(name=k, value=v["value"], source=v["source"]) for k, v in s["thresholds"].items()], sortable=False))
        w.append("<h3>Samples</h3><p>Click a heading to sort; type to filter. <input id=\"filter\" placeholder=\"filter: a sample, a flag, a batch\" onkeyup=\"filterTable(this,'samples')\"></p>")
        rows = [self.qc[x] for x in self.C.samples if x in self.qc]
        cols = [c for c in (list(rows[0].keys()) if rows else []) if c not in ("y_lrr_raw", "y_call_rate", "lrr_sd_diff", "gc_corr_before", "gc_corr_after")]
        w.append(table(cols, rows, id_="samples", link=("sample", lambda r: self.sample_href(r["sample"]))))
        return "\n".join(w)

    def section_sex(self):
        w = ["<h2 id=\"sex\">Sex chromosomes</h2>"]
        if not self.sex:
            return "\n".join(w + ['<p class="muted">The sex step has not run.</p>'])
        s = self.sex_summary
        xm, ym = s.get("x_modes", {}), s.get("y_modes", {})
        w.append("<p>The complement of every sample from the X outside the pseudoautosomal regions and the Y-specific probes, on the cohort's own copy scale: one X at LRR %s, two at %s (%s per copy); a Y present at %s, absent at %s. %s. "
                 "Complements: %s. Checks against the reported sex: %s.</p>" % (_f(xm.get("one_copy")), _f(xm.get("two_copies")), _f(xm.get("scale_per_copy")), _f(ym.get("present")), _f(ym.get("absent")), esc(s.get("y_probes", "")),
                                                                                 ", ".join("%s %d" % kv for kv in sorted(s.get("complements", {}).items())), ", ".join("%s %d" % kv for kv in sorted(s.get("checks", {}).items()))))
        w.append(figure(os.path.join(self.out, "sex", "figures", "x_vs_y.png"), "The X level against the Y level of every sample, coloured by the X heterozygosity relative to the autosomes; samples not plainly XX or XY are ringed and labelled."))
        odd = [r for r in self.sex.values() if r["mosaic"] != NA or r["complement"] not in ("XX", "XY") or r["check"] not in ("agrees", NA, "not read")]
        w.append("<h3>Every sample not plainly XX or XY, or disagreeing with its reported sex</h3>")
        if odd:
            w.append(table(["sample", "reported_sex", "complement", "karyotype", "mosaic", "f", "f_lrr", "f_baf", "x_copies", "y_copies", "x_het_rel", "check", "note"], odd, link=("sample", lambda r: self.sample_href(r["sample"]))))
        else:
            w.append('<p class="ok">none</p>')
        return "\n".join(w)

    def section_pedigree(self):
        w = ["<h2 id=\"pedigree\">Relatedness and families</h2>"]
        if not self.ped_summary:
            return "\n".join(w + ['<p class="muted">The pedigree step has not run.</p>'])
        s = self.ped_summary
        w.append("<p>KING-robust kinship and IBS0 of every pair from %d autosomal probes; %d related pairs (%s); families: %s. The parent-offspring / full-sibling cut on IBS0 is %.4f (the error floor %.4f); %d unrelated passing samples; %d probes fail Hardy-Weinberg among them.</p>" % (
            s["probes"], s["pairs_related"], ", ".join("%s %d" % kv for kv in sorted(s["classes"].items())) or "none", ", ".join("%s %d" % kv for kv in sorted(s["families"].items())) or "none", s["po_fs_cut"], s["ibs0_floor"], s["unrelated_passing"], s["hwe_failed_probes"]))
        w.append(figure(os.path.join(self.out, "pedigree", "figures", "kinship_vs_ibs0.png"), "Kinship against IBS0 for every pair; the class boundaries dotted, the parent-offspring / full-sibling cut dashed."))
        w.append("<h3>Families found</h3>")
        w.append(table(["family", "kind", "members", "kid", "dad", "mom", "kid_sex", "mendelian_error_rate", "note"], self.families) if self.families else '<p class="muted">none</p>')
        if self.check is not None:
            bad = [r for r in self.check if r["status"] != "agrees"]
            w.append("<h3>The check against the given pedigree</h3><p>%d claimed relationships agree; %d rows disagree or add something.</p>" % (sum(1 for r in self.check if r["status"] == "agrees"), len(bad)))
            w.append(table(["kid", "relative", "claimed", "found", "kinship", "status", "note"], bad) if bad else '<p class="ok">every claimed relationship agrees and nothing undisclosed was found</p>')
        return "\n".join(w)

    def section_events(self):
        w = ["<h2 id=\"events\">Near-chromosomal events</h2>"]
        if not self.ev_summary:
            return "\n".join(w + ['<p class="muted">The events step has not run.</p>'])
        s = self.ev_summary
        cal = s["calibration"]
        w.append("<p>%d events on %d of %d samples read (%s). Cell fractions from the BAF split first, from the LRR second on the cohort's copy scale: one copy at LRR %.3f (%s), three copies at %.3f (%s). "
                 "Floors: %g Mb and %d probes. The event-count flag cut is %s events.</p>" % (s["events"], s["samples_with_events"], s["called"], ", ".join("%s %d" % kv for kv in sorted(s["by_type"].items())) or "none",
                                                                                            cal["one_copy_lrr"], esc(cal["one_copy_source"]), cal["three_copies_lrr"], esc(cal["three_copies_source"]), s["params"]["min_mb"], s["params"]["min_probes"], _f(s.get("event_count_cut"))))
        w.append(figure(os.path.join(self.out, "report", "figures", "landscape.png"), "Events per chromosome by type, and one row per sample with an event."))
        w.append("<p>Click a heading to sort; type to filter. <input id=\"filter2\" placeholder=\"filter: a sample, a chromosome, a type, a word of a note\" onkeyup=\"filterTable(this,'events')\"></p>")
        rows = [dict(e, event=what(e), size_mb="%.1f" % (int(e["size_bp"]) / 1e6)) for e in self.events]
        cols = ["sample", "event", "chrom", "start", "end", "size_mb", "span", "bands", "type", "f", "f_lrr", "f_baf", "constitutional", "n_probes", "n_het", "lrr", "d", "d_llr", "het_ratio", "source", "artefact", "note", "sample_flags"]
        w.append(table(cols, rows, id_="events", link=("sample", lambda r: self.sample_href(r["sample"]))))
        return "\n".join(w)

    def section_roh(self):
        w = ["<h2 id=\"roh\">Runs of homozygosity</h2>"]
        if not self.kar:
            return "\n".join(w + ['<p class="muted">The events step has not run.</p>'])
        fr = np.array([num(r["f_roh"]) for r in self.kar.values()])
        fr = fr[np.isfinite(fr)]
        w.append("<p>%d runs of 1.5 Mb or more on %d samples (PLINK's sliding window); the genome fraction in runs, F_ROH, has a median of %.4f (robust SD %.4f). Runs of 10 Mb or more are events; a whole chromosome without heterozygosity is an isodisomic uniparental disomy or consanguinity.</p>" % (
            len(self.rohs), len({r["sample"] for r in self.rohs}), float(np.median(fr)) if fr.size else np.nan, float(1.4826 * np.median(np.abs(fr - np.median(fr)))) if fr.size else np.nan))
        w.append(figure(os.path.join(self.out, "report", "figures", "froh.png"), "The genome fraction in runs of homozygosity (F_ROH) per sample."))
        long = [r for r in self.rohs if int(r["size_bp"]) >= 5_000_000]
        w.append("<h3>Runs of 5 Mb or more</h3>")
        w.append(table(["sample", "chrom", "start", "end", "size_bp", "n_probes", "n_het", "bands", "lrr", "whole_chromosome", "is_event"], sorted(long, key=lambda r: -int(r["size_bp"])), link=("sample", lambda r: self.sample_href(r["sample"])), max_rows=300) if long else '<p class="muted">none</p>')
        return "\n".join(w)

    def section_aside(self):
        w = ["<h2 id=\"aside\">Samples set aside</h2>"]
        hard = [r for r in self.qc.values() if r.get("hard_fail") == "1"]
        outl = [r for r in self.qc.values() if r.get("hard_fail") != "1" and any(f in r.get("flags", "") for f in ("event_count", "lrr_sd", "duplicate"))]
        w.append("<p>Not read (a hard QC failure): %d. Read with a warning (noise above the mosaic flag, an event count far above the cohort's, a duplicate): %d. Their calls stay in the tables with the flag beside them, as triokaryo's cohort guards do.</p>" % (len(hard), len(outl)))
        rows = hard + outl
        w.append(table(["sample", "call_rate", "lrr_sd", "baf_sd", "wf", "n_events", "hard_fail", "flags"], rows, link=("sample", lambda r: self.sample_href(r["sample"]))) if rows else '<p class="ok">none</p>')
        return "\n".join(w)

    def write_index(self):
        w = ['<!doctype html><html><head><meta charset="utf-8"><title>arraykaryo: %d samples</title><style>%s</style>%s</head><body>' % (self.C.n_samples, CSS, JS)]
        if self.mock:
            w.append('<div class="mock">%s</div>' % esc(self.mock))
        w.append("<h1>arraykaryo: %d samples, %d events</h1>" % (self.C.n_samples, len(self.events)))
        w.append('<p class="muted">arraykaryo %s, %s. Export: %s (%s, %d probes). Every call is a reading of the array\'s LRR, BAF and genotypes, to be confirmed by a clinical karyotype or another assay before it is reported.</p>' % (
            __version__, time.strftime("%Y-%m-%d"), esc(os.path.basename(self.C.manifest.get("export", ""))), esc(self.C.manifest.get("genome_label", "")), sum(self.C.n_probes.values())))
        w.append('<p class="nav"><a href="#qc">quality</a><a href="#sex">sex chromosomes</a><a href="#pedigree">families</a><a href="#events">events</a><a href="#roh">runs of homozygosity</a><a href="#aside">set aside</a><a href="guide.html">the guide</a><a href="digest.html">the digest</a></p>')
        w.append(self.tiles())
        w.append(self.section_qc())
        w.append(self.section_sex())
        w.append(self.section_pedigree())
        w.append(self.section_events())
        w.append(self.section_roh())
        w.append(self.section_aside())
        w.append('<h2>Tables</h2><p class="muted">qc/probes.tsv, qc/samples.tsv, sex/sex.tsv, pedigree/relatedness.tsv, pedigree/families.tsv, pedigree/pedigree.inferred.tsv, pedigree/check.tsv, events/events.all.tsv, events/roh.all.tsv, events/karyotypes.tsv; every column is explained in <a href="guide.html">the guide</a>.</p></body></html>')
        with open(os.path.join(self.out, "report", "index.html"), "w") as fh:
            fh.write("\n".join(w))

    # ----------------------------------------------------------------------------------------------- per-sample pages

    def write_sample_page(self, s):
        q = self.qc.get(s, {})
        sx = self.sex.get(s, {})
        k = self.kar.get(s, {})
        evs = [e for e in self.events if e["sample"] == s]
        rohs = [r for r in self.rohs if r["sample"] == s]
        fams = [f for f in self.families if s in f["members"].split(",")]
        w = ['<!doctype html><html><head><meta charset="utf-8"><title>arraykaryo %s</title><style>%s</style>%s</head><body>' % (esc(s), CSS, JS)]
        if self.mock:
            w.append('<div class="mock">%s</div>' % esc(self.mock))
        w.append('<p class="nav"><a href="../index.html">the cohort</a><a href="../guide.html">the guide</a></p>')
        w.append("<h1>%s: %s</h1>" % (esc(s), esc(k.get("karyotype", NA) if k else NA)))
        flags = q.get("flags", NA)
        w.append("<p>%s. QC: %s%s. Sex chromosomes: %s%s.</p>" % (
            "Not to be read (a hard QC failure); the calls below are shown for the record" if q.get("hard_fail") == "1" else ("Flagged: " + flags) if flags != NA else "Passes QC",
            "call rate %s, LRR SD %s (cohort median %s), BAF SD %s, wave factor %s" % (q.get("call_rate"), q.get("lrr_sd"), _f(self.qc_summary.get("cohort", {}).get("lrr_sd", {}).get("median")), q.get("baf_sd"), q.get("wf")),
            ", batch %s" % q["batch"] if q.get("batch", NA) != NA else "",
            esc(sx.get("karyotype", NA)) if sx else "not read", (" (" + esc(sx["check"]) + ")") if sx and sx.get("check") not in (NA, "") else ""))
        if fams:
            w.append("<p>Families: %s.</p>" % "; ".join("%s (%s: %s)" % (f["family"], f["kind"], f["members"]) for f in fams))
        w.append("<h2>Events</h2>")
        if evs:
            w.append("<ul>" + "".join("<li><b>%s</b> - %s%s%s</li>" % (esc(what(e)), "%s:%s-%s" % (e["chrom"], e["start"], e["end"]), (", " + e["bands"]) if e["bands"] not in (NA, "") else "",
                                                                      "; LRR %s, d %s, %s probes, %s heterozygous calls; cell fraction %s (LRR %s, BAF %s)%s%s" % (e["lrr"], e["d"], e["n_probes"], e["n_het"], e["f"], e["f_lrr"], e["f_baf"],
                                                                                                                                                             ("; " + esc(e["note"])) if e["note"] != NA else "", ("; artefact-prone region: " + e["artefact"]) if e["artefact"] != NA else "")) for e in evs) + "</ul>")
        else:
            w.append('<p class="ok">no event</p>')
        if rohs:
            w.append("<p>Runs of homozygosity: %d (%s); F_ROH %s.</p>" % (len(rohs), ", ".join("%s %.1f Mb" % (r["chrom"], int(r["size_bp"]) / 1e6) for r in rohs[:12]), k.get("f_roh", NA)))
        fdir = os.path.join(self.out, "events", "figures", s)
        if os.path.exists(os.path.join(fdir, "genome.png")):
            w.append(figure(os.path.join(fdir, "genome.png"), "Every chromosome: corrected LRR, BAF and the heterozygous ratio; events as labelled bars."))
            for e in sorted(evs, key=lambda e: (self.G.chrom_order().get(e["chrom"], 99), int(e["start"]))):
                p = os.path.join(fdir, "chr_%s.png" % e["chrom"])
                if os.path.exists(p):
                    w.append(figure(p, "%s: %s." % (e["chrom"], esc(what(e)))))
        p = os.path.join(self.out, "sex", "figures", "tracks_%s.png" % s)
        if os.path.exists(p):
            w.append(figure(p, "The X and Y tracks: corrected LRR and BAF along the X, LRR along the Y."))
        w.append("<h2>QC row</h2>" + table(list(q.keys()), [q], sortable=False) if q else "")
        w.append("</body></html>")
        with open(os.path.join(self.out, "report", "samples", "%s.html" % s), "w") as fh:
            fh.write("\n".join(w))

    # ------------------------------------------------------------------------------------------------------- the guide

    def write_guide(self):
        from .plots import fig_patterns
        fig = fig_patterns(self.out)
        w = ['<!doctype html><html><head><meta charset="utf-8"><title>arraykaryo guide</title><style>%s</style></head><body>' % CSS]
        w.append('<p class="nav"><a href="index.html">the cohort</a><a href="digest.html">the digest</a></p><h1>The guide: how to read an arraykaryo report</h1>')
        w.append('<p>arraykaryo reads a genotyping array\'s export - per probe and sample the genotype (GT), the log R ratio (LRR: the probe\'s total intensity against the cohort, on a log2 scale; 0 is two copies, about -0.5 one copy and about +0.35 three on the array\'s compressed scale) and the B-allele frequency (BAF: the share of the B allele in the intensity; 0, 0.5 and 1 for the three genotypes) - and reports, per sample, the array\'s quality, the sex-chromosome complement, the relatives in the cohort, and the near-chromosomal events: whole chromosomes, arms and segments of 2 Mb or more gained, lost or made homozygous, in every cell (constitutional) or in part of the cells (mosaic, with the cell fraction).</p>')
        w.append("<h2>What each event looks like</h2>")
        w.append(figure(fig, "Pattern cards (simulated, typical array noise): a gain lifts the LRR and splits the heterozygous BAF band towards 1/3 and 2/3 (fully for a trisomy, less for a mosaic); a loss lowers the LRR and splits the band outward; a copy-neutral LOH leaves the LRR flat and splits the band; a run of homozygosity has no heterozygous calls and a flat LRR; 47,XXY shows an X at two copies with heterozygous calls and a Y present; a mosaic loss of Y sits between the male and the absent Y level."))
        w.append("<p>The cell fraction f follows from the band split d (the BAF's distance from 0.5 at heterozygous calls): gain f = 4d/(1-2d), loss f = 4d/(1+2d), copy-neutral LOH f = 2d. The LRR gives a second estimate on the cohort's copy scale (one copy calibrated from the male X, three copies from any constitutional trisomy found). An event at or above 0.9 is constitutional.</p>")
        w.append("<h2>Colours</h2><ul><li><span style=\"color:#D55E00\">&#9632;</span> gain</li><li><span style=\"color:#0072B2\">&#9632;</span> loss</li><li><span style=\"color:#CC79A7\">&#9632;</span> copy-neutral loss of heterozygosity</li><li><span style=\"color:#E69F00\">&#9632;</span> run of homozygosity</li>"
                 "<li><span style=\"color:#009E73\">&#9632;</span> BAF at heterozygous calls, the heterozygous rate</li><li><span style=\"color:#444444\">&#9632;</span> LRR per probe</li><li><span style=\"color:#9a9a9a\">&#9632;</span> BAF at homozygous calls</li></ul>")
        w.append("<h2>The karyotype string</h2><p>ISCN-like: the modal number and the sex-chromosome complement, then the constitutional terms in chromosome order (<code>+21</code>, <code>del(1)(p36.32p36.21)</code>, <code>roh(1)(q21.2q24.1)</code>, <code>upd(15)(iso)?</code> for a whole chromosome without heterozygosity), and one line per mosaic with its cell fraction in brackets: <code>mos 47,XY,+12[0.28]/46,XY</code>, <code>mos 46,XY,del(7)(q11.21q36.3)[0.23]/46,XY</code>, <code>mos 45,X[0.51]/46,XY</code> for a mosaic loss of Y.</p>")
        w.append("<h2>Calling rules</h2><ul><li>Windows of about 100 probes give three tracks per sample: the median corrected LRR, the band split d at heterozygous calls, the heterozygous rate against the Hardy-Weinberg expectation. They are segmented jointly (circular binary segmentation) and every boundary refined to probe resolution.</li>"
                 "<li>A segment is a run of homozygosity when its heterozygous rate sits at the error floor (0.1 of expectation); with LRR at one copy it is a constitutional deletion.</li><li>A gain or loss needs an LRR shift of at least 0.03 and 5 standard errors, with the band split too unless the shift reaches 0.10.</li><li>A copy-neutral LOH is a flat LRR with a band split (d at least 0.02, log-likelihood ratio at least 30).</li>"
                 "<li>Floors: 2 Mb and 40 probes. A segment covering 90% of a chromosome's probed span is whole, of an arm's p or q, else a stretch.</li><li>Runs of homozygosity from the genotypes with PLINK's sliding window (50 probes, one heterozygous call allowed; runs of 1.5 Mb and 100 probes or more); runs of 10 Mb or more are events.</li>"
                 "<li>Artefact-prone regions (the MHC, the 8p23.1 and 17q21.31 inversions, the centromere flanks, the chromosome ends, acrocentric short arms) tag an event; they never remove it.</li></ul>")
        w.append("<h2>Every column</h2>")
        for title, cols in COLUMNS:
            w.append("<h3>%s</h3>" % esc(title) + table(["column", "meaning"], [dict(column=c, meaning=m) for c, m in cols], sortable=False))
        if self.qc_summary:
            w.append("<h2>Thresholds and their sources</h2>" + table(["name", "value", "source"], [dict(name=k, value=v["value"], source=v["source"]) for k, v in self.qc_summary["thresholds"].items()], sortable=False))
        w.append("</body></html>")
        with open(os.path.join(self.out, "report", "guide.html"), "w") as fh:
            fh.write("\n".join(w))

    # ------------------------------------------------------------------------------------------------------ the digest

    def write_digest(self, top=25, min_mb=10.0, min_f=0.1, min_mb_roh=20.0, sex_mosaics=True):
        listed, left = [], []
        for e in self.events:
            f = num(e.get("f", ""))
            size_mb = int(e["size_bp"]) / 1e6
            whole = e["span"] == "whole"
            hard = self.qc.get(e["sample"], {}).get("hard_fail") == "1"
            why = ""
            if hard:
                why = "the sample fails QC (not to be read)"
            elif e["artefact"] != NA and e["type"] in ("gain", "loss") and not finite(num(e["f_baf"])):
                why = "a copy change read from the LRR alone in an artefact-prone or recurrent region (%s)" % e["artefact"]
            elif e["type"] in ("gain", "loss") and e["source"] == "segmentation" and finite(f) and f < 0.9 and not finite(num(e["f_baf"])):
                why = "a mosaic copy change read from the LRR alone (no BAF split to confirm it)"
            elif e["source"] == "sex" and finite(f) and f < 0.9 and not sex_mosaics:
                why = "an X or Y lost or gained in part of the cells (culture or age: --sex-mosaics lists them)"
            elif e["type"] == "roh" and not whole and size_mb < min_mb_roh:
                why = "a run of homozygosity under %g Mb" % min_mb_roh
            elif e["span"] == "stretch" and size_mb < min_mb:
                why = "a stretch under %g Mb" % min_mb
            elif not whole and finite(f) and f < min_f:
                why = "in under %.0f%% of the cells" % (100 * min_f)
            (left if why else listed).append((e, why))
        listed.sort(key=lambda t: rank_key(t[0]))
        shown, beyond = listed[:top], listed[top:]
        rows = [dict(rank=i + 1, sample=e["sample"], event=what(e), chrom=e["chrom"], start=e["start"], end=e["end"], size_mb="%.1f" % (int(e["size_bp"]) / 1e6), span=e["span"], type=e["type"], f=e["f"], f_lrr=e["f_lrr"], f_baf=e["f_baf"],
                     karyotype=self.kar.get(e["sample"], {}).get("karyotype", NA), note=e["note"], left_out="") for i, (e, _) in enumerate(shown)]
        cols = ["rank", "sample", "event", "chrom", "start", "end", "size_mb", "span", "type", "f", "f_lrr", "f_baf", "karyotype", "note"]
        write_tsv(os.path.join(self.out, "report", "digest.tsv"), cols, rows)
        write_tsv(os.path.join(self.out, "report", "digest_left_out.tsv"), ["left_out"] + cols[1:], [dict(left_out=why, sample=e["sample"], event=what(e), chrom=e["chrom"], start=e["start"], end=e["end"], size_mb="%.1f" % (int(e["size_bp"]) / 1e6), span=e["span"], type=e["type"], f=e["f"], f_lrr=e["f_lrr"], f_baf=e["f_baf"], karyotype=self.kar.get(e["sample"], {}).get("karyotype", NA), note=e["note"]) for e, why in left + [(e, "beyond the first %d" % top) for e, _ in beyond]])
        w = ['<!doctype html><html><head><meta charset="utf-8"><title>arraykaryo digest</title><style>%s</style></head><body>' % CSS]
        if self.mock:
            w.append('<div class="mock">%s</div>' % esc(self.mock))
        w.append('<p class="nav"><a href="index.html">the cohort</a><a href="guide.html">the guide</a></p><h1>The largest chromosomal events in %d samples</h1>' % self.C.n_samples)
        w.append('<p class="muted">arraykaryo %s digest, %s. Every event is a reading of the array\'s LRR, BAF and genotypes, to be confirmed by a clinical karyotype or another assay before it is reported.</p>' % (__version__, time.strftime("%Y-%m-%d")))
        w.append('<div class="box"><p><b>What is on this page.</b> The %d largest of %d events, ranked by impact: a whole chromosome first (an aneuploidy, an isodisomy), then an arm, then a stretch by size; within each, an event in every cell before one in part of the cells, and a change of copy number before a copy-neutral LOH before a run of homozygosity. '
                 'Left out and counted in digest_left_out.tsv: stretches under %g Mb, mosaic segments in under %.0f%% of the cells, runs of homozygosity under %g Mb, the calls of samples that fail QC, mosaic copy changes read from the LRR alone (no BAF split to confirm them), copy changes read from the LRR alone in an artefact-prone region%s (the centromere flanks, the chromosome ends, the MHC and the inversions: tagged in events.all.tsv, never removed there). A whole chromosome is never left out for its share of cells.</p>'
                 '<p><b>How to read a figure.</b> <b>LRR</b> is the probe intensity against the cohort on a log2 scale: 0 is two copies, about -0.5 one, about +0.35 three; a mosaic sits in between in proportion to the share of cells. <b>BAF</b> is the B-allele share at the sample\'s heterozygous calls: one band at 0.5 means two balanced copies; two bands mean an imbalance (a gain, a loss, or a copy-neutral LOH when the LRR stays at 0); no heterozygous calls mean one copy or two identical copies.</p></div>' % (
                     len(shown), len(self.events), min_mb, 100 * min_f, min_mb_roh, "" if sex_mosaics else ", and an X or Y lost or gained in part of the cells (the commonest change of a cultured line, and of ageing blood; --sex-mosaics lists them)"))
        w.append("<h2>The events</h2>" + table(["rank", "sample", "event", "size_mb", "f", "karyotype", "note"], rows, sortable=False, link=("sample", lambda r: self.sample_href(r["sample"]))) if rows else "<p class=\"muted\">no event passes the floors</p>")
        if beyond:
            w.append('<p class="muted">%d more events pass the floors but lie beyond the first %d; %d are left out by the floors (digest_left_out.tsv).</p>' % (len(beyond), top, len(left)))
        genome_shown = {}
        for i, (e, _) in enumerate(shown):
            s = e["sample"]
            q = self.qc.get(s, {})
            w.append('<div class="ev" id="ev%d"><h2>%d. %s: %s<span class="tag">%s</span></h2>' % (i + 1, i + 1, esc(s), esc(what(e)), esc({0: "a whole chromosome", 1: "an arm", 2: "a stretch"}[TIER.get(e["span"], 2)] + (", in every cell" if (finite(num(e["f"])) and num(e["f"]) >= 0.9) else ", mosaic"))))
            facts = [("cells carrying it", "%s (by LRR %s, by BAF %s)" % (e["f"], e["f_lrr"], e["f_baf"])), ("LRR over the event", "%s (SE %s) over %s probes" % (e["lrr"], e["lrr_se"], e["n_probes"])),
                     ("band split", "d %s (log-likelihood ratio %s) over %s heterozygous calls; heterozygous ratio %s" % (e["d"], e["d_llr"], e["n_het"], e["het_ratio"])), ("bands", e["bands"]),
                     ("the sample", "%s; QC %s%s" % (self.kar.get(s, {}).get("karyotype", NA), "passes" if q.get("pass") == "1" else "flagged: " + q.get("flags", NA), ("; sex check: " + self.sex[s]["check"]) if s in self.sex and self.sex[s]["check"] not in ("agrees", NA) else "")),
                     ("note", (e["note"] if e["note"] != NA else "-") + (("; artefact-prone region: " + e["artefact"]) if e["artefact"] != NA else ""))]
            w.append('<table class="facts">' + "".join("<tr><td>%s</td><td>%s</td></tr>" % (esc(k), esc(v)) for k, v in facts) + "</table>")
            fdir = os.path.join(self.out, "events", "figures", s)
            if s in genome_shown:
                w.append('<p class="muted">Genome-wide figure of this sample: see event <a href="#ev%d">%d</a>.</p>' % (genome_shown[s], genome_shown[s]))
            else:
                genome_shown[s] = i + 1
                w.append(figure(os.path.join(fdir, "genome.png"), "%s, every chromosome. The event: %s." % (esc(s), esc(what(e))), max_width=1300))
            p = os.path.join(fdir, "chr_%s.png" % e["chrom"]) if e["source"] != "sex" else os.path.join(self.out, "sex", "figures", "tracks_%s.png" % s)
            w.append(figure(p, "%s of %s." % (esc(e["chrom"]), esc(s)), max_width=1300))
            w.append("</div>")
        w.append('<h2>Left out</h2><p class="muted">%s</p></body></html>' % esc("; ".join("%d %s" % (n, why) for why, n in sorted({why: sum(1 for _, x in left if x == why) for _, why in left}.items(), key=lambda kv: -kv[1])) or "nothing"))
        with open(os.path.join(self.out, "report", "digest.html"), "w") as fh:
            fh.write("\n".join(w))
        return len(shown), len(left)


COLUMNS = [
    ("qc/samples.tsv", [("call_rate", "called genotypes over the autosomal SNV probes kept for genotypes"), ("het_rate", "heterozygous share of the called autosomal genotypes"), ("lrr_mean", "mean raw autosomal LRR"),
                        ("lrr_sd_raw / lrr_sd", "SD of the autosomal LRR before / after the GC-wave and reference corrections"), ("baf_sd", "SD of BAF - 0.5 at heterozygous calls, the median over chromosomes"),
                        ("baf_drift", "share of homozygous-called probes with BAF in 0.2-0.25 or 0.75-0.8"), ("wf", "PennCNV's waviness factor: the median absolute deviation of the 1-Mb window LRR medians, signed by the GC correlation"),
                        ("sex / sex_check", "the sex the complement implies (M with a Y, F without) and its check against the reported sex"), ("n_events", "events called on the sample"), ("pass / hard_fail / flags", "no flag / not to be read / the flags (call_rate, lrr_sd, lrr_sd_fail, baf_sd, baf_drift, wave, lrr_mean, sex_mismatch, duplicate, unexpected_relative, event_count)")]),
    ("sex/sex.tsv", [("x_lrr / y_lrr", "the median corrected LRR on the X outside the PARs / over the Y-specific probes"), ("x_het_rel", "X heterozygosity relative to the autosomal rate (0 one X, 1 two X's)"), ("x_copies / y_copies", "copies on the cohort's own scale"),
                     ("base / complement", "the complement the Y implies (XX or XY) / the constitutional reading (XX, XY, X, XXY, XXX, XYY ...)"), ("karyotype / mosaic", "the sex-chromosome karyotype string / the mosaic in words"), ("f / f_lrr / f_baf", "the cell fraction of the mosaic, from the level and from the BAF split"), ("check", "against the reported sex: agrees, or what a disagreement means")]),
    ("pedigree/relatedness.tsv and families.tsv", [("kinship", "KING-robust kinship (0.5 identical, 0.25 first degree, 0.125 second, 0.0625 third)"), ("ibs0", "share of probes at which the pair is opposite homozygous (near 0 for parent-offspring)"), ("class", "duplicate or monozygotic; parent-offspring; full sibling; second degree; third degree"),
                                                  ("kind", "trio; duo; parent-offspring pair unoriented; full siblings; duplicate or monozygotic twins"), ("mendelian_error_rate", "the median over chromosomes of the Mendelian-error rate in the trio's orientation (genome-wide and the other orientations beside)")]),
    ("events/events.all.tsv", [("span", "whole (90% of the chromosome's probed span), p or q (90% of an arm), stretch"), ("type", "gain, loss, loh (copy-neutral loss of heterozygosity), roh (run of homozygosity)"), ("f / f_lrr / f_baf", "the cell fraction read (BAF first) / from the LRR on the cohort's copy scale / from the band split"),
                               ("constitutional", "f at or above 0.9"), ("lrr / lrr_se", "the segment's median corrected LRR and its standard error"), ("d / d_llr", "the band split at heterozygous calls and its log-likelihood ratio against none"), ("het_ratio", "heterozygous rate over the Hardy-Weinberg expectation"),
                               ("source", "segmentation; roh_scan (a run of 10 Mb or more); sex (the sex step's whole-chromosome call)"), ("artefact", "artefact-prone regions the event mostly lies in"), ("note", "the evidence in words and every caveat")]),
    ("events/roh.all.tsv and karyotypes.tsv", [("size_bp / n_probes / n_het", "the run's span, probes and heterozygous calls (errors allowed)"), ("is_event", "10 Mb or more"), ("f_roh", "the genome fraction in runs of homozygosity (consanguinity: a parental first-cousin union gives about 0.06)"), ("karyotype", "the ISCN-like string")]),
]


def _f(v, nd=2):
    return ("%.*f" % (nd, v)) if (v is not None and finite(v)) else NA


def fig_froh(out, kar):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .plots import MM, PAL, save
    v = np.array([num(r["f_roh"]) for r in kar.values()])
    v = v[np.isfinite(v)]
    fig, ax = plt.subplots(figsize=(120 * MM, 70 * MM))
    if v.size:
        ax.hist(v, bins=30, color=PAL["roh"], edgecolor="white", linewidth=0.3)
        ax.axvline(np.median(v), color=PAL["fit"], lw=0.8, label="median %.4f" % np.median(v))
        ax.axvline(0.0625, color=PAL["ref"], lw=0.8, ls=":", label="0.0625: a first-cousin union")
        ax.legend(frameon=False)
    ax.set_xlabel("F_ROH: genome fraction in runs of homozygosity")
    ax.set_ylabel("samples")
    ax.set_title("Runs of homozygosity per sample")
    return save(fig, os.path.join(out, "report", "figures", "froh"), "F_ROH per sample", "The genome fraction in runs of homozygosity of 1.5 Mb or more per sample; the dotted line is the expectation for the child of first cousins.", [("orange", "samples")])


def run_report(out, cache_dir=None, pages="flagged", top=25, digest_min_mb=10.0, digest_min_f=0.1, sex_mosaics=False, threads=1, log_fn=log):
    """pages: 'flagged' (flagged samples, samples with an event, samples not plainly XX or XY), 'digest' (the digest's samples, the
    sex disagreements and the hard QC failures: a cohort of thousands where most samples carry some call stays readable), or
    'all'. Figures missing for a page's sample are drawn here."""
    t0 = time.time()
    R = Report(out, cache_dir)
    rdir = os.path.join(out, "report")
    os.makedirs(os.path.join(rdir, "samples"), exist_ok=True)
    os.makedirs(os.path.join(rdir, "figures"), exist_ok=True)
    from .plots import draw_sample_figures, fig_landscape
    if R.events:
        fig_landscape(out, R.events, R.G, R.C.samples)
    if R.kar:
        fig_froh(out, R.kar)
    shown, left = R.write_digest(top=top, min_mb=digest_min_mb, min_f=digest_min_f, sex_mosaics=sex_mosaics)
    odd_sex = {s for s, r in R.sex.items() if r["mosaic"] != NA or r["check"] not in ("agrees", NA, "not read")}
    if pages == "all":
        want = set(R.C.samples)
    elif pages == "digest":
        digest_samples = {r["sample"] for r in read_tsv(os.path.join(rdir, "digest.tsv"))}
        want = digest_samples | {s for s, q in R.qc.items() if q.get("hard_fail") == "1"} | {s for s, r in R.sex.items() if r["check"] not in ("agrees", NA, "not read", "")}
    else:
        want = {e["sample"] for e in R.events} | {s for s, q in R.qc.items() if q.get("flags", NA) != NA} | odd_sex
    R.pages = {s for s in R.C.samples if s in want}
    for f in os.listdir(os.path.join(rdir, "samples")):  # the pages of an earlier run with another --pages choice
        if f.endswith(".html") and f[:-5] not in R.pages:
            os.remove(os.path.join(rdir, "samples", f))
    missing = [s for s in sorted(R.pages) if s in {e["sample"] for e in R.events} and not os.path.exists(os.path.join(out, "events", "figures", s, "genome.png"))]
    if missing and os.path.exists(os.path.join(out, "events", "windows")):
        n = draw_sample_figures(out, R.cache_dir, missing, threads, log_fn)
        log_fn("figures drawn for %d samples of the pages that had none" % n, "report")
    for s in sorted(R.pages):
        R.write_sample_page(s)
    R.write_index()
    R.write_guide()
    shown, left = R.write_digest(top=top, min_mb=digest_min_mb, min_f=digest_min_f, sex_mosaics=sex_mosaics)
    log_fn("report/index.html, %d sample pages, the guide, the digest (%d events listed, %d left out) in %.1f s" % (len(R.pages), shown, left, time.time() - t0), "report")
    mark_done(out, "report", dict(seconds=time.time() - t0, pages=len(R.pages)))
    return dict(pages=len(R.pages), digest=shown)
