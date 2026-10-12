"""The karyotype string per sample in the ISCN-like grammar triokaryo uses: the main line is the modal number and the
sex-chromosome complement plus the constitutional autosomal terms in chromosome order (+21, -18, upd(15)(iso)? for a
whole-chromosome run of homozygosity, roh(1)(q21.3q24.1), loh(11)(p15.5p11.12), dup(13)(q31.2q34), del(1)(p36.32p36.22));
each mosaic event (cell fraction under 0.9) forms its own line relative to the base complement with its cell fraction in
brackets: mos 47,XY,+12[0.30]/46,XY; mos 46,XX,del(7)(q11.21q36.3)[0.25]/46,XX; the sex step's own mosaics (mos 45,X[0.50]/46,XY)
are kept as it wrote them. Segments are given by their cytogenetic bands from the vendored tables, else in Mb."""
from .util import finite

NAME = {"XXY": "47,XXY", "XYY": "47,XYY", "XXX": "47,XXX", "X": "45,X", "XXYY": "48,XXYY", "XXXY": "48,XXXY", "XXXX": "48,XXXX", "XY": "46,XY", "XX": "46,XX"}


def _coords(e):
    if e.get("bands"):
        return e["bands"]
    return "%.1f-%.1fMb" % (e["start"] / 1e6, e["end"] / 1e6)


def karyotype_string(events, complement, sex_karyotype, chrom_order):
    """events: dicts with chrom, start, end, type (gain/loss/loh/roh), span (whole/p/q/stretch), f, bands, constitutional (bool),
    source. complement: the constitutional sex-chromosome complement (XX, XY, XXY ...); sex_karyotype: the sex step's string
    (46,XY or mos 45,X[0.50]/46,XY)."""
    comp = complement or "?"
    base_n = 44 + len(comp) if complement else 46
    autos = [e for e in events if e["chrom"] not in ("chrX", "chrY") or e.get("source") != "sex"]
    autos = sorted(autos, key=lambda e: (chrom_order.get(e["chrom"], 99), e["start"]))
    is_whole_cn = lambda e: e["span"] == "whole" and e["type"] in ("gain", "loss") and e["chrom"] not in ("chrX", "chrY")  # noqa: E731
    modal = base_n + sum((1 if e["type"] == "gain" else -1) for e in autos if is_whole_cn(e) and e.get("constitutional"))
    const, mosaic = [], []
    for e in autos:
        n = e["chrom"][3:]
        mos = not e.get("constitutional")
        doubt = "?" if e.get("doubt") else ""
        if e["type"] == "roh" and e["span"] == "whole":
            term = "upd(%s)(iso)?" % n
        elif e["type"] == "roh":
            term = "roh(%s)(%s)" % (n, _coords(e))
        elif e["type"] == "loh":
            term = "loh(%s)(%s)" % (n, _coords(e))
        elif e["span"] == "whole":
            term = "%s%s" % ("+" if e["type"] == "gain" else "-", n)
        else:
            term = "%s(%s)(%s)" % ("dup" if e["type"] == "gain" else "del", n, _coords(e))
        term += doubt
        if mos and finite(e.get("f")):
            m2 = modal + (1 if (is_whole_cn(e) and e["type"] == "gain") else -1 if is_whole_cn(e) else 0)
            mosaic.append(("%d,%s,%s" % (m2, comp, term), float(e["f"])))
        else:
            const.append(term)
    main = "%d,%s" % (modal, comp) + ("," + ",".join(const) if const else "")
    sex_mos = ""
    if sex_karyotype and sex_karyotype.startswith("mos "):
        sex_mos = sex_karyotype[4:].rsplit("/", 1)[0]                       # the sex step's mosaic line(s), without its base
    elif sex_karyotype and sex_karyotype not in ("NA", "") and not const and not mosaic and len(comp) <= 2:
        pass
    lines = ([sex_mos] if sex_mos else []) + ["%s[%.2f]" % (line, f) for line, f in mosaic]
    if lines:
        return "mos " + "/".join(lines) + "/" + main
    return main
