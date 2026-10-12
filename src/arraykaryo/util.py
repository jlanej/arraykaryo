"""Shared helpers: logging with timing, TSV and JSON I/O in the house conventions (header, booleans 1/0, missing NA), robust
statistics, the done-markers that make every command idempotent, and the readers of a sample sheet and a pedigree (by header
names, never by column position where a header exists)."""
import csv
import json
import math
import os
import sys
import time

import numpy as np

from . import __version__

NA = "NA"


def log(msg, step=""):
    print("[arraykaryo%s] %s" % ((" " + step) if step else "", msg), file=sys.stderr, flush=True)


class Timer:
    """with Timer('extract') as t: ...  -> logs the elapsed time at the end; t.seconds afterwards."""

    def __init__(self, step, what=""):
        self.step, self.what, self.seconds = step, what, 0.0

    def __enter__(self):
        self.t0 = time.time()
        return self

    def __exit__(self, *a):
        self.seconds = time.time() - self.t0
        log("%sdone in %.1f s" % ((self.what + " ") if self.what else "", self.seconds), self.step)


def fmt(x, nd=6):
    """A value for a TSV cell: booleans 1/0, None and NaN as NA, floats with nd significant digits, everything else as str."""
    if x is None:
        return NA
    if isinstance(x, (bool, np.bool_)):
        return "1" if x else "0"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        if math.isnan(x):
            return NA
        if float(x).is_integer() and abs(x) < 1e15:
            return str(int(x))
        return "%.*g" % (nd, x)
    s = str(x)
    return s if s != "" else NA


def write_tsv(path, cols, rows, nd=6):
    """rows: dicts (missing keys -> NA) or sequences in `cols` order."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            if isinstance(r, dict):
                fh.write("\t".join(fmt(r.get(c), nd) for c in cols) + "\n")
            else:
                fh.write("\t".join(fmt(v, nd) for v in r) + "\n")


def read_tsv(path):
    """A TSV with a header -> list of dicts of strings (a leading '#' on the header is dropped)."""
    with open(path) as fh:
        head = fh.readline()
        if not head:
            return []
        cols = head.lstrip("#").rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def num(v, default=np.nan):
    """A TSV cell back to a float (NA, '' -> default)."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def write_json(obj, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_json_default)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def robust_sd(x):
    """1.4826 x the median absolute deviation of the finite values (NaN with fewer than 3)."""
    v = np.asarray(x, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return float("nan")
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def nanmedian(x):
    v = np.asarray(x, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else float("nan")


def finite(x):
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def done_path(out, step):
    return os.path.join(out, step, "done.json")


def is_done(out, step):
    return os.path.exists(done_path(out, step))


def mark_done(out, step, info=None):
    d = dict(step=step, version=__version__, time=time.strftime("%Y-%m-%d %H:%M:%S"))
    if info:
        d.update(info)
    write_json(d, done_path(out, step))


def clear_done(out, step):
    p = done_path(out, step)
    if os.path.exists(p):
        os.remove(p)


def parse_sex(v):
    """M/F/1/2/male/female (any case) -> 'M', 'F'; anything else ''."""
    s = str(v).strip().lower()
    if s in ("m", "1", "male", "man", "xy"):
        return "M"
    if s in ("f", "2", "female", "woman", "xx"):
        return "F"
    return ""


def _split(line):
    line = line.rstrip("\r\n")
    if "\t" in line:
        return line.split("\t")
    if "," in line and " " not in line.strip():
        return line.split(",")
    return line.split()


def _find(cols, names):
    low = {c.strip().lower(): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None


SAMPLE_COL = ("sample", "sample_id", "sampleid", "iid", "id", "individual id", "name")
SEX_COL = ("reported_sex", "sex", "gender", "reported sex")
BATCH_COL = ("batch", "plate", "chip", "run", "sentrixbarcode_a", "sentrix barcode", "array")


def read_sample_sheet(path):
    """A sample sheet (TSV or CSV; header required) -> ({sample: {column: value, 'reported_sex': M/F/'', 'batch': str}}, columns).
    Columns are found by name: the sample in sample/sample_id/IID/ID, the sex in reported_sex/sex/gender, the batch in
    batch/plate/chip/run; every other column is carried through."""
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("##")]
    if not lines:
        raise ValueError("empty sample sheet: " + path)
    cols = [c.lstrip("#").strip() for c in _split(lines[0])]
    sc = _find(cols, SAMPLE_COL)
    if sc is None:
        raise ValueError("the sample sheet %s has no sample column (one of %s)" % (path, ", ".join(SAMPLE_COL)))
    xc, bc = _find(cols, SEX_COL), _find(cols, BATCH_COL)
    out = {}
    for line in lines[1:]:
        vals = _split(line)
        row = dict(zip(cols, vals))
        s = row.get(sc, "").strip()
        if not s:
            continue
        row["reported_sex"] = parse_sex(row.get(xc, "")) if xc else ""
        row["batch"] = row.get(bc, "").strip() if bc else ""
        out[s] = row
    return out, cols


KID_COL = ("kid", "child", "proband", "individual id", "iid", "sample_id", "sample", "id")
DAD_COL = ("dad", "father", "paternal id", "paternal_id", "pat", "fid_father")
MOM_COL = ("mom", "mother", "maternal id", "maternal_id", "mat")
FAM_COL = ("family id", "family_id", "fid", "family")


def read_pedigree(path):
    """A pedigree in the house form (`#kid dad mom kid_sex dad_sex mom_sex`), a PLINK .fam/.ped (the standard six columns, no
    header), or any table with a header naming the individual, the paternal and the maternal id and the sex (the IGSR table,
    peddy's ped) -> (rows, sex_of, family_of): rows = [{kid, dad, mom, kid_sex, dad_sex, mom_sex, family}] for every
    individual with a parent named ('' for an unnamed parent), sex_of = {sample: M/F/''} for every individual named, and
    family_of = {sample: family id or ''}."""
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    if not lines:
        return [], {}, {}
    first = _split(lines[0])
    head = [c.lstrip("#").strip() for c in first]
    low = [c.lower() for c in head]
    if low[:3] == ["kid", "dad", "mom"]:
        mode = "house"
    elif lines[0].startswith("#") or _find(head, DAD_COL) and _find(head, MOM_COL) and _find(head, KID_COL):
        mode = "header"
    else:
        mode = "plink"
    rows, sex_of, family_of = [], {}, {}
    if mode == "plink":
        for line in lines:
            v = _split(line)
            if len(v) < 5:
                continue
            fam, kid, dad, mom, sex = v[0], v[1], v[2], v[3], v[4]
            _add(rows, sex_of, family_of, kid, dad, mom, parse_sex(sex), "", "", fam)
        return rows, sex_of, family_of
    body = lines[1:]
    if mode == "house":
        ix = {c: i for i, c in enumerate(low)}
        for line in body:
            v = _split(line)
            get = lambda c: v[ix[c]] if c in ix and ix[c] < len(v) else ""  # noqa: E731
            _add(rows, sex_of, family_of, get("kid"), get("dad"), get("mom"), parse_sex(get("kid_sex")), parse_sex(get("dad_sex")), parse_sex(get("mom_sex")), get("family"))
        return rows, sex_of, family_of
    kc, dc, mc = _find(head, KID_COL), _find(head, DAD_COL), _find(head, MOM_COL)
    xc, fc = _find(head, ("sex", "gender")), _find(head, FAM_COL)
    if not (kc and dc and mc):
        raise ValueError("pedigree %s: cannot find the individual, paternal and maternal id columns in its header" % path)
    for line in body:
        v = _split(line)
        row = dict(zip(head, v))
        _add(rows, sex_of, family_of, row.get(kc, ""), row.get(dc, ""), row.get(mc, ""), parse_sex(row.get(xc, "")) if xc else "", "", "", row.get(fc, "") if fc else "")
    return rows, sex_of, family_of


def _add(rows, sex_of, family_of, kid, dad, mom, kid_sex, dad_sex, mom_sex, fam):
    dad = "" if dad in ("0", ".", "", "NA") else dad
    mom = "" if mom in ("0", ".", "", "NA") else mom
    kid = kid.strip()
    if not kid:
        return
    if kid_sex or kid not in sex_of:
        sex_of[kid] = kid_sex or sex_of.get(kid, "")
    if fam:
        family_of[kid] = fam
    for p, s, expect in ((dad, dad_sex, "M"), (mom, mom_sex, "F")):
        if p:
            sex_of.setdefault(p, s or expect)
            if fam:
                family_of.setdefault(p, fam)
    if dad or mom:
        rows.append(dict(kid=kid, dad=dad, mom=mom, kid_sex=sex_of.get(kid, ""), dad_sex=sex_of.get(dad, "") if dad else "", mom_sex=sex_of.get(mom, "") if mom else "", family=fam))


def write_pedigree(path, rows):
    """The house form."""
    with open(path, "w") as fh:
        fh.write("#kid\tdad\tmom\tkid_sex\tdad_sex\tmom_sex\n")
        for r in rows:
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (r["kid"], r.get("dad") or "0", r.get("mom") or "0", _sex_code(r.get("kid_sex")), _sex_code(r.get("dad_sex")), _sex_code(r.get("mom_sex"))))


def _sex_code(s):
    return {"M": "1", "F": "2"}.get(s or "", "0")


def chunks(n, size):
    """[(start, stop)] covering range(n) in blocks of `size`."""
    return [(i, min(i + size, n)) for i in range(0, n, size)]


def csv_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))
