"""arraykaryo command line.

  arraykaryo inspect cohort.bcf [--records 500] [--samples 200] [--json out.json]   what the export carries; refuses one without GT/LRR/BAF
  arraykaryo extract | qc | sex | pedigree | events | report | mock | run           specified in PROMPT.md; not yet implemented
"""
import argparse
import sys

from . import __version__

TODO = {
    "extract": "one pass over the export into per-chromosome arrays of GT, LRR and BAF (PROMPT.md, section 4.1)",
    "qc": "probe and sample QC, GC-wave correction, the cohort reference (PROMPT.md, section 4.2)",
    "sex": "the sex-chromosome complement per sample, aneuploidies and mosaics included (PROMPT.md, section 4.3)",
    "pedigree": "relatedness from the genotypes and the pedigree it implies, against any pedigree given (PROMPT.md, section 4.4)",
    "events": "near-chromosomal gains, losses, copy-neutral LOH and runs of homozygosity, mosaic or not (PROMPT.md, section 4.5)",
    "report": "the cohort report and the per-sample pages (PROMPT.md, section 4.6)",
    "mock": "a simulated export with planted events, families and QC failures, with its truth (PROMPT.md, section 4.7)",
    "run": "extract, qc, sex, pedigree, events and report in one go (PROMPT.md, section 4.8)",
}


def _log(msg):
    print("[arraykaryo] " + msg, file=sys.stderr, flush=True)


def cmd_inspect(a):
    from .inspect import describe, inspect, write_json
    d = inspect(a.path, records=a.records, samples=None if a.samples <= 0 else a.samples)
    print(describe(d))
    if a.json:
        write_json(d, a.json)
        _log("-> " + a.json)
    return 2 if d["required_missing"] else 0


def cmd_todo(a):
    sys.exit("arraykaryo %s is not implemented in this starting point: %s" % (a.cmd, TODO[a.cmd]))


def build_parser():
    p = argparse.ArgumentParser(prog="arraykaryo", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version="arraykaryo %s" % __version__)
    sub = p.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("inspect", help="what an array VCF/BCF export carries (samples, contigs, assembly, fields, conventions); refuses one without GT/LRR/BAF")
    i.add_argument("path")
    i.add_argument("--records", type=int, default=500, help="records read, in file order, for the conventions (default 500)")
    i.add_argument("--samples", type=int, default=200, help="samples read for the conventions (default 200; 0 for all)")
    i.add_argument("--json", default="", help="also write the inspection as JSON")
    i.set_defaults(func=cmd_inspect)
    for name, what in TODO.items():
        t = sub.add_parser(name, help="(not yet implemented) " + what)
        t.set_defaults(func=cmd_todo)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
