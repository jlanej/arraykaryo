"""arraykaryo command line.

  arraykaryo inspect   cohort.bcf [--json out.json]            what the export carries; refuses one without GT/LRR/BAF
  arraykaryo mock      --out DIR [--samples 60 --seed 1 ...]   a simulated export with families, planted events and QC failures, with truth.json
  arraykaryo extract   --bcf cohort.bcf --out DIR              one pass into per-chromosome arrays (GT int8, LRR float32, BAF float16), a probe table, a manifest
  arraykaryo qc        --out DIR [--sample-sheet ...]          probe QC, sample QC, GC-wave correction, the cohort reference
  arraykaryo sex       --out DIR                               the sex-chromosome complement per sample, mosaics included, against the reported sex
  arraykaryo pedigree  --out DIR [--pedigree ...]              relatedness, families, the check against a given pedigree
  arraykaryo events    --out DIR                               near-chromosomal events and runs of homozygosity per sample; karyotype strings; figures
  arraykaryo report    --out DIR                               the cohort page, per-sample pages, the guide, the digest
  arraykaryo run       --bcf cohort.bcf --out DIR              all of the above

Every command reads and writes under one output directory (--out; the cache is DIR/cache), skips what is done unless --redo,
and logs what it did and how long it took."""
import argparse
import os
import sys

from . import __version__
from .util import log

TODO = {
    "report": "the cohort report and the per-sample pages (PROMPT.md, section 4.6)",
    "run": "extract, qc, sex, pedigree, events and report in one go (PROMPT.md, section 4.8)",
}


def cmd_inspect(a):
    from .inspect import describe, inspect, write_json
    d = inspect(a.path, records=a.records, samples=None if a.samples <= 0 else a.samples)
    print(describe(d))
    if a.json:
        write_json(d, a.json)
        log("-> " + a.json)
    return 2 if d["required_missing"] else 0


def cmd_mock(a):
    from .mock import build
    build(a.out, n_samples=a.samples, seed=a.seed, genome=a.genome, contigs=a.contigs, full=a.full, density_kb=a.density_kb, trios=a.trios, duos=a.duos,
          events=not a.no_events, families=not a.no_families, qc_failures=not a.no_qc_failures, noise=a.noise, cell_fraction=a.cell_fraction)
    return 0


def cache_dir(a):
    return a.cache if getattr(a, "cache", "") else os.path.join(a.out, "cache")


def cmd_extract(a):
    from .extract import extract
    samples = None
    if a.samples:
        with open(a.samples) as fh:
            samples = [l.strip().split()[0] for l in fh if l.strip()]
    chroms = [c for c in a.chroms.split(",") if c] if a.chroms else None
    extract(a.bcf, cache_dir(a), threads=a.threads, route=a.route, baf_dtype=a.baf_dtype, redo=a.redo, samples=samples, genome=a.genome, chroms=chroms)
    return 0


def cmd_qc(a):
    from .qc import run_qc
    from .util import is_done
    if is_done(a.out, "qc") and not a.redo:
        log("qc is done (--redo to repeat)", "qc")
        return 0
    run_qc(a.out, cache_dir=cache_dir(a), sample_sheet=a.sample_sheet, exclude_probes=a.exclude_probes, exclude_samples=a.exclude_samples, threads=a.threads)
    return 0


def cmd_sex(a):
    from .sex import run_sex
    from .util import is_done
    if is_done(a.out, "sex") and not a.redo:
        log("sex is done (--redo to repeat)", "sex")
        return 0
    run_sex(a.out, cache_dir=cache_dir(a), threads=a.threads)
    return 0


def cmd_pedigree(a):
    from .pedigree import run_pedigree
    from .util import is_done
    if is_done(a.out, "pedigree") and not a.redo:
        log("pedigree is done (--redo to repeat)", "pedigree")
        return 0
    run_pedigree(a.out, cache_dir=cache_dir(a), pedigree=a.pedigree, sample_sheet=a.sample_sheet, threads=a.threads)
    return 0


def cmd_events(a):
    from .events import run_events
    from .util import is_done
    if is_done(a.out, "events") and not a.redo:
        log("events is done (--redo to repeat)", "events")
        return 0
    params = dict(min_mb=a.min_mb, window=a.window)
    run_events(a.out, cache_dir=cache_dir(a), threads=a.threads, params=params, figures=a.figures)
    return 0


def cmd_todo(a):
    sys.exit("arraykaryo %s is not implemented yet: %s" % (a.cmd, TODO[a.cmd]))


def add_common(p, bcf=False):
    p.add_argument("--out", required=True, help="the output directory (the cache lives in DIR/cache unless --cache names another)")
    p.add_argument("--cache", default="", help="the cache directory (default DIR/cache)")
    p.add_argument("--redo", action="store_true", help="redo the step even if it is marked done")
    p.add_argument("--threads", type=int, default=1, help="chromosomes (extract) or samples (events, figures) in parallel")
    if bcf:
        p.add_argument("--bcf", required=True, help="the indexed VCF/BCF export (FORMAT GT, LRR, BAF)")
        p.add_argument("--genome", default="", help="the assembly (grch37, grch38, chm13) or a band table, when the contig lengths do not identify it")
        p.add_argument("--samples", default="", help="a file of sample names to restrict to")
        p.add_argument("--route", default="auto", choices=["auto", "bcftools", "cyvcf2"], help="the extraction route (auto: cyvcf2 when installed, else bcftools)")
        p.add_argument("--baf-dtype", default="float16", choices=["float16", "float32"], help="the BAF storage type (default float16; METHODS records the precision check)")
        p.add_argument("--chroms", default="", help="a subset of chromosomes (chr21,chrX) for a quick pass")


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
    m = sub.add_parser("mock", help="a simulated export with families, planted events and QC failures, and truth.json")
    m.add_argument("--out", required=True)
    m.add_argument("--samples", type=int, default=60)
    m.add_argument("--seed", type=int, default=1)
    m.add_argument("--genome", default="chm13", help="chm13 (default), grch38 or grch37")
    m.add_argument("--contigs", default="chr1,chr7,chr11,chr12,chr15,chr21,chrX,chrY", help="the chromosomes written (default: the ones the planted events need)")
    m.add_argument("--full", action="store_true", help="all 24 chromosomes")
    m.add_argument("--density-kb", type=float, default=12.0, help="mean probe spacing in kb (default 12; the example array is about 1.3)")
    m.add_argument("--trios", type=int, default=4)
    m.add_argument("--duos", type=int, default=2)
    m.add_argument("--no-events", action="store_true"), m.add_argument("--no-families", action="store_true"), m.add_argument("--no-qc-failures", action="store_true")
    m.add_argument("--noise", type=float, default=1.0, help="scales every sample's LRR and BAF noise")
    m.add_argument("--cell-fraction", type=float, default=1.0, help="scales the mosaic cell fractions")
    m.set_defaults(func=cmd_mock)
    e = sub.add_parser("extract", help="one pass over the export into per-chromosome arrays (GT, LRR, BAF), a probe table and a manifest")
    add_common(e, bcf=True)
    e.set_defaults(func=cmd_extract)
    q = sub.add_parser("qc", help="probe QC, sample QC, the GC-wave correction and the cohort reference; qc/ tables and figures")
    add_common(q)
    q.add_argument("--sample-sheet", default="", help="sample, reported_sex, batch (by header name); other columns carried to the report")
    q.add_argument("--exclude-probes", default="", help="a file of probe ids to drop")
    q.add_argument("--exclude-samples", default="", help="a file of samples to mark not to be read")
    q.set_defaults(func=cmd_qc)
    x = sub.add_parser("sex", help="the sex-chromosome complement per sample, mosaics included, against the reported sex; sex/ tables and figures")
    add_common(x)
    x.set_defaults(func=cmd_sex)
    pd = sub.add_parser("pedigree", help="relatedness from the genotypes (KING-robust kinship, IBS0), the families it implies, the check against a given pedigree")
    add_common(pd)
    pd.add_argument("--pedigree", default="", help="a pedigree to check: the house form (#kid dad mom kid_sex dad_sex mom_sex), a PLINK .fam, or any table naming the ids")
    pd.add_argument("--sample-sheet", default="", help="the sample sheet (an age column orients duos)")
    pd.set_defaults(func=cmd_pedigree)
    ev = sub.add_parser("events", help="near-chromosomal events and runs of homozygosity per sample, with cell fractions; karyotype strings; figures")
    add_common(ev)
    ev.add_argument("--min-mb", type=float, default=2.0, help="the smallest segment reported, in Mb (default 2)")
    ev.add_argument("--window", type=int, default=100, help="probes per window (default 100)")
    ev.add_argument("--figures", default="flagged", choices=["flagged", "all", "none"], help="figures for samples with an event or a flag (default), for all, or none")
    ev.set_defaults(func=cmd_events)
    for name, what in TODO.items():
        t = sub.add_parser(name, help="(not yet implemented) " + what)
        add_common(t, bcf=(name == "run"))
        t.set_defaults(func=cmd_todo)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
