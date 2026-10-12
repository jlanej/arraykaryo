# arraykaryo

Large chromosomal events in a genotyping-array cohort, from the array's VCF/BCF export alone (FORMAT `GT`, `LRR`, `BAF`):
array QC, a pedigree inferred from the genotypes, the sex-chromosome complement of every sample with its aneuploidies
and mosaics, and the near-chromosomal events - whole chromosomes, arms and long segments gained, lost or made
homozygous, in every cell or in part of them - with runs of homozygosity, each with its cell fraction; one
self-contained report. The array sibling of [triokaryo](https://github.com/jlanej/triokaryo), which reads the same
events from whole-genome sequencing of parent-offspring trios.

**Status: a starting point.** `arraykaryo inspect` works and gates the input; everything else is specified in
[PROMPT.md](PROMPT.md), the implementation brief, and is to be built from it. The brief is the contract; this README
will follow the implementation.

## What goes in

One VCF or BCF, indexed, with every sample of the cohort and, per call, `FORMAT/GT`, `FORMAT/LRR` and `FORMAT/BAF`, as
`bcftools +gtc2vcf` writes an Illumina export (other fields are ignored). Any assembly: GRCh37, GRCh38 and T2T-CHM13
are recognised from the contig lengths. Optionally a sample sheet (reported sex, batch) and a pedigree in the
`#kid dad mom kid_sex dad_sex mom_sex` form; both are checked against what the genotypes say, never trusted over them.

```bash
arraykaryo inspect cohort.bcf          # samples, contigs and assembly, fields, the conventions seen; refuses an export without GT/LRR/BAF
```

## What will come out (see PROMPT.md)

`extract` (one pass into per-chromosome arrays), `qc` (probe and sample QC, GC-wave correction, the cohort reference),
`sex`, `pedigree`, `events`, `report`, a `mock` with planted truth for the tests, and `run` for the whole chain.

## Development

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e '.[test]'
pytest -q
```

CI runs the tests on Python 3.9 and 3.12, builds the image and publishes it to GHCR with the tags `main` and
`sha-<commit>` (no `latest`):

```bash
apptainer build arraykaryo.sif docker://ghcr.io/<owner>/arraykaryo:main
apptainer run arraykaryo.sif inspect cohort.bcf
```

## The example

A public 1000 Genomes array export (2,141 samples, about 2.4 million probes) is the example and the validation set:
its known trios, its reported sexes, and the samples whose large events triokaryo read from sequencing. How to run it
and what to check: [docs/example/README.md](docs/example/README.md). No data and no data locations live in this
repository; the example is pointed at by environment variables.

## Conventions

- No real data in the repository, ever; the mock writes what the tests need, the example is public and fetched by its owner.
- No local or remote paths in code, docs or tests; inputs come from arguments and environment variables.
- Figures carry no legend or title inside the image: a text sidecar and an editable legend image go beside each.
- Reports are one self-contained HTML file per cohort (figures embedded), readable without a server.
- MIT licence.
