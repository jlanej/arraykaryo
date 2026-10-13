# arraykaryo

Large chromosomal events in a genotyping-array cohort, from the array's VCF/BCF export alone (FORMAT `GT`, `LRR`, `BAF`):
array QC, a pedigree inferred from the genotypes, the sex-chromosome complement of every sample with its aneuploidies
and mosaics, and the near-chromosomal events - whole chromosomes, arms and long segments gained, lost or made
homozygous, in every cell or in part of them - with runs of homozygosity, each with its cell fraction; one
self-contained report. The array sibling of [triokaryo](https://github.com/jlanej/triokaryo), which reads the same
events from whole-genome sequencing of parent-offspring trios.

**Status: 1.0.0 - every milestone of [PROMPT.md](PROMPT.md).** The whole chain works and is tested on the mock (every planted
event, family, complement and QC failure recovered; no event on a null sample), its detection curve is measured at the
example array's probe spacing ([docs/calibration](docs/calibration/README.md)), and the public example - a 1000 Genomes
export of 2,141 lymphoblastoid lines - is run, validated against its known pedigree and sexes and against the events
triokaryo read from sequencing, and documented with the platform's caveats ([docs/example](docs/example/README.md)). Methods and every threshold with its source:
[docs/METHODS.md](docs/METHODS.md); every file and column: [docs/OUTPUT.md](docs/OUTPUT.md); timings: [benchmarks/](benchmarks/README.md).

## What goes in

One VCF or BCF, indexed, with every sample of the cohort and, per call, `FORMAT/GT`, `FORMAT/LRR` and `FORMAT/BAF`, as
`bcftools +gtc2vcf` writes an Illumina export (other fields are ignored). Any assembly: GRCh37, GRCh38 and T2T-CHM13
are recognised from the contig lengths. Optionally a sample sheet (reported sex, batch) and a pedigree in the
`#kid dad mom kid_sex dad_sex mom_sex` form; both are checked against what the genotypes say, never trusted over them.

```bash
arraykaryo run --bcf cohort.bcf --out out --threads 8 --sample-sheet samples.tsv --pedigree trios.tsv   # everything below, one log
#   add --pages digest on a cohort of thousands: pages and figures only for the digest's samples, the sex disagreements and the hard failures
#   add --sex-mosaics to list an X or Y lost or gained in part of the cells on the digest (left out by default: cultured lines, ageing blood)
arraykaryo inspect cohort.bcf                                   # samples, contigs and assembly, fields, the conventions seen; refuses an export without GT/LRR/BAF
arraykaryo extract --bcf cohort.bcf --out out --threads 8       # one pass into out/cache: per-chromosome arrays of GT, LRR and BAF, a probe table, a manifest
#   add --fasta ref.fa (indexed) when the export carries no INFO/GC: the probe GC is computed from the reference (gtc2vcf's window) for the GC-wave correction
arraykaryo qc --out out --sample-sheet samples.tsv --threads 8  # probe QC, sample QC (pass or flagged with the reason), the GC-wave correction, the cohort reference
arraykaryo sex --out out                                        # the sex-chromosome complement of every sample, mosaics included, against the reported sex
arraykaryo pedigree --out out --pedigree trios.tsv              # kinship and IBS0 of every pair, the families they imply, the check against the given pedigree
arraykaryo events --out out --threads 8                         # near-chromosomal events and runs of homozygosity per sample, cell fractions, karyotype strings, figures
arraykaryo report --out out                                     # the cohort page, per-sample pages, the guide, the digest (self-contained HTML)
arraykaryo mock --out mock_out                                  # a simulated export with families, planted events and QC failures, and truth.json
```

## The example, on one page

[jlanej.github.io/arraykaryo](https://jlanej.github.io/arraykaryo/) shows what one run found on the public 1000 Genomes
cohort: every event on the chromosomes, the mosaics checked against sequencing, the sex chromosomes of every sample, the
families, the digest. It is `docs/index.html`, built from the run's tables by `scripts/build_site.py RUN_DIR docs` (GitHub
Pages serves `docs/`).

## What comes out

`out/cache/` (the extraction), `out/qc/` (probes.tsv, samples.tsv, the correction, figures), `out/sex/` (sex.tsv, the X-against-Y
figure, X and Y tracks of every sample not plainly XX or XY), `out/pedigree/` (relatedness.tsv, families.tsv, the inferred
pedigree, check.tsv, the kinship-against-IBS0 figure), `out/events/` (events.all.tsv, roh.all.tsv, karyotypes.tsv, the
genome-wide and per-chromosome figures of every sample with an event or a flag), `out/report/` (index.html, samples/, guide.html,
digest.html). Every file and
column: [docs/OUTPUT.md](docs/OUTPUT.md).

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

A public 1000 Genomes array export (2,141 samples, 2.4 million probes, CHM13 coordinates) is the example and the validation
set: 331 known trios and 47 duos recovered from the genotypes alone with 779 agreeing relationships and no claimed parent
unrelated; 2,126 of 2,130 reported sexes agree, the four that do not being two XX "males" (a swap or an XX male) and two
47,XXY, with the lymphoblastoid lines' partial losses of Y (12 men) and of an X (29 women) listed with their cell fractions;
HG01103's mosaic losses of 2q and 14q and NA12748's mosaic gain of 13q and terminal loss of 21q found where whole-genome
sequencing put them. The platform's own numbers - a one-copy level of -0.40, a baseline band split of 0.01-0.03, mismapped
probes on CHM13's heterochromatin and acrocentric arms - are measured and written down: [docs/example/README.md](docs/example/README.md).
No data and no data locations live in this repository; the example is pointed at by environment variables.

## Conventions

- No real data in the repository, ever; the mock writes what the tests need, the example is public and fetched by its owner.
- No local or remote paths in code, docs or tests; inputs come from arguments and environment variables.
- Figures carry no legend or title inside the image: a text sidecar and an editable legend image go beside each.
- Reports are one self-contained HTML file per cohort (figures embedded), readable without a server.
- MIT licence.
