# arraykaryo: the implementation brief

You are implementing `arraykaryo`, a Python tool that reads a genotyping-array cohort's VCF/BCF export and reports,
for every sample, the array's quality, the sex-chromosome complement, the relatedness and pedigree the genotypes
imply, and the large chromosomal events - whole chromosomes, arms and long segments gained, lost or made homozygous,
in every cell or in part of them, with runs of homozygosity - each with its cell fraction, in one self-contained,
interpretable report. It is the array sibling of triokaryo (whole-genome sequencing of trios); borrow its grammar
(event types, karyotype strings, the plain-words digest, figures without legends) and its discipline (mock first,
every claim tested), not its code.

This file is the contract. Read it whole before writing code. Where it says VERIFY, do not trust the number given:
check the source named and the data.

## 0. How to work

- **Mock first, always.** No real data is ever committed, downloaded into the repository, or read by a test. The
  tests run on `arraykaryo mock`, a simulator you write first (section 4.7), with a `truth.json` every test reads.
  The only real data you may look at is the public 1000 Genomes example (section 9), and only through environment
  variables - never a path in code, docs, tests or commit messages, local or remote.
- **Every sample is called on its own.** QC decides which samples are reliable; every reliable sample receives
  sex, events and ROH calls from its own LRR, BAF and genotypes. Families (trios, duos) are an addition - parent of
  origin, inheritance, Mendelian-error evidence - never a requirement, and a sample outside any family loses nothing
  but those additions.
- **Flags, not drops.** A QC failure marks a sample (and its calls) as not to be read; it does not remove the sample
  from the tables. The cohort reference and the pedigree are built from the passing samples only.
- **Small commits, each with tests green, each with a one-paragraph message that says what changed and why.** Bump
  the version (one place: `pyproject.toml` and `src/arraykaryo/__init__.py`) on every change that alters an output.
- **Write for a scientist who did not run the tool.** Names and sentences in plain words; numbers in tables; every
  column and colour explained in a guide page; no acronym without its expansion the first time.
- **Measure before choosing.** Where this brief offers two routes (section 4.1), time both on the mock at cohort
  scale and keep the faster; record the numbers in `docs/METHODS.md`.
- Dependencies: numpy, pysam, matplotlib; `bcftools` on PATH (and in the image). Add a dependency only when you have
  measured that it earns its place (cyvcf2 for extraction, ruptures for segmentation are candidates), and keep the
  package installable with `pip install -e '.[test]'` on Python 3.9 and 3.12.
- Figures: no legend and no title inside the image; a text sidecar (`<figure>.txt`: suggested title, caption, key)
  and an editable legend image beside each, as triokaryo does.
- Reports: one self-contained HTML file per cohort (figures embedded as base64), readable from a file:// URL.
- Licence MIT. CI publishes the image to GHCR with the tags `main` and `sha-<commit>`; there is no `latest`.

## 1. Mission and the first phase's scope

Phase 1 covers, in this order of priority:

1. Array QC per probe and per sample, the corrections the calls need (GC wave, the cohort reference), and a verdict
   per sample: pass, or flagged with the reason.
2. The sex-chromosome complement of every sample from the X and Y signals - XX, XY, XXY, XYY, X, XXX and the rest,
   constitutional or mosaic (loss of Y in part of the cells, loss or gain of an X in part of the cells) - checked
   against any reported sex.
3. Relatedness between every pair of samples from the genotypes alone, the pedigree it implies (duplicates or
   monozygotic twins, parent-offspring, full siblings, second and third degree; trios and duos assembled), and the
   comparison with any pedigree given.
4. Near-chromosomal events on every reliable sample: whole chromosomes and arms, and segments of 2 Mb or more,
   gained, lost or copy-neutral loss of heterozygosity, in every cell or in part of them, with the cell fraction;
   runs of homozygosity; a karyotype string per sample.
5. The report, and the digest of the largest events for hand-off.
6. The public example validated and documented.

Out of scope for phase 1 (write the hooks, not the features): sub-megabase CNVs (PennCNV territory), haplotype
phasing from a reference panel for very low cell fractions (MoChA territory), imputation, ancestry. Phase 2 will add
family-aware readings where a trio or duo exists: parent of origin of a gain or a uniparental disomy from the
parents' opposite homozygotes, inheritance, Mendelian-error rates, and transmission phasing as triokaryo does.

## 2. The input contract

One VCF or BCF, indexed (.tbi or .csi), all samples of the cohort, as `bcftools +gtc2vcf` writes an Illumina export.
What such an export carries, observed on the public example (section 9) - design for this and tolerate the rest:

- `FORMAT/GT`, `FORMAT/LRR`, `FORMAT/BAF` are required; `arraykaryo inspect` refuses an export without them. Other
  FORMAT fields (GQ, IGC, NORMX, NORMY, R, THETA, X, Y) are ignored.
- GT is written diploid everywhere, on a male's X and Y too (`0/0`, `0/1`, `1/1`); a no-call is `./.`.
- A missing LRR or BAF is `nan`. BAF is not clipped: values slightly below 0 and above 1 occur; clip to [0, 1] before
  any band statistic.
- On the Y, a female's probes are no-calls with LRR around -3 to -5 (the Y intensity is reported against the male
  cluster); a male's Y probes are called, LRR around 0. This makes the Y level a clean sex signal.
- `INFO/INTENSITY_ONLY` flags probes with intensity but no genotype (the record's ALT is `.`): LRR usable, BAF and
  GT not. `INFO/GC` is the GC content around the probe (the GC-wave covariate). `INFO/ALLELE_A`, `ALLELE_B`,
  `GenTrain_Score`, `Cluster_Sep` and the cluster statistics may be present; use them for probe QC when present,
  never require them.
- Probes are SNVs almost everywhere; skip multiallelic and indel records for genotype statistics, keep their LRR.
- Contigs may or may not carry the `chr` prefix; GRCh37, GRCh38 and T2T-CHM13v2.0 are recognised from the contig
  lengths (`genome.py`); refuse an unknown assembly unless `--genome` names a table.
- Scale: the example is 2,141 samples by about 2.4 million probes in a 132-GB BCF. `bcftools query` streams one
  chromosome (chr21, 32,823 records, GT:LRR:BAF for every sample) in 11 s on a laptop, so the whole export is a
  single pass of roughly 15 minutes; the design must not read the export more than once (section 4.1).

Optional inputs, each read by header names, never by column position:

- a sample sheet: `sample`, optionally `reported_sex` (M/F/1/2/male/female), `batch` (plate, chip, run), any other
  columns carried through to the report;
- a pedigree in the house form `#kid dad mom kid_sex dad_sex mom_sex` (sexes 1/2/0), the same file the other
  pipelines take; also accepted: a PLINK-style `.fam`/ped with the standard six columns. A pedigree is compared
  with what the genotypes say and never trusted over them;
- `--samples` to restrict to a subset; `--exclude-probes` and `--exclude-samples` lists.

## 3. Genome assemblies

`src/arraykaryo/genome.py` recognises GRCh37, GRCh38 and CHM13 by chr1/chrX/chrY lengths and holds the
pseudoautosomal regions. VERIFY every coordinate there against the assembly's own annotation (UCSC for hg19/hg38,
the T2T consortium for CHM13v2.0) and empirically: a male's heterozygous X calls are frequent inside a PAR and near
absent beyond it (on the example, 145 of 610 calls in PAR1 against 0 of 567 on distal Xq). Vendor per assembly, with
a `data/` file each and its provenance in `docs/METHODS.md`: centromere positions (arm boundaries), cytogenetic
bands, the acrocentric short arms (no probes: an arm with no probes is not a loss), and the PARs. Allow a custom
table by `--genome-table` for any other assembly.

## 4. Architecture and commands

```
arraykaryo inspect   cohort.bcf                       what the export carries (done)
arraykaryo mock      --out DIR [...]                  a simulated export with families, planted events and QC failures, with truth.json
arraykaryo extract   --bcf cohort.bcf --out CACHE     one pass into per-chromosome arrays (GT int8, LRR, BAF), a probe table, a manifest
arraykaryo qc        --cache CACHE [--sample-sheet]   probe QC, sample QC, GC-wave correction, the cohort reference; qc/ tables + figures
arraykaryo sex       --cache CACHE                    the complement per sample, mosaics included, against the reported sex
arraykaryo pedigree  --cache CACHE [--pedigree]       relatedness, families, the check against a given pedigree
arraykaryo events    --cache CACHE [--pedigree]       near-chromosomal events and ROH per sample; karyotype strings; figures
arraykaryo report    --cache CACHE                    the cohort page, per-sample pages for flagged samples, the guide, the digest
arraykaryo run       --bcf cohort.bcf --out DIR       all of the above
```

Every command reads and writes under one output directory (`--out`), is idempotent (skips what is done unless
`--redo`), and logs what it did and how long it took. `--threads N` parallelises over chromosomes (extract) or
samples (events, figures).

### 4.1 extract

One pass over the export into a cache the rest of the pipeline reads: per chromosome, memory-mapped arrays of
shape probes x samples - `gt` int8 (0 hom-ref, 1 het, 2 hom-alt, -1 no-call), `lrr` float32, `baf` float16 or
float32 (decide by a measured precision check: BAF noise at a heterozygous site is about 0.03, so float16's 0.001
resolution costs nothing; keep LRR float32) - plus `probes.tsv` (chrom, pos, id, ref, alt, intensity_only, gc, and
the cluster fields when present), `samples.txt` in the export's order, and `manifest.json` (export path as given,
size, modification time, md5 of the header, sample count, per-chromosome record counts, the arraykaryo version).
Resumable per chromosome; `--threads` runs chromosomes in parallel.

Two routes to time on the mock scaled to 2,000 samples x 200,000 probes (one chromosome) and keep the faster:
(a) `bcftools query -r CHROM -f '%POS\t%ID\t%REF\t%ALT\t%INFO/INTENSITY_ONLY\t%INFO/GC[\t%GT\t%LRR\t%BAF]\n'` piped
into a vectorised parser (numpy on the bytes; GT is one of four 3-byte tokens, the floats parse in bulk); (b) cyvcf2
with its per-record numpy arrays (`gt_types`, `format('LRR')`, `format('BAF')`). pysam's per-sample record access is
too slow at this scale and is for `inspect` only. Target: the full example (2,141 x 2.4M) in 20 minutes single
stream, under 6 minutes with 8 threads; the cache about 26 GB with float16 BAF; peak memory under 4 GB.

### 4.2 qc

**Probe QC** (cohort-wide, written to `qc/probes.tsv` with a `keep` column and the reason when not kept):
intensity-only probes keep their LRR and lose BAF/GT; duplicate positions; non-SNV; genotype call rate across
samples below 0.95; per-probe LRR median beyond +-0.3 or LRR robust SD in the top 1% (mis-normalised or noisy);
a heterozygous-call BAF whose cohort median departs from 0.5 by more than 0.05 (a shifted cluster); Hardy-Weinberg
failure at p < 1e-10 among unrelated passing samples (after the pedigree step; two passes are acceptable); and,
when present, `GenTrain_Score` below 0.6 or `Cluster_Sep` below 0.3. PAR probes are kept but tagged `par`, and are
used for the autosomal-like X only. Document each threshold with its source (PennCNV's and MoChA's practice,
Illumina's GenTrain guidance) and make each a parameter.

**Sample QC** (`qc/samples.tsv`; each metric with the cohort median and robust SD, a `pass` column and `flags`):
call rate (flag below 0.97); autosomal LRR SD after GC correction (flag above 0.20 for mosaic calling, a hard
fail above 0.30 - PennCNV's and MoChA's values; VERIFY against the cohort's distribution and state the platform);
BAF SD at heterozygous calls (flag above 0.06); BAF drift, PennCNV's share of homozygous-called probes whose BAF
lies in (0.2, 0.25) or (0.75, 0.8) (flag above 0.002); the wave factor, the slope of LRR on GC before correction
(flag |WF| above 0.05); LRR mean (flag beyond +-0.1); the sex check result (section 4.3); a duplicate or an
unexpected relative (section 4.4); and, after the events step, an event count beyond 3 robust SD above the median.
A sample with a hard fail is "not to be read"; a flagged sample is read with the flag shown beside every call.

**Corrections**, in this order, per sample, on the autosomes and the X: (1) the GC wave - regress LRR on the
probe's `INFO/GC` and on the GC of its 1-Mb window with a low-order polynomial (Diskin 2008's model; MoChA's
`--adjust-BAF-LRR` is the modern reference) and subtract the fit; (2) the cohort reference - subtract, per probe,
the median corrected LRR over the passing samples, per batch when a batch column gives 20 or more passing samples
per batch, else cohort-wide (this removes the shared array structure exactly as triokaryo's panel does for
sequencing); (3) BAF: clip to [0, 1]; re-centre heterozygous-call BAF per probe by the cohort's median at that
probe (MoChA's cluster adjustment) when the shift exceeds 0.01. Write the before/after SD per sample so the
correction's gain is visible, and a figure of LRR against GC before and after for a typical sample.

### 4.3 sex

Per sample, from the X outside the PARs and the Y: the heterozygous rate on the X relative to the sample's
autosomal rate; the median corrected LRR on the X; the median LRR on the Y; the Y call rate. The cohort's two X
modes (one copy, two copies) and two Y modes (present, absent) are found by mixture fitting on the passing samples
and give the copy scale - do not assume the sequencing values (-1 for one copy): arrays compress LRR, a one-copy X
sits near -0.5, and the scale is the cohort's own. Copies of X = 1 + (LRR_X - m1) / (m2 - m1) on that scale; a Y in
part of the cells from the Y level between the male and female modes (Forsberg 2014's mLRR-Y reading, calibrated
on the cohort); an X lost or gained in part of the cells of an XX sample from the BAF split at X heterozygous calls
(the same cell-fraction formulas as section 4.5) and the X level together; an extra X in part of the cells of an
XY sample from heterozygous-looking X calls appearing (their BAF sits away from 0.5 by (1 - f) / (2 (1 + f))) and
the X level rising. Name the complement (XX, XY, XXY, XYY, X, XXX, XXYY, ...) and every mosaic with its fraction,
exactly as triokaryo's sexchrom does, and compare with the reported sex when given: agree; a Y in a reported female
(XY female, or a swap); no Y in a reported male (XX male, or a swap); 47,XXY and 45,X named for what they are; a
sex inferred elsewhere from X heterozygosity alone (peddy's) reads 47,XXY as female and 45,X as male - say so in
the check's wording, as triokaryo does. Outputs: `sex/sex.tsv`, a scatter of X level against Y level coloured by
X heterozygosity with the complements labelled in the sidecar, and per-sample X and Y tracks for every sample not
plainly XX or XY.

### 4.4 pedigree

Relatedness from the genotypes alone, on autosomal SNV probes passing QC with cohort call rate 0.98 or more and
minor allele frequency 0.05 or more, thinned to one probe per 50 kb (about 50,000 probes; a distance thinning is
an adequate stand-in for LD pruning at this density - VERIFY on the example that the kinship estimates of the
known trios do not move when the thinning is halved). Estimator: KING-robust kinship (Manichaikul 2010), the
between-family estimator as PLINK 2 `--make-king` implements it; implement from the paper's equation and test
against PLINK 2 on the mock when PLINK 2 is on PATH (optional test), and on the known 1000 Genomes families.
Compute every pair at once with matrix products of the 0/1 indicator matrices (heterozygous, homozygous-reference,
homozygous-alternate; missing masked): 2,141 samples by 50,000 probes is seconds of BLAS. Alongside kinship,
IBS0 (the share of probes where the pair is opposite homozygous) separates parent-offspring (IBS0 near zero, the
genotyping-error floor learned from the cohort) from full siblings (IBS0 clearly above it).

Classes by kinship: above 0.354 duplicate or monozygotic; 0.177 to 0.354 first degree, split by IBS0 into
parent-offspring and full sibling; 0.0884 to 0.177 second degree; 0.0442 to 0.0884 third degree; the boundaries
are KING's and each is a parameter. Assemble families from the first-degree graph: a trio is a sample with two
parent-offspring relatives who are unrelated to each other, of opposite sex by section 4.3, and whose three
genotypes show a Mendelian-error rate at the error floor in that orientation and far above it in any other - the
orientation test resolves who is the child without ages. A parent-offspring pair with no third member is a duo,
oriented by the given pedigree or a sample sheet age when available, else reported as an unoriented pair.
Outputs: `pedigree/relatedness.tsv` (every pair of third degree or closer: kinship, IBS0, class), `pedigree/
pedigree.inferred.tsv` in the house form (with sexes from section 4.3), `pedigree/families.tsv`, and, when a
pedigree was given, `pedigree/check.tsv` with one row per claimed relationship and per discovered one: agrees;
a claimed parent unrelated (non-paternity or a sample swap - the sex check and the duplicate table say which);
an undisclosed first- or second-degree pair; a duplicate; a sex contradiction. A kinship-against-IBS0 scatter
with the classes in the sidecar.

### 4.5 events

Every passing sample, every chromosome, from the corrected LRR and the BAF, with the genotypes for homozygosity:

1. **Tracks per window** of about 100 probes (pool by probe count so that a window is a similar amount of
   evidence everywhere): median LRR; the BAF deviation d at heterozygous calls, estimated per window as the
   location of the symmetric two-band mixture 0.5 +- d (fit by maximum likelihood with a shared band SD; the
   per-window mean of |BAF - 0.5| is a biased stand-in that must not be used for the cell fraction); the observed
   heterozygous rate against the rate expected from the cohort's allele frequencies at those probes (Hardy-Weinberg
   expectation; the ratio is the homozygosity signal).
2. **Segmentation** of the three tracks jointly, chromosome by chromosome: circular binary segmentation on the
   window means (implement it in numpy with a permutation significance, or measure `ruptures`' PELT and keep it if
   it is faster and as accurate on the mock), then refine every boundary to probe resolution. Minimum segment
   2 Mb (`--min-mb`), minimum 40 probes.
3. **Typing** by the LRR step and the BAF split: gain (LRR up, d > 0), loss (LRR down, d > 0), copy-neutral loss of
   heterozygosity (LRR flat, d > 0 at heterozygous calls - a mosaic), and, where the heterozygous rate falls to the
   error floor with no split, a run of homozygosity (LRR flat: autozygosity or a constitutional isodisomy; LRR down:
   a constitutional deletion, one copy). A segment covering 90% or more of a chromosome's probed span is `whole`,
   90% or more of an arm `p` or `q`, else `stretch`.
4. **Cell fraction** from the BAF split, the primary estimate: gain f = 4d / (1 - 2d), loss f = 4d / (1 + 2d),
   copy-neutral LOH f = 2d (these are the band positions (1 + f) / (2 + f), 1 / (2 - f) and (1 + f) / 2 solved for
   f); and from LRR on the cohort's copy scale, the secondary estimate: calibrate the one-copy level from the male
   X (section 4.3) and the three-copy level from any constitutional trisomy found, else from the array literature
   (PennCNV's defaults, VERIFY), and state the calibration in the report. Report both and the disagreement.
   A whole chromosome or arm with f at or above 0.9 is constitutional.
5. **Runs of homozygosity**: from the genotypes with a sliding window (PLINK's `--homozyg` defaults as the
   starting point: 1 Mb, 50 probes, a few heterozygous calls allowed as errors) or an HMM on the per-probe
   homozygosity with the cohort's allele frequencies (bcftools roh's model); keep runs of 1.5 Mb or more; report
   per sample the genome fraction in ROH (F_ROH) and every run of 10 Mb or more as an event; a whole-chromosome ROH
   is an isodisomic uniparental disomy or consanguinity - say both, and let the trio addition (phase 2) decide.
6. **The sex chromosomes** enter as whole-chromosome events from section 4.3 (an X gained or lost, in every cell
   or in part; a Y lost or gained), with the same cell fractions.
7. **Floors and flags**: the target sensitivity is a cell fraction of 0.10 for a whole chromosome or arm and 0.20
   for a 2-10 Mb stretch at an LRR SD of 0.15 (VERIFY on the mock and record the detection curve in
   `docs/calibration/`); a null sample must receive, on the mock, no event in 19 of 20 samples. Every event carries
   its evidence (LRR step, d, probes, heterozygous calls, the two cell fractions) and a note where the evidence is
   thin or the sample is flagged. A known-artefact list per assembly (the MHC, the 8p23 and 17q21.31 inversions,
   pericentromeric and subtelomeric probe-sparse regions) tags events there; it does not remove them.
8. **Karyotype string** per sample in the ISCN-like grammar triokaryo uses (`47,XXY`; `mos 47,XY,+12[0.35]/46,XY`;
   `46,XX,del(2)(q22.1q31.1)[0.67]`; `46,XX,upd(X)(iso)?` for a whole-chromosome ROH; `46,XY,roh(6)(p22.1p21.1)`),
   with cytogenetic bands from the vendored tables.

Outputs: `events/events.all.tsv` (one row per event, every number above), `events/roh.all.tsv`,
`events/karyotypes.tsv`, and figures for every sample with an event or a flag: the genome-wide LRR/BAF/homozygosity
figure (every chromosome in one image, the karyogram-like view people expect) and one figure per chromosome with
an event. Figures are drawn only on demand for the rest (`--figures all`).

### 4.6 report

`report/index.html`, self-contained: tiles (samples; passing; flagged by reason; families found; samples with
events; events by type and class); the QC section with the distributions (LRR SD, BAF SD, call rate, wave factor,
drift) as histograms with the thresholds marked in the sidecar, and the sortable sample table; the sex section
(the X-against-Y figure, the table of every sample not plainly XX or XY, the comparison with the reported sex); the
pedigree section (families, the check against the given pedigree, the kinship-against-IBS0 figure); the events
section (the landscape - events per chromosome, one row per sample with an event - and the sortable table with
links to the per-sample pages); the ROH section (F_ROH distribution; the long runs); a "trios set aside" section
like triokaryo's cohort guards (relationship or noise outliers); and the guide page explaining every column and
colour, with pattern cards of what each event kind looks like on an array (LRR and BAF bands for a constitutional
trisomy, a mosaic gain, a mosaic loss, mosaic copy-neutral LOH, a ROH, 47,XXY, mosaic loss of Y). Per-sample pages
(`report/samples/<sample>.html`) for flagged samples and samples with events, each with its figures, its QC row,
its complement and its events in plain words. A digest (`report/digest.html`) of the largest events across the
cohort for hand-off, ranked as triokaryo's digest ranks (whole chromosome, then arm, then stretch by size; every
cell before part of the cells), each with its figures embedded.

### 4.7 mock

`arraykaryo mock --out DIR --samples N --seed S [...]` writes an export in the gtc2vcf style (`.vcf.gz` and `.bcf`,
indexed; diploid GT, `nan` for missing, BAF unclipped, intensity-only probes with ALT `.`, `INFO/GC`, PAR probes,
the assembly's contig lengths - CHM13 by default, `--genome grch38` as the alternative), a sample sheet with
reported sexes and two batches, a pedigree, and `truth.json`. Probe density about one per 1.2 kb on a reduced
contig set for speed (`--contigs chr1,chr7,chr12,chr15,chr21,chrX,chrY` default; `--full` for all), with
realistic LRR compression (one copy near -0.55, three copies near +0.35 on the LRR scale), per-sample noise
(LRR SD drawn between 0.10 and 0.25, BAF SD at heterozygous calls between 0.03 and 0.06), a per-sample GC wave with
amplitude drawn around 0.03, a per-batch per-probe offset, missingness, and families (`--trios K --duos D`) with
Mendelian transmission. Planted by default, each in its own sample and recorded in truth.json with sample,
chromosome, span, type and cell fraction: a constitutional trisomy 21; a mosaic gain of chromosome 12 at 0.30; a
mosaic loss of 7q at 0.25; a mosaic copy-neutral LOH of 11p at 0.40; a constitutional 10-Mb deletion of 1p36; a
15-Mb ROH; a whole-chromosome isodisomy of 15 in a trio child; 47,XXY; 45,X; mosaic loss of Y at 0.5 in a man;
mosaic loss of an X at 0.3 in a woman; a sample swap between a man and a woman (so a sex mismatch and a broken
trio); a non-paternity; a duplicate pair; two samples with LRR SD above the hard fail; and null samples for the
false-positive count. Every default can be switched off or scaled (`--no-events`, `--noise`, `--cell-fraction`).

### 4.8 run

The chain with one command and one log; refuses to proceed past a step that failed; prints where every table and
page is at the end.

## 5. Performance

Targets on a laptop (8 cores, 16 GB) and the example cohort: extract under 6 minutes with 8 threads; qc under
3 minutes; sex under 1 minute; pedigree under 3 minutes; events under 10 minutes for 2,141 samples with 8 threads
(vectorised windows per sample; samples in parallel; no Python loop over probes); report under 3 minutes with
figures for flagged samples only. Peak memory under 4 GB (memory-mapped cache; one chromosome of one sample block
in memory at a time). The mock cohort of 60 samples runs end to end, with figures, in under 3 minutes in CI. Keep a
`benchmarks/` note with the measured times per version.

## 6. Outputs and conventions

Every table is TSV with a header; booleans 1/0; missing `NA`; positions 1-based inclusive; sizes in Mb to one
decimal in the pages, in bp in the tables; sample identifiers exactly as the export names them. Every figure has
its sidecar and legend image. `docs/OUTPUT.md` lists every file and column; `docs/METHODS.md` the methods with the
references (Peiffer 2006 for LRR and BAF; Wang 2007 PennCNV; Colella 2007 QuantiSNP; Diskin 2008 for the GC wave;
Manichaikul 2010 KING; Jacobs 2012 and Laurie 2012 for mosaic detection from arrays; Forsberg 2014 for mosaic loss
of Y; Machiela 2016 for mosaic loss of X; Loh 2018 MoChA; Purcell 2007 PLINK for ROH; Narasimhan 2016 bcftools
roh; Pedersen 2017 peddy for sex and relatedness checks) and every threshold with its source.

## 7. Tests and CI

pytest, on the mock only: extraction round-trips the mock exactly (GT, LRR, BAF per sample and probe); the probe
and sample QC flag the planted failures and nothing else; the GC-wave correction removes the planted wave (the
LRR-GC correlation falls under 0.05); the sex calls recover every planted complement and mosaic within 0.1 of the
fraction, and name the swap; the pedigree recovers every planted family and flags the non-paternity, the swap and
the duplicate, with no spurious first-degree pair; every planted event is found with its type, class and cell
fraction within 0.08, and the ROH; null samples receive no event in 19 of 20; the karyotype strings match the
truth; the report and digest contain what the tables contain; the CLI runs end to end; the image runs the mock in
CI. Keep tests under about 10 minutes in total; mark the slow ones.

CI (`.github/workflows/ci.yml`, in place): tests on Python 3.9 and 3.12 with bcftools installed; the image built,
run on the mock, and pushed to GHCR with the tags `main` and `sha-<commit>` on main.

## 8. Milestones

- M1: `mock` and `extract`, with the timing study of the two extraction routes. Version 0.1.0.
- M2: `qc` and `sex`. 0.2.0.
- M3: `pedigree`. 0.3.0.
- M4: `events` and ROH, karyotype strings, figures, the calibration curve. 0.4.0.
- M5: `report`, per-sample pages, the guide, the digest, `run`. 0.5.0.
- M6: the public example run, validated and documented in `docs/example/` (section 9). 1.0.0.

## 9. Validation on the public example

A public 1000 Genomes array export (2,141 samples, about 2.4 million probes, CHM13 coordinates in the copy at
hand; other cohorts are GRCh38 or GRCh37) is pointed at by `ARRAY_BCF`, and the IGSR 3,202-sample pedigree table
by `KNOWN_PED`; neither lives in the repository, and no location is written anywhere in it. Check, and write the
counts into `docs/example/README.md`:

- the inferred pedigree recovers the known trios and duos present (for instance NA12878 with NA12891 and NA12892;
  NA19240 with NA19238 and NA19239; HG00733 with HG00731 and HG00732; HG00702 with HG00703 and HG00704; NA19675
  with NA19678 and NA19679), and every first-degree pair found is in the known pedigree or explained;
- the sex-chromosome complements agree with the known sexes except where an aneuploidy or a mosaic is named with
  its evidence (lymphoblastoid lines lose the X and the Y in part of their cells: expect some);
- the samples triokaryo read from whole-genome sequencing agree: HG01103 carries mosaic losses on 2q (about 0.67)
  and 14q (about 0.64), NA12748 a mosaic gain of 13q (about 0.69) and a mosaic terminal loss of 21q (about 0.59);
  both, with their relatives, are on the array;
- culture artefacts of lymphoblastoid lines (gains of chromosome 12, X losses in women, other mosaics) are listed
  as findings with their cell fraction, not treated as errors;
- the QC distributions and the thresholds make sense for this platform, written down.

## 10. Assumptions made here, to be confirmed by the owner

- Pure Python plus bcftools; MoChA, PennCNV and PLINK are references and optional cross-checks, never dependencies.
- Reported sex comes from a sample sheet column or the pedigree's sex columns; absent both, the genotypes alone.
- The event floor is 2 Mb; the digest lists whole chromosomes, arms and stretches of 10 Mb or more by default.
- Per-sample pages are written for flagged samples and samples with events only, so that a cohort of thousands
  yields a report of manageable size; `--pages all` for every sample.
- Families improve calls (parent of origin, inheritance, Mendelian errors) in phase 2; phase 1 calls every sample
  alone and reports families as a finding of the pedigree step.

## 11. Definition of done, per milestone

Tests green on both Python versions; the image builds and runs the mock; `docs/METHODS.md` and `docs/OUTPUT.md`
describe what the milestone added, with references and thresholds; the version bumped; the README's status line
updated; no path, no real data, no identifier of a non-public sample anywhere in the repository; and a commit
message that a colleague can read in a minute.
