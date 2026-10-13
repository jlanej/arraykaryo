# Methods

What arraykaryo computes, how, and where every threshold comes from. Sections follow the pipeline. Numbers marked
*measured* were measured on the mock at the stated scale (`arraykaryo mock`), never on real data in this repository;
numbers marked *example* come from the public 1000 Genomes array export described in `docs/example/README.md`, which lives
outside the repository.

## Input

One VCF or BCF, indexed, as `bcftools +gtc2vcf` writes an Illumina export: `FORMAT/GT`, `FORMAT/LRR` and `FORMAT/BAF` per
call (the log R ratio and B-allele frequency of Peiffer 2006). Everything else is ignored. GT is diploid everywhere, on a
man's X and Y too; a no-call is `./.`; a missing LRR or BAF is `nan`; BAF is not clipped (values slightly outside [0, 1]
occur) and is clipped to [0, 1] before any band statistic. `INFO/INTENSITY_ONLY` (ALT `.`) marks probes with intensity but
no genotype: their LRR is used, their BAF and GT are not. Multiallelic and indel records keep their LRR and lose their
genotype. `arraykaryo inspect` reports what an export carries and refuses one without the three FORMAT fields.

## Genome tables (`genome.py`, `data/`)

The assembly is recognised from the header's chr1, chrX and chrY lengths: GRCh37 (hg19/b37, with or without the `chr`
prefix), GRCh38 (hg38) or T2T-CHM13v2.0; `--genome` names one, or a band table in the same format, for anything else.
Per assembly the package ships the cytogenetic bands, from which the chromosome lengths, the centromere (the p/q boundary of
the two `acen` bands; the end of the short arm) and the acrocentric short arms (13, 14, 15, 21, 22: ribosomal DNA, no probes -
an arm with no probes is not a loss) follow:

| assembly | source | file |
| --- | --- | --- |
| GRCh37 | UCSC hg19 `cytoBand.txt.gz` (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/) | `data/cytoband.grch37.tsv` |
| GRCh38 | UCSC hg38 `cytoBand.txt.gz` (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/) | `data/cytoband.grch38.tsv` |
| CHM13 | T2T consortium `chm13v2.0_cytobands_allchrs.bed` (https://github.com/marbl/CHM13, annotation) | `data/cytoband.chm13.tsv` |

`scripts/vendor_genome_tables.py` rewrites them from the downloaded files and records a checksum in each file's header; it
checks that every chromosome is tiled from 0 to its length without a gap and that each carries a p and a q `acen` band.

The pseudoautosomal regions (1-based, inclusive; X and Y) are typed into `genome.py` from UCSC's hg19 and hg38 `par`
tables and the T2T consortium's `chm13v2.0_PAR.bed` (analysis set):

| assembly | PAR1 | PAR2 |
| --- | --- | --- |
| GRCh37 | chrX:60,001-2,699,520; chrY:10,001-2,649,520 | chrX:154,931,044-155,260,560; chrY:59,034,050-59,363,566 |
| GRCh38 | chrX:10,001-2,781,479; chrY:10,001-2,781,479 | chrX:155,701,383-156,030,895; chrY:56,887,903-57,217,415 |
| CHM13 | chrX:1-2,394,410; chrY:1-2,458,320 | chrX:153,925,835-154,259,566; chrY:62,122,810-62,460,029 |

Empirical check (*example*, three men, one record in twenty): heterozygous X calls in PAR1 26 of 96 calls, on distal Xq
(150 Mb to the end) 1 of 170. The export carries PAR probes on the Y as well, with diploid calls in both sexes; they are
tagged `par` and never read as the Y. Beyond the PAR, the Y of the example shows, at 3.0-4.8 Mb, probes that are
heterozygous in men and called at LRR 0 in women: the X-transposed region. The Y level is therefore read over the probes
that are empirically Y-specific (section "sex"), not over a table.

Regions where array signal is known to mislead are tagged on every event, never removed: the MHC (6p22.1-p21.32), the 8p23.1
and 17q21.31 inversion polymorphisms (by band, so that they resolve on any assembly with a band table), the centromere
with 2 Mb each side, the first and last 1 Mb of every chromosome, and the acrocentric short arms.

## Extraction (`extract.py`, `cache.py`)

One pass over the export into a cache the rest of the pipeline reads: per chromosome, `.npy` arrays of shape samples x
probes - `gt` int8 (0 hom-ref, 1 het, 2 hom-alt, -1 no-call), `lrr` float32, `baf` float16 - plus `probes.tsv` (and the
binary `probes.npz`), `samples.txt` and `manifest.json`. Sample-major layout: one sample's chromosome is one contiguous
read (the per-sample figures and calls); a block of probes over all samples is a strided read of contiguous runs (QC, the
cohort reference). Records are written in blocks of 8,192 with positioned writes, so the process holds one block (peak
RSS 0.49 GB *measured* at 2,000 samples; a memory-mapped writer reached 3-4 GB of dirty pages). Resumable per chromosome
(`done/<chrom>.json`), chromosomes in parallel with `--threads`.

**The two routes, timed** (*measured*: the mock scaled to 2,000 samples x 196,008 probes on one chromosome, one thread,
a laptop; `--route` selects either):

| route | wall time | peak RSS | note |
| --- | --- | --- | --- |
| (a) `bcftools query` streamed into a vectorised numpy parser | 164 s | 4.0 GB (memory-mapped writer) | `bcftools query` alone takes 41 s to print the 3.9 GB of text; the float parse dominates |
| (b) cyvcf2 (`gt_types`, `format('LRR')`, `format('BAF')` filled in C) | 33 s (41 s with positioned writes) | 0.49 GB | **kept as the default** |

cyvcf2 earns its place as a dependency (five times faster); route (a) remains the fallback when it is not installed (pure
`bcftools` + numpy). Both produce bit-identical caches (tested). Projected for the example (2,141 x 2.4 million):
about 7 minutes in one stream, under 2 minutes with 8 threads.

**BAF as float16** (*measured*): the largest rounding difference against the float32 values is 0.00049 (mean 0.00009),
against a band SD of 0.03-0.05 at a heterozygous call; the band deviation d fitted on 200 sets of heterozygous calls
(d = 0, 0.03, 0.07, 0.17) differs by 0.0000 between the two storages. LRR stays float32. The cache is 7 bytes per cell:
2.6 GB for the 2,000 x 196,008 mock; about 36 GB for the example.

## QC (`qc.py`)

**Probe QC**, cohort-wide (`qc/probes.tsv`: `keep` = the LRR is usable, `keep_gt` = the BAF and genotype are usable,
`reason`):

| rule | threshold | source |
| --- | --- | --- |
| intensity-only probe (`INFO/INTENSITY_ONLY`, ALT `.`) | LRR kept, BAF and GT not | the export's definition |
| non-SNV record (indel, multiallelic) | genotype not used | the brief: genotype statistics on SNV probes |
| duplicate position | second record dropped | - |
| genotype call rate across samples | < 0.95 | PennCNV/MoChA practice (MoChA plots variants at 0.97) |
| per-probe LRR median (raw) | beyond +-0.3 | a mis-normalised probe |
| per-probe LRR robust SD | top 1% | a noisy probe |
| heterozygous-call BAF median | more than 0.05 from 0.5 | a shifted cluster; MoChA re-centres shifts above 0.01 |
| Hardy-Weinberg exact test (Wigginton 2005), with an excess of heterozygotes | p < 1e-10 | PLINK `--hwe` practice for arrays; over the passing samples here, over the unrelated passing samples after the pedigree step. Only an excess of heterozygotes marks a bad probe (one reading two loci); a deficit is what population structure does to a cohort of several populations (the Wahlund effect: *example*, 221,583 of 2.4 million probes fail at 1e-10 in both directions, most by deficit) and is left alone |
| `GenTrain_Score` | < 0.6 | Illumina's GenTrain guidance |
| `Cluster_Sep` | < 0.3 | Illumina: overlapping clusters |
| PAR probe | kept, tagged `par` | used for the autosomal-like X only |
| acrocentric short arm (13p, 14p, 15p, 21p, 22p) | not used | an array's probes lifted into an assembly that carries those arms (CHM13 does) are mismapped: *example*, 1,226 spurious arm events before the rule |

On the Y the call rate is taken among the samples with a Y (its PAR probes among everyone). Every threshold is a parameter
(`THRESHOLDS` in `qc.py`, written with its source to `qc/summary.json`).

**Sample QC** (`qc/samples.tsv`; each metric's cohort median and robust SD in `qc/summary.json`):

| metric | definition | flag | source |
| --- | --- | --- | --- |
| call rate | called / autosomal SNV probes kept for genotypes | < 0.97 | Illumina's and MoChA's recommendation |
| LRR SD | SD of the corrected autosomal LRR | > 0.20 flag; > 0.30 hard fail | MoChA's mosaic practice; PennCNV `filter_cnv --qclrrsd 0.3` |
| BAF SD | SD of BAF - 0.5 at heterozygous calls, the median over the autosomes of the per-chromosome value (so one chromosome's split bands do not flag a sample) | > 0.06 | the band SD that blurs the bands (MoChA's `baf_auto` flags contamination at 0.03) |
| BAF SD, robust | the median of \|BAF - 0.5\| at heterozygous calls over 0.6745 (a half-normal's median): the band SD the fits use, immune to mis-called genotypes, which inflate the plain SD (0.056 for a planted 0.040 on the mock) | - | - |
| BAF drift | share of homozygous-called probes with BAF in (0.2, 0.25) or (0.75, 0.8) | > 0.002 | PennCNV's BAF_DRIFT (over every probe, filtered at 0.01) |
| wave factor | PennCNV's waviness factor: the median absolute deviation of the 1-Mb window medians of the raw autosomal LRR from their overall median, signed by the correlation of the window medians with the window GC (PennCNV signs it against a reference wave) | abs > 0.05 | PennCNV `filter_cnv --qcwf 0.05` (`detect_cnv` warns at 0.04) |
| LRR mean | mean raw autosomal LRR | beyond +-0.1 | a poorly normalised sample |
| sex check | section "sex" | `sex_mismatch` when an XX/XY complement contradicts the reported sex | - |
| relatedness | section "pedigree" | a duplicate or an undisclosed first-degree relative | - |
| event count | section "events" | > 3 robust SD above the cohort median | - |

A hard failure marks the sample "not to be read"; a flag is shown beside every call. Flags, not drops: every sample stays
in every table; the cohort reference and the pedigree are built from the passing samples only. *Measured* on the default
mock (60 samples): the six planted failures (two LRR SD above the hard fail, two strong waves, one high BAF SD, one low
call rate) are flagged with their reason and nothing else is.

**Corrections**, in this order, per sample, on the autosomes and the X:

1. *The GC wave* (Diskin 2008; MoChA's `--adjust-BAF-LRR` is the modern reference): LRR regressed on the probe's
   `INFO/GC` and on the mean GC of the probes within 500 kb (the 1-Mb window), each standardised and entered as a
   quadratic, fitted by least squares over the autosomes on the LRR clipped to +-1 (robust to outliers) and on
   *within-arm deviations* - every chromosome arm's own mean of the LRR and of the covariates is removed before the
   products are formed - so that a whole chromosome or arm gained or lost carries no leverage on the slopes (fitted on all
   probes, the regression absorbed a third of a planted arm gain whose arm differed in GC from the rest). The fit is
   subtracted and the corrected LRR centred per sample on its median over the autosomes (`lrr_offset`). *Measured* on the
   mock: the correlation of LRR with the window GC falls from a median of 0.04-0.09 before to 0.000 after on the null
   samples; the autosomal LRR SD from 0.175 to 0.150, the planted per-sample noise.
2. *The cohort reference*: the median corrected LRR per probe over the passing samples, per batch when every batch holds
   20 or more passing samples, else cohort-wide, subtracted. This removes the array's shared structure (the per-probe offsets
   every sample carries) exactly as triokaryo's panel does for sequencing. On the X the reference is built from the samples
   without a Y, so that one X reads as one copy on the cohort's scale; the Y keeps its GC-corrected raw level against the
   male cluster (a reference over a mixed cohort would sit between the sexes).
3. *BAF*: clipped to [0, 1] at use; the heterozygous-call BAF is re-centred per probe by the cohort's median at that probe
   when its shift exceeds 0.01 and 30 or more heterozygous calls support the median (MoChA's cluster adjustment; the median
   of fewer calls is noise and was adding it), and a probe whose shift exceeds 0.05 loses its genotype.

The corrected LRR is never written: `qc/correction.npz` holds the parameters (five coefficients per sample, the reference
per probe and per batch) and `Correction.apply` reproduces it from the cache for any block.

## Sex (`sex.py`)

Per sample: the heterozygous rate on the X outside the PARs relative to the sample's autosomal rate; the median corrected
LRR on the X outside the PARs (reference from the samples without a Y); the median GC-corrected LRR on the Y over the
Y-specific probes; the Y call rate. *Y-specific probes*: the cohort is split by its raw Y level (the midpoint between the
two modes), and a non-PAR Y probe is Y-specific when its median LRR is above -0.6 in the samples with a Y and below -2 in
the samples without, at least 2 apart. This excludes the X-transposed region and other X-homologous probes by what they do,
not by a table (*example*: see "Genome tables"); with fewer than five samples in one group every non-PAR probe is used and the
report says so.

The cohort's two X modes (one copy, two copies) are found by a two-component Gaussian mixture (EM) on the X levels of the
passing samples, initialised from the Y split, and the two Y modes (present, absent) by the medians of the two groups. The
copy scale is the cohort's own: copies of X = 1 + (LRR_X - m1) / (m2 - m1); arrays compress LRR, a one-copy X sits near
-0.5 (*measured* on the mock, which simulates -0.55: m2 - m1 = 0.54-0.57). Copies of Y = 2^((LRR_Y - mY) / (m2 - m1)),
the Y level read on the same compression - Forsberg 2014's mLRR-Y reading, calibrated on the cohort (a Y at half the male
level reads as 0.5 copies; the absent mode reads as 0). With no sample of one sex, the missing mode is set 0.55 from the
other (PennCNV's one-copy level is -0.66) and the summary says so.

The complement the Y implies is the anchor, as in triokaryo's sexchrom: XY with a Y (0.5 copies or more; between 0.25 and
0.5 the nearer complement), XX without - except that one X without heterozygosity (the heterozygous rate under 0.15 of the
autosomal rate) is read against XY whatever the Y: a man whose cells lost the Y and a 45,X are both one X, and the Y's share
of cells says which (*example*: two men whose lymphoblastoid lines lost the Y in 85% and 78% of cells read as 45,X/46,XY
mosaics, not as women with an X loss). **The mosaic floors are the cohort's own**: a deviation in copies must exceed three
robust SD of the levels of the sample's cluster (the X and the Y of the XY samples, the X of the XX samples), never under
0.10 - on the example the X and Y levels of men spread by 0.072 and 0.075 copies (robust SD) against 0.023 for the X of
women, so the floors are 0.21, 0.22 and 0.10; a fixed 0.10 read that spread as mosaic Y gains in 119 men. An X mosaic in an
XX sample needs the BAF split too (d at least 0.03, the null's 95th percentile on the example, with a log-likelihood ratio of
30 or more) unless the level change is large (0.25 or 5 robust SD); an extra X in an XY sample needs heterozygous-looking
calls (the X heterozygous rate above 0.15 of the autosomal rate), and two fully heterozygous X's with a Y are a 47,XXY in
every cell; a single X lost in part of the cells is not a viable state and is noted, not called. So called: an X deviating
from the anchor's count is a mosaic (`mos 45,X[0.30]/46,XX`, `mos 47,XXY[0.20]/46,XY`), at 0.90 or more a constitutional
aneuploidy (47,XXY, 45,X, 47,XXX, 48,XXYY ...); a Y under one copy in a man is a mosaic loss of Y (`mos 45,X[0.50]/46,XY`; a
Y gone from every cell makes the complement X), above one a gain (47,XYY) - tentative when read from the level alone, since
the Y level varies between Y haplogroups (probes deleted or duplicated on some); a Y at 0.25 copies or more in an XX genome
is "a Y in part of the cells (46,XX/47,XXY mosaic, or an XX/XY mixture)". The cell fraction read is the BAF's where it is
significant, else the level's; the two are reported with a note when they disagree by more than 0.2 (on real X losses the
BAF's runs above the level's: the one-copy scale calibrated from men underestimates a partial loss). The band deviation d at the X heterozygous calls (`baf.py`: the maximum-likelihood location of the
two-band mixture 0.5 +- d with the sample's own band SD) gives the second estimate of the fraction: loss f = 4d/(1+2d),
gain f = 4d/(1-2d); an extra X in part of the cells of an XY sample shows heterozygous-looking calls with
d = (1-f)/(2(1+f)). *Measured* on the mock: 47,XXY, 45,X, a mosaic loss of Y at 0.50 (read 0.49-0.51) and a mosaic loss of an
X at 0.30 (read 0.27-0.28 from the level, 0.31 from the BAF) are recovered; every other sample reads XX or XY as planted.

The check against the reported sex (sample sheet or pedigree) says what a disagreement means: "XX in a reported male (an XX
male, or a sample swap)", "XY in a reported female (a sample swap, or an XY female)", "XXY in a reported female (47,XXY; a
sex called from X heterozygosity reads two X's as female)", "X in a reported male (45,X; a sex called from X heterozygosity
reads one X as male)" - the last two because a sex inferred from X heterozygosity alone (peddy's) mislabels them. Only an
XX/XY contradiction raises the `sex_mismatch` flag (a swap candidate); an aneuploidy is a finding.

## Pedigree (`pedigree.py`)

**Probes**: autosomal SNV probes with `keep_gt`, a cohort call rate of 0.98 or more and a minor allele frequency of 0.05 or
more, thinned to one per 50 kb (*measured* on the default mock: 11,370 of 82,000; projected about 50,000 on the example).
A distance thinning stands in for LD pruning at this density: halving it to 25 kb moves no related pair's kinship by
more than 0.02 (tested on the mock; to be repeated on the example's known trios).

**Estimator**: KING-robust kinship, the between-family estimator of Manichaikul 2010 (equation 11), as PLINK 2's
`--make-king` implements it: phi = (2 N_AaAa - 4 N_AAaa - N_Aa(i) - N_Aa(j)) / (4 min(N_Aa(i), N_Aa(j))) + 1/2, with
N_AaAa the probes at which both samples are heterozygous, N_AAaa the probes at which they are opposite homozygotes (IBS0)
and N_Aa(i) sample i's heterozygous calls, all over the probes called in both. *Measured* against `plink2 --make-king-table`
on the same probes of the mock: kinship within 0.0001 on every related pair, IBS0 identical. The simplified form
(N_AaAa - 2 N_AAaa) / (2 min N_Aa), which some tutorials quote, differs by up to 0.03 when the two heterozygote counts differ
(a child with a uniparental disomy has fewer). Every pair at once by matrix products of the 0/1 indicator matrices
(heterozygous, homozygous-reference, homozygous-alternate, called), accumulated over blocks of 8,000 probes (BLAS;
memory one block). IBS0, the share of probes at which the pair is opposite homozygous, separates parent-offspring (at the
genotyping-error floor: the lowest IBS0 among first-degree pairs, 0.0009 on the mock) from full siblings (a quarter of the
unrelated expectation, the mean of 2 p^2 q^2 over the probes); the cut is their midpoint.

**Classes** (KING's boundaries, each a parameter): kinship above 0.354 duplicate or monozygotic; 0.177-0.354 first degree,
split into parent-offspring and full sibling by IBS0; 0.0884-0.177 second degree; 0.0442-0.0884 third degree.

**Families**: a trio is a sample with two relatives (parent-offspring class; a full-sibling- or second-degree-classed pair
is tried too, see below) who are unrelated to each other and, where both sexes are known, of opposite sex, and whose three
genotypes show a Mendelian-error rate at the error floor in that orientation and far above it in the other two (at most
0.02, the others at least three times it, and errors standing out on at most two chromosomes) - the orientation test
resolves who is the child without ages. The floor is then the cohort's own: the median rate of the trios accepted, and a
trio above five times it (never under 0.003) is dropped - on the example true trios sit at 0.0003 and a grandfather in a
father's place gave 0.019, under the absolute cap but sixty times the floor. The rate is taken
per chromosome and its median used, so that one chromosome's errors do not hide a true trio: the mock's isodisomy-15 child
reads as a second-degree relative of its father by kinship (0.159) and IBS0 (0.020), because chromosome 15 carries no
paternal copy; the orientation test passes on the other chromosomes (0.0035), the pair is reclassified parent-offspring with a
note, and the trio's note names chromosome 15 (0.18 errors against 0.0035 elsewhere) as a uniparental disomy, a deletion in
the child, or a chromosome not inherited as expected - the family-aware reading phase 2 will make. A parent-offspring pair
with no third member is a duo, oriented by the given pedigree or a sample sheet age column when available, else reported as
an unoriented pair; full-sibling pairs and duplicate pairs are listed as families of their own kind.

**The check against a given pedigree** (`pedigree/check.tsv`; a table with `Siblings`, `Second Order` and `Third Order`
columns, as the IGSR's, has those pairs counted as disclosed): per claimed parent, agrees / the claimed parent is unrelated
(a sample swap when the parent's sex check says so, or the parent has a duplicate; else non-paternity or non-maternity, or a
swap elsewhere) / the claimed parent is a relative of another degree; per discovered pair not in the pedigree, an undisclosed
parent-offspring, full-sibling or second-degree pair, or an undisclosed duplicate (siblings with the same given parents
agree); per member, a sex contradiction between the pedigree's sex and the complement found. *Measured* on the mock: the
swap (the claimed father unrelated and XX, the true father found as an undisclosed parent-offspring pair under the swapped
name, and the trio assembled with him), the non-paternity, the duplicate and all three recoverable trios and three duos, with
no spurious first-degree pair. QC flags added: `duplicate` on both members of an undisclosed duplicate pair,
`unexpected_relative` on an undisclosed first- or second-degree pair when a pedigree was given.

**Hardy-Weinberg, second pass**: after the pedigree, one member of every related pair (third degree or closer) is dropped and
the exact test recomputed over the unrelated passing samples; `qc/probes.tsv` is updated.

## Events (`events.py`, `windows.py`, `segment.py`, `roh.py`, `karyotype.py`)

Every sample that does not fail QC hard, every chromosome, from the corrected LRR, the BAF and the genotypes; a flagged
sample is read with its flags beside every call.

**Tracks per window** (`windows.py`): windows of about 100 probes pooled by probe count over the probes kept for LRR
(a similar amount of evidence everywhere; 1,939 windows on a 196,000-probe chromosome). Per window and sample: the median
corrected LRR; the band deviation d at heterozygous calls, the location of the symmetric two-band mixture 0.5 +- d fitted by
maximum likelihood with the sample's own robust band SD (an EM with a known SD: the E-step gives each folded value's
probability of the upper band, the M-step d; `baf.d_em_windows`) and its log-likelihood ratio against no split; and the
heterozygous rate against the Hardy-Weinberg expectation at the called probes, 2 p (1 - p) from the cohort's allele
frequencies. The mixture carries a 2% uniform outlier component on [0, 0.5]: a mis-called genotype (a homozygous probe
called heterozygous, its BAF near 0 or 1) would otherwise be explained only by a huge split and give one window a
likelihood ratio in the tens, which fragmented every chromosome. The per-window mean of |BAF - 0.5| is biased upward by
noise and is not used. Tracks are computed for blocks of 64 samples at once and stored (`events/windows/`).

**Segmentation** (`segment.py`): the three tracks, standardised by the sample's own window-level noise (LRR from the
first differences of the window medians, robust to events; the band track as sqrt(2 x the window's log-likelihood ratio),
at the null a half-normal, centred on its median and scaled by its robust spread - the deviation d itself has a point mass
at 0 and no usable spread; the heterozygous ratio by its robust spread), are segmented jointly by circular binary
segmentation (Olshen 2004): for every candidate segment (i, j] the statistic is the
sum over tracks of the squared t-statistic of the segment's mean against the rest of the chromosome, computed for all
pairs at once; the best segment is accepted when its statistic exceeds the 3-df chi-square quantile at 0.01 over the number
of candidate pairs (never below 20) - a Bonferroni-style bound in place of permutations, which the cohort scale forbids - and
the pieces are scanned again. Adjacent segments whose means do not differ (half that bound) are merged. Every boundary is
refined to probe resolution by the least-squares split of the probe-level LRR and folded BAF around the window boundary,
with a tie-break towards the largest probe gap (a boundary in a centromere or a probe desert is placed at the gap). On the X of
an XY sample the baseline is the sample's own X level and the heterozygous track is not used. *Measured*: 100 null
chromosomes of 400 windows receive no breakpoint; a planted 300-window segment (1.5 and 2 SD on two tracks) is found
within two windows.

**The sample's baseline**. A real array gives every sample a baseline band split: the heterozygous cluster of each probe
sits a little off 0.5 in its own way, so the symmetric two-band mixture fitted to any long segment finds a d of 0.01-0.03
even where nothing is there (*example*: 0.014 for a typical sample), and thousands of heterozygous calls make that
"significant" against d = 0. So before typing, d is fitted per chromosome arm over the genome (most arms carry no event):
the median is the sample's baseline split d0, the robust spread its noise; a segment's split must exceed d0 by 0.02 and by
three times that spread, its log-likelihood ratio is taken against d0, and the cell fraction comes from the split beyond
the baseline, sqrt(d^2 - d0^2). Likewise the median corrected LRR per chromosome and the robust spread of those medians
across the chromosomes (0.005-0.008 on the example) is the noise a whole-chromosome or arm level has to beat.

**Typing** from the probe-level evidence of each segment against the baseline (LRR 0, d at the sample's d0, heterozygous
ratio 1):

| reading | rule |
| --- | --- |
| run of homozygosity | heterozygous ratio at or under 0.1 over 40 or more called probes (a true run sits at 0.003-0.03 on the example; 0.2 let two runs join across normal sequence between them) (the calls left are errors, whatever their BAF says); LRR at one copy (under half the one-copy level) makes it a constitutional deletion, LRR up a gain without heterozygosity; a run under 10 Mb is a run (roh.all.tsv), not an event |
| gain / loss | the median LRR at least 0.025 and 5 standard errors from 0 (the median's SE from the sample's own probe noise), with the band deviation significant too (beyond the baseline as above, a likelihood ratio of 30 or more over 20 or more heterozygous calls) unless the shift reaches 0.10 on its own - residual GC waves of a few hundredths survive the correction and must not become events - or the segment is a whole chromosome or arm whose level stands four robust SD out among the sample's chromosomes (never under 0.03) |
| copy-neutral LOH | LRR flat, and - with no LRR change to back it - a split of at least 0.04 beyond the baseline and four times its spread (a cell fraction of 0.08), a likelihood ratio of 50 or more over 50 or more heterozygous calls; the brief's floors (0.10 for a whole chromosome, 0.20 for a 2-10 Mb stretch) lie above it |
| nothing | everything else (a segment the noise cut) |

On the single X of a man the same rules hold with no BAF reading; the X-transposed region (Xq21.31-q21.32, 99% identical
to a stretch of Yp, so that its probes read the X and the Y copies in a man: two copies, heterozygous-looking calls) is an
artefact region and a man's segment lying mostly inside it is left alone.

Adjacent segments typed as nothing are joined without condition (a boundary between two nothings is noise), so that a
chromosome whose every piece is a little high reads as one segment and is typed as a whole; a typed segment takes in an
untyped neighbour at the same level when the union types the same (a sliver the scan cut off its event); adjacent segments
typed alike, over any short untyped piece between them (under 5 Mb or 15% of the joined span), are joined when their levels
agree (within 0.03, or within 0.1 when the BAF puts both pieces in the same share of cells, within 0.1 of each other: on a
real array the LRR level of one 26-Mb mosaic loss wanders by 0.05 along it) and the joined segment types the same. Floors: 2 Mb (a segment's extent counts half of each gap around its first and
last probe, where the change of copy number lies) and 40 probes (`--min-mb`). A segment covering 90% or more of a
chromosome's probed span is `whole`, of an arm's `p` or `q`, else `stretch`; a whole-X gain or loss is the sex step's call and
is not repeated.

**Cell fraction**: the segment's d is fitted on all its heterozygous calls (a subsample of 4,000 beyond that) with the
sample's robust band SD - the band SD is not fitted jointly with d, because two bands close together and one wider band are
nearly the same mixture and a joint fit on a long segment reports small splits that are not there - and read beyond the
sample's baseline split. From the band deviation, the primary estimate - gain f = 4d/(1-2d), loss f = 4d/(1+2d), copy-neutral LOH
f = 2d (the band positions (1+f)/(2+f), 1/(2-f) and (1+f)/2 solved for f); and from the LRR on the cohort's copy scale, the
secondary - the one-copy level from the male X against the female X (sex step; -0.53 to -0.55 on the mock, which simulates
-0.55), the three-copy level from the median LRR of any constitutional trisomy found (0.33 on the mock, which simulates
0.35), else PennCNV's +0.40 scaled by the cohort's compression. Both are reported; a disagreement above 0.2 is noted. An
event with f at or above 0.9 is constitutional. *Measured* on the default mock (60 samples, 12-kb spacing, two seeds): the
constitutional trisomy 21 (f 0.94-0.96 from the BAF), the mosaic gain of 12 at 0.30 (0.27-0.28), the mosaic loss of 7q at
0.25 (0.23), the mosaic copy-neutral LOH of 11p at 0.40 (0.39-0.40), the 10-Mb 1p36 deletion (one copy, no heterozygosity),
the 15-Mb run of homozygosity, the whole-chromosome isodisomy of 15 (`upd(15)(iso)?`) and the four sex-chromosome
events are each found with their type and class, and no null sample receives an event. The detection curve over cell
fraction and size at the example array's probe spacing is in `docs/calibration/`.

**Runs of homozygosity** (`roh.py`): PLINK's sliding window (Purcell 2007, `--homozyg` defaults): over the called SNV
probes a window of 50 is homozygous with at most one heterozygous call; a probe is in a run when 5% or more of the windows
covering it are; runs of 1.5 Mb or more with 100 or more probes, no gap above 1 Mb and at most 50 kb per probe are kept
(`events/roh.all.tsv`); the genome fraction in runs (F_ROH, over the autosomal probed span) is reported per sample; runs of
10 Mb or more are events (joined to a segmentation call of the same place when there is one). A whole chromosome without
heterozygosity is an isodisomic uniparental disomy or consanguinity: both are said, and the trio addition of phase 2 will
decide (the pedigree step already names the chromosome whose Mendelian errors stand out). The ROH scan is not run on a
single X (an XY sample, or a 45,X).

**Sex chromosomes**: a whole X gained or lost and a Y lost or gained, in every cell or in part, enter from the sex step with
its cell fractions; a constitutional aneuploidy is read against the reported sex (a 45,X is an X lost in a reported
female, a Y lost in a reported male), else against the complement the Y implies. **Flags and tags**: every event carries
its evidence (LRR, its SE, d, its log-likelihood ratio, probes, heterozygous calls, the two cell fractions) and a note where
the evidence is thin (under 20 heterozygous calls, under 100 probes), where the sample is flagged (its QC flags travel with
the event), and the names of the artefact-prone regions it mostly lies in (the MHC, the 8p23.1 and 17q21.31 inversions, the
X-transposed region, the centromere flanks, the chromosome ends, an acrocentric short arm) - tagged, never removed from
the table; the digest leaves out a copy change read from the LRR alone in such a region (on the example these are the
commonest residual calls: a few hundred probes at the centromere flanks, +0.10 to +0.13, no heterozygous calls). **The cohort's own artefact map**: a 1-Mb bin (on an autosome or the X) where an event of the same type is called in 2% or
more of the samples read (three at least) is recurrent - a probe cluster the array reads oddly, or a common copy-number polymorphism (*example*: the
centromere flanks of chromosomes 5 and 10 in 14% of samples, 16p11.2 in 13%, the 9q12 heterochromatin block of CHM13 in
2.3%) - and an event whose bins are mostly recurrent is tagged `recurrent (N samples)`, never removed from the table. An
event whose probes lie five times further apart than the chromosome's median spacing is tagged `probe-sparse` (*example*:
450 probes over the 42 Mb of CHM13's 9p11-q21 heterochromatin and segmental duplications read as a 44% "gain" in 48
samples, from the LRR alone). The digest leaves out a copy change read from the LRR alone in a recurrent, probe-sparse or
artefact-prone region, and any mosaic copy change the BAF split does not confirm. A sample whose event count exceeds the cohort median by 3 robust SD
is flagged `event_count`.

**Karyotype string** (`karyotype.py`): the events with evidence (a copy change read from the LRR alone in a tagged region
stays in the table only; a plain XX or XY against the reported sex is the sex check's finding, not a term), in triokaryo's grammar - the modal number and complement, the constitutional terms in
chromosome order (`+21`, `del(1)(p36.32p36.21)`, `roh(1)(q21.2q24.1)`, `upd(15)(iso)?`), and one line per mosaic with its
cell fraction: `mos 47,XY,+12[0.28]/46,XY`, `mos 46,XY,del(7)(q11.21q36.3)[0.23]/46,XY`, `mos 46,XX,loh(11)(p15.5p11.12)[0.40]/46,XX`,
`mos 45,X[0.51]/46,XY` for a mosaic loss of Y; bands from the vendored tables.

**Figures**: for every sample with an event or a flag (`--figures all` for every sample), the genome-wide figure (every
chromosome: corrected LRR with the one- and three-copy levels and the called segments, BAF, the heterozygous ratio with the
runs of homozygosity; events as labelled bars) and one figure per chromosome with an event (ideogram, LRR with the
segment's level, BAF with the bands expected at the called cell fraction, the heterozygous ratio). Each image carries its
title and legend; a text sidecar goes beside it.

**Performance** (*measured*, 2,000 samples x 196,008 probes, one chromosome, 8 threads): window tracks 65 s (one
chromosome, so one thread); calls 0.30 s of one core per sample-chromosome of 1,939 windows after vectorising the
segmentation scan and the merge and fitting d on at most 4,000 heterozygous calls (the SE of d is then under 0.001).
Projected for the example: about 10 minutes with 8 threads.

## Report (`report.py`)

`report/index.html`, one self-contained page per cohort (every figure embedded as base64, readable from a file:// URL):
tiles (samples; passing; flagged by reason; families; samples with events; samples not plainly XX or XY), the QC section
(the metric distributions with the thresholds, the GC wave before and after, the thresholds with their sources, the sortable
and filterable sample table), the sex section (the X-against-Y figure, the table of every sample not plainly XX or XY or
disagreeing with its reported sex), the pedigree section (families, the check against the given pedigree with the
disagreeing rows, the kinship-against-IBS0 figure), the events section (the landscape - events per chromosome, one row per
sample with an event - and the sortable table with every number and the event in words, linked to the sample pages), the ROH
section (F_ROH, the runs of 5 Mb or more), and the samples set aside (hard QC failures; noise, event-count and duplicate
warnings), like triokaryo's cohort guards. Per-sample pages (`report/samples/<sample>.html`) for flagged samples and samples
with events (`--pages all` for every sample): the karyotype string, the QC row against the cohort medians, the complement
and its check, the families, every event in plain words with its evidence, the figures. The guide (`report/guide.html`)
explains every column and colour, the calling rules and the thresholds with their sources, with pattern cards of what each
event kind looks like on an array (a constitutional trisomy, a mosaic gain, a mosaic loss, a mosaic copy-neutral LOH, a run of
homozygosity, 47,XXY, a mosaic loss of Y), drawn from idealised tracks with array noise. The digest (`report/digest.html`,
`digest.tsv`, `digest_left_out.tsv`) lists the largest events for hand-off, ranked as triokaryo's digest ranks: a whole
chromosome, then an arm, then a stretch by size; an event in every cell before one in part of the cells; a change of copy
number before a copy-neutral LOH before a run of homozygosity; the floors (stretches of 10 Mb or more, mosaics of 10% or
more for an arm or a stretch, runs of homozygosity of 20 Mb or more, samples that pass QC) are parameters and every event
left out is counted with its reason; each listed event carries its facts and the sample's genome-wide and chromosome
figures. `arraykaryo run` chains extract, qc, sex, pedigree, events and report with one log, stops at a failed step, and
prints where every table and page is. *Measured*: the default mock (60 samples, 82,000 probes) runs end to end, with
figures and the report, in 17 s on four threads.

## References

- Peiffer DA et al. (2006) High-resolution genomic profiling of chromosomal aberrations using Infinium whole-genome genotyping. Genome Res 16:1136 (LRR and BAF).
- Wang K et al. (2007) PennCNV: an integrated hidden Markov model designed for high-resolution copy number variation detection in whole-genome SNP genotyping data. Genome Res 17:1665.
- Colella S et al. (2007) QuantiSNP: an Objective Bayes hidden-Markov model to detect and accurately map copy number variation using SNP genotyping data. Nucleic Acids Res 35:2013.
- Diskin SJ et al. (2008) Adjustment of genomic waves in signal intensities from whole-genome SNP genotyping platforms. Nucleic Acids Res 36:e126.
- Manichaikul A et al. (2010) Robust relationship inference in genome-wide association studies. Bioinformatics 26:2867 (KING).
- Jacobs KB et al. (2012) Detectable clonal mosaicism and its relationship to aging and cancer. Nat Genet 44:651.
- Laurie CC et al. (2012) Detectable clonal mosaicism from birth to old age and its relationship to cancer. Nat Genet 44:642.
- Forsberg LA et al. (2014) Mosaic loss of chromosome Y in peripheral blood is associated with shorter survival and higher risk of cancer. Nat Genet 46:624.
- Machiela MJ et al. (2016) Female chromosome X mosaicism is age-related and preferentially affects the inactivated X chromosome. Nat Commun 7:11843.
- Loh PR et al. (2018) Insights into clonal haematopoiesis from 8,342 mosaic chromosomal alterations. Nature 559:350 (MoChA).
- Purcell S et al. (2007) PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 81:559 (`--homozyg`, `--hwe`).
- Narasimhan V et al. (2016) BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data. Bioinformatics 32:1749.
- Pedersen BS, Quinlan AR (2017) Who's who? Detecting and resolving sample anomalies in human DNA sequencing studies with peddy. Am J Hum Genet 100:406.
- Wigginton JE, Cutler DJ, Abecasis GR (2005) A note on exact tests of Hardy-Weinberg equilibrium. Am J Hum Genet 76:887.
