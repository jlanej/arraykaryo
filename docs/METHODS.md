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
| Hardy-Weinberg exact test (Wigginton 2005) | p < 1e-10 | PLINK `--hwe` practice for arrays; over the passing samples here, over the unrelated passing samples after the pedigree step |
| `GenTrain_Score` | < 0.6 | Illumina's GenTrain guidance |
| `Cluster_Sep` | < 0.3 | Illumina: overlapping clusters |
| PAR probe | kept, tagged `par` | used for the autosomal-like X only |

On the Y the call rate is taken among the samples with a Y (its PAR probes among everyone). Every threshold is a parameter
(`THRESHOLDS` in `qc.py`, written with its source to `qc/summary.json`).

**Sample QC** (`qc/samples.tsv`; each metric's cohort median and robust SD in `qc/summary.json`):

| metric | definition | flag | source |
| --- | --- | --- | --- |
| call rate | called / autosomal SNV probes kept for genotypes | < 0.97 | Illumina's and MoChA's recommendation |
| LRR SD | SD of the corrected autosomal LRR | > 0.20 flag; > 0.30 hard fail | MoChA's mosaic practice; PennCNV `filter_cnv --qclrrsd 0.3` |
| BAF SD | SD of BAF - 0.5 at heterozygous calls, the median over the autosomes of the per-chromosome value (so one chromosome's split bands do not flag a sample) | > 0.06 | the band SD that blurs the bands (MoChA's `baf_auto` flags contamination at 0.03) |
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
   quadratic, fitted by least squares over the autosomes on the LRR clipped to +-1 (robust to events and outliers), and the
   fit subtracted. *Measured* on the mock: the correlation of LRR with the window GC falls from a median of 0.04-0.09 before
   to 0.000 after; the autosomal LRR SD from 0.175 to 0.150, the planted per-sample noise.
2. *The cohort reference*: the median corrected LRR per probe over the passing samples, per batch when every batch holds
   20 or more passing samples, else cohort-wide, subtracted. This removes the array's shared structure (the per-probe offsets
   every sample carries) exactly as triokaryo's panel does for sequencing. On the X the reference is built from the samples
   without a Y, so that one X reads as one copy on the cohort's scale; the Y keeps its GC-corrected raw level against the
   male cluster (a reference over a mixed cohort would sit between the sexes).
3. *BAF*: clipped to [0, 1] at use; the heterozygous-call BAF is re-centred per probe by the cohort's median at that probe
   when its shift exceeds 0.01 (MoChA's cluster adjustment), and a probe whose shift exceeds 0.05 loses its genotype.

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
0.5 the nearer complement), XX without. An X deviating from the anchor's count by a cell fraction of 0.10 or more is a
mosaic (`mos 45,X[0.30]/46,XX`, `mos 47,XXY[0.20]/46,XY`), by 0.90 or more a constitutional aneuploidy (47,XXY, 45,X, 47,XXX,
48,XXYY ...); a Y under one copy in a man by 0.10 or more is a mosaic loss of Y (`mos 45,X[0.50]/46,XY`), above one by
0.10 or more a gain (47,XYY); a Y at 0.25 copies or more in an XX genome is "a Y in part of the cells (46,XX/47,XXY mosaic, or
an XX/XY mixture)". The band deviation d at the X heterozygous calls (`baf.py`: the maximum-likelihood location of the
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
0.02, the others at least three times it) - the orientation test resolves who is the child without ages. The rate is taken
per chromosome and its median used, so that one chromosome's errors do not hide a true trio: the mock's isodisomy-15 child
reads as a second-degree relative of its father by kinship (0.159) and IBS0 (0.020), because chromosome 15 carries no
paternal copy; the orientation test passes on the other chromosomes (0.0035), the pair is reclassified parent-offspring with a
note, and the trio's note names chromosome 15 (0.18 errors against 0.0035 elsewhere) as a uniparental disomy, a deletion in
the child, or a chromosome not inherited as expected - the family-aware reading phase 2 will make. A parent-offspring pair
with no third member is a duo, oriented by the given pedigree or a sample sheet age column when available, else reported as
an unoriented pair; full-sibling pairs and duplicate pairs are listed as families of their own kind.

**The check against a given pedigree** (`pedigree/check.tsv`): per claimed parent, agrees / the claimed parent is unrelated
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
