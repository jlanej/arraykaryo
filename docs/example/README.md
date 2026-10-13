# The public example: a 1000 Genomes genotyping-array export

The example is a public 1000 Genomes Project cohort genotyped on an Illumina array, exported with `bcftools +gtc2vcf` to
one BCF (2,141 samples, 2,421,634 probes, T2T-CHM13v2.0 coordinates in the copy used here; other projects are on GRCh38 or
GRCh37, which arraykaryo recognises from the contig lengths). The data are not in this repository and no location is written
here: point `ARRAY_BCF` at your copy. The known 1000 Genomes pedigree (the 3,500-sample table published by the IGSR, with
its `Siblings`, `Second Order` and `Third Order` columns) is the reference for the pedigree step; point `KNOWN_PED` at your
copy. A sample sheet with the reported sex (`sample`, `reported_sex`) was made from the same table.

```bash
export ARRAY_BCF=/path/to/your/1000g_array_export.bcf       # indexed (.csi beside it)
export KNOWN_PED=/path/to/your/1000g_pedigree.txt            # the IGSR pedigree table
arraykaryo inspect "$ARRAY_BCF" --json inspect.json
arraykaryo run --bcf "$ARRAY_BCF" --out out/1000g --threads 8 --sample-sheet samples.tsv --pedigree "$KNOWN_PED"
```

What `inspect` reports on this export, for the record: FORMAT GT, GQ, IGC, BAF, LRR, NORMX, NORMY, R, THETA, X, Y (only
GT, LRR and BAF are used); INFO ALLELE_A/ALLELE_B, INTENSITY_ONLY, GC, GenTrain_Score, Cluster_Sep and the cluster
statistics; diploid GT everywhere, including a man's X and Y; a no-call is `./.`; a missing LRR or BAF is `nan`; BAF is
not clipped (values slightly below 0 or above 1 occur); on the Y a woman's probes are `./.` with LRR around -3 to -5 outside
the pseudoautosomal regions, and the PAR probes of the Y carry diploid calls in both sexes.

## What was checked (version 1.0.0; the numbers below are from this run)

**Timing** (a laptop, 8 threads): extract 390 s (a 35-GB cache, 1.5 GB resident), qc 373 s, sex about 90 s, pedigree
134 s, events 35 minutes without figures (the window tracks 532 s of it, the calls 1,552 s on 8 threads over 23
chromosomes; 5 GB resident), report with `--pages digest` 18 s (44 sample pages, the figures of 11 samples drawn; the digest
page 18 MB with its figures embedded, the cohort page 5.5 MB). `benchmarks/`.

**Quality** (`qc/`). 1,788 of 2,141 samples pass with no flag; 17 fail hard (LRR SD above 0.30); flagged: 257 for a BAF SD
above 0.06 (this platform's noisy tail: their LRR SD is 0.18 against 0.14 for the rest, and the BAF SD is continuous from
its median of 0.038 to 0.084, so a fixed 0.06 flags 12% of this cohort - they are read, with the flag beside every call),
130 for an LRR SD above 0.20, 15 for a wave factor above 0.05, 6 for an LRR mean beyond 0.1. Cohort medians: call rate
0.997, LRR SD 0.149 raw and 0.142 corrected, BAF SD 0.038 (robust band SD 0.035), waviness factor -0.007; the correlation
of LRR with the 1-Mb GC falls from a median absolute 0.02 (robust SD 0.18 across samples) to 0.007 (SD 0.004). Probes:
2,044,352 of 2,421,634 carry no flag; 108,400 lose their genotype for a cluster separation under 0.3, 11,452 for a GenTrain
score under 0.6, 11,576 for a call rate under 0.95, 9,907 as duplicate positions, 640 for a shifted heterozygous cluster,
345 as non-SNV records, 39 for an LRR median beyond 0.3, 23,632 (the top 1%) for a noisy LRR; 221,583 fail Hardy-Weinberg at
1e-10 in a first pass over every direction, which is the Wahlund effect of a cohort of many populations - only an excess of
heterozygotes (a probe reading two loci) counts, and over the 1,623 unrelated passing samples 400 probes fail that way.

**Sex chromosomes** (`sex/`). The X modes sit at -0.40 (one copy) and 0.01 (two copies): this platform compresses one
copy to -0.40, not the -0.55 of the mock or PennCNV's -0.66, and the cohort's own scale is used. The Y modes sit at 0.01
(present) and -3.92 (absent); 1,209 of the 1,918 non-PAR Y probes are Y-specific (the X-transposed region at 3-5 Mb of
the CHM13 Y is heterozygous in men and called in women at LRR 0, and is excluded by the empirical test). The X and Y
levels of men spread by 0.072 and 0.075 copies (robust SD), the X of women by 0.023, so the mosaic floors are 0.21, 0.22
and 0.10 copies. Complements: XX 1,059, XY 1,033, 47,XXY 2, and the mosaics of lymphoblastoid culture: 12 men with a loss of
Y in 23-85% of cells (two of them, with the Y gone from 78% and 85% of cells and one X without heterozygosity, read as
45,X/46,XY, not as women), 29 women with an X lost in 15-53% of cells (median 24%) confirmed by the BAF split at X
heterozygous calls, 6 women with an X gained in 13-45% of cells, 2 men with a Y gain read from the level alone (tentative:
the Y level varies between Y haplogroups). Against the reported sex, 2,126 of 2,130 agree; the four that do not: two
reported men are XX with fully heterozygous X's (HG01465, HG01896: an XX male or a sample swap), one reported man is 47,XXY
(HG02372) and one reported woman is 47,XXY (HG01683: a sex called from X heterozygosity reads two X's as female).

**Relatedness and families** (`pedigree/`). KING-robust kinship and IBS0 of every pair from 49,639 autosomal probes
(call rate 0.98, MAF 0.05, one per 50 kb): 719 parent-offspring pairs, 26 full-sibling, 64 second-degree, 66 third-degree;
331 trios and 47 duos assembled (the Mendelian-error floor of the trios is 0.00014; a grandfather in a father's place, at
0.019, is rejected by the floor rule). The known trios the brief names are recovered: NA12878 with NA12891 and NA12892;
NA19240 with NA19239 and NA19238; HG00733 with HG00731 and HG00732; NA19675 with NA19679 and NA19678 (and NA19677 with
NA19676 and NA19675); HG01103 with HG01101 and HG01102. The brief's "HG00702 with HG00703 and HG00704" is in fact the
three-generation family HG00656 and HG00657 -> HG00702, HG00701 and HG00702 -> HG00703, which is what the table says and
what the genotypes find. NA12739, whose father NA12748 the brief names, is not on the array; NA12748 and NA12749 are.
Against the IGSR table: 779 claimed relationships agree; 370 children, 12 fathers and 13 mothers of the table are not in
the export; 7 full-sibling pairs and 22 second-degree pairs are found that the table does not list (the CHS samples
HG00501, HG00512, HG00524 and their relatives, HG00578/HG00581/HG00635, HG02429/HG02479, NA18914-NA19240 and others: the
cryptic relatedness of the 1000 Genomes collection), one pair listed as third order is full siblings (NA19331, NA19334),
and one sample's sex differs between the IGSR table and the sample sheet (HG02300). No claimed parent is unrelated: the
collection's pedigree holds.

**Events** (`events/`, `report/`). 7,094 events on 1,795 of the 2,124 samples read (median 2 per sample; 179 samples
with more than 10 are flagged `event_count`), 23,102 runs of homozygosity of 1.5 Mb or more (F_ROH median 0.007; 10 samples
above 0.0625, the level of a first-cousin union), 356 runs of 10 Mb or more. Of the 5,343 autosomal events 4,539 carry a
tag: 4,355 recurrent (the same type in 2% or more of the samples over the same 1-Mb bins: the centromere flanks of
chromosomes 5 and 10 in 14% of samples, 16p11.2 in 13%, the 9p11-q21 heterochromatin block of CHM13 in 48 samples),
4,482 probe-sparse, 4,187 in a listed artefact region - 4,390 of them read from the LRR alone, with no heterozygous calls to
split; they stay in `events.all.tsv` with their tags and off the digest. On the X, 1,508 of the 1,703 segments carry a tag
(the centromere flank of the X in men, read from the LRR alone); 48 events are the sex step's. The 804 untagged autosomal
events on 371 samples: 428 with the BAF split behind them (139 gains, median 22% of cells; 77 losses, median 17%; 212
copy-neutral LOH, median 12%), 23 copy changes in every cell or near it read from the LRR alone (the 10.9-Mb deletion of
18q22.2q23 in NA06984 among them) and 353 runs of homozygosity of 10 Mb or more. Seven whole autosomes: a gain of
chromosome 12 in 66%, 84% and 49% of the cells of HG01357, HG02322 and NA20759 and of chromosome 9 in 72% and 14% of
NA19208 and NA21143 (the lymphoblastoid culture's commonest gains), a whole 21 without heterozygosity in NA20509 (an
isodisomy or consanguinity: F_ROH says which) and a copy-neutral LOH of 21 in 14% of the cells of HG01461, plus the sex
chromosomes.

The samples triokaryo read from whole-genome sequencing (docs of that project): **HG01103** carries a mosaic loss of 2q
from q22.1 to q31.2 (141-180 Mb) in 47% of cells and a loss at q33.1 (201-203 Mb) in 55%, and a mosaic loss of
14q13.2-q21.1 (29-35 Mb) in 55% of cells; triokaryo read the 2q loss reaching q33.1 with a subclonal extension (67%) and
14q at 64% from sequencing of a different DNA preparation. **NA12748** carries a mosaic gain of 13q31.2-q34 (87-113 Mb)
in 19% of cells and a mosaic terminal loss of 21q22.3 (41-45 Mb) in 20%, where triokaryo read 69% and 59%: the same
events at the same places, at a third of the fraction - the array's DNA and the sequenced DNA are different preparations
of a lymphoblastoid line, and a clone's share changes with passage. His two pericentromeric calls (chromosomes 5 and
10) carry the recurrent and probe-sparse tags. HG01103's parents HG01101 and HG01102 and NA12748's partner NA12749 carry
none of these events.

Lymphoblastoid culture artefacts are listed as findings with their cell fraction, not treated as errors: the X losses in
women and the losses of Y in men above, gains of whole or near-whole chromosomes in part of the cells (chromosome 12 in
three lines at 49-84% and 12q or 12p in four more; chromosome 9 in NA19208 at 72% and NA21143 at 14%; 13q in NA20533 at
78%), and copy-neutral LOH of stretches of 4q (median 9 Mb, varying ends) in 8-14% of the cells of 45 lines - at the floor
of the LOH rule (a band split of 0.05 against a baseline split three times the cohort's median, in lines mostly flagged for
their BAF SD), a pattern of this platform until one is confirmed; the digest's 10-Mb and 10% floors keep most of them off
it. The digest (`report/digest.html`) opens with the 47,XXY and 45,X samples, the whole-chromosome gains, the isodisomic 21,
the arm-sized gains of 11q, 12p and 12q and the loss of 18p, the arm-sized copy-neutral LOH of 5q, 7q and 12p, and the
constitutional deletion of 18q22.2q23.

Caveats of this platform, written down: the one-copy level is -0.40 (not -0.55), so mosaic fractions read from the LRR
run below those from the BAF; the baseline band split is 0.01-0.03 per sample and sets the floor of what the BAF can
confirm (a gain of a whole chromosome needs about 14% of cells here); the Y level varies between Y haplogroups; CHM13's
assembled heterochromatin and acrocentric arms carry mismapped probes (not used, or tagged); the BAF SD threshold of 0.06
flags 12% of this cohort, which is read with the flag; a run of homozygosity is called at a heterozygous ratio of 0.1 or
under (a true run sits at 0.001-0.05 here), so a run interrupted by a short heterozygous stretch is two runs.
