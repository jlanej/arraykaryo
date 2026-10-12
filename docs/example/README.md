# The public example: a 1000 Genomes genotyping-array export

The example is a public 1000 Genomes Project cohort genotyped on an Illumina array, exported with `bcftools +gtc2vcf` to
one BCF (2,141 samples, about 2.4 million probes, T2T-CHM13v2.0 coordinates in the copy used here; other projects are on
GRCh38 or GRCh37, which arraykaryo recognises from the contig lengths). The data are not in this repository and no
location is written here: point `ARRAY_BCF` at your copy. The known 1000 Genomes pedigree (the 3,202-sample pedigree
table published by the IGSR) is the reference for the pedigree step; point `KNOWN_PED` at your copy of it.

```bash
export ARRAY_BCF=/path/to/your/1000g_array_export.bcf       # indexed (.csi beside it)
export KNOWN_PED=/path/to/your/1000g_pedigree.txt            # IGSR's 3,202-sample pedigree table
arraykaryo inspect "$ARRAY_BCF" --json inspect.json
```

What `inspect` reports on this export, for the record: FORMAT GT, GQ, IGC, BAF, LRR, NORMX, NORMY, R, THETA, X, Y (only
GT, LRR and BAF are used); INFO ALLELE_A/ALLELE_B, INTENSITY_ONLY, GC, GenTrain_Score, Cluster_Sep and the cluster
statistics; diploid GT everywhere, including a male's X and Y; a no-call is `./.`; a missing LRR or BAF is `nan`; BAF is
not clipped (values slightly below 0 or above 1 occur); on the Y a female's probes are `./.` with LRR around -3 to -5.

Once the pipeline exists (PROMPT.md), the example run is:

```bash
arraykaryo run --bcf "$ARRAY_BCF" --pedigree "$KNOWN_PED" --out out/1000g --threads 8
```

and the things to check in `out/1000g/report/index.html`:

- the pedigree inferred from the genotypes recovers the known 1000 Genomes trios and duos present on the array (for
  instance the CEU trio NA12878 with NA12891 and NA12892, the YRI trio NA19240, the PUR trio HG00733, the CHS trio
  HG00702, the MXL trio NA19675), with no undisclosed first-degree pair unexplained;
- the sex-chromosome complement agrees with the pedigree's sex except where it names an aneuploidy or a mosaic, each
  of which is listed with its evidence (X heterozygosity, X and Y levels);
- the events found on samples that triokaryo read from whole-genome sequencing agree: HG01103 (mosaic losses on 2q and
  14q) and NA12748 (a mosaic gain of 13q, a terminal loss of 21q), whose parents or children are also on the array;
- lymphoblastoid-cell artefacts (gains of chromosome 12, losses of the X in women, other mosaics acquired in culture)
  are listed as findings with their cell fraction, not treated as errors.
