"""arraykaryo: large chromosomal events in a genotyping-array cohort from its VCF/BCF export (GT, LRR, BAF).

inspect (what an export carries), mock (a simulated export with planted truth), extract (one pass into a cache), qc (probe
and sample QC, the GC wave, the cohort reference), sex (the sex-chromosome complement), pedigree (relatedness, families, the
check against a given pedigree), events (near-chromosomal events, runs of homozygosity, karyotype strings, figures); report
follows
(PROMPT.md at the repository root is the brief)."""
__version__ = "0.4.0"
