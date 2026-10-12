# Output

Every table is TSV with a header; booleans 1/0; missing `NA`; positions 1-based inclusive; sizes in bp in the tables (Mb to
one decimal in the pages); sample identifiers exactly as the export names them. Every figure (`.png`) has a text sidecar
(`.txt`: title, caption, key) beside it. Every command writes under one output directory (`--out DIR`) and marks itself done
in `DIR/<step>/done.json`.

## DIR/cache/ (`arraykaryo extract`)

| file | content |
| --- | --- |
| `manifest.json` | the export path as given, its size and modification time, the md5 of its header, the sample count, per-chromosome record counts, the assembly, the contig names as the file spells them, the INFO fields found, the BAF storage type, the route, the arraykaryo version, the time taken |
| `samples.txt` | the samples in the export's order (one per line) |
| `probes.tsv` | chrom, pos, id, ref, alt, snv (a biallelic single-base change), intensity_only, par, gc, gentrain, cluster_sep, allele_a, allele_b (NA when the export lacks them); every record of the export, in file order within each chromosome |
| `probes.npz` | the same, binary (pos, snv, intensity_only, par, gc, gentrain, cluster_sep) |
| `gt.<chrom>.npy`, `lrr.<chrom>.npy`, `baf.<chrom>.npy` | samples x probes: int8 (0 hom-ref, 1 het, 2 hom-alt, -1 no-call), float32, float16 |
| `done/<chrom>.json` | the chromosome is complete (resumability) |

## DIR/qc/ (`arraykaryo qc`)

`probes.tsv`: chrom, pos, id, keep (the LRR is usable), keep_gt (the BAF and genotype are usable), reason (comma-joined:
intensity_only, non_snv, duplicate_position, excluded, lrr_median, lrr_sd, call_rate, baf_shift, gentrain, cluster_sep,
hwe), par, intensity_only, snv, call_rate (across samples; on the Y among the samples with a Y), lrr_median (raw, across
samples), lrr_sd (robust, across samples), het_baf_median, hwe_p (autosomes), gentrain, cluster_sep, maf (across the
samples), n_called. `probes.npz` holds keep, keep_gt, maf, call_rate, hwe and het_baf_median in cache order.

`samples.tsv`: sample, reported_sex, batch, call_rate, het_rate (autosomal), lrr_mean (raw), lrr_sd_raw, lrr_sd (after
correction), lrr_sd_diff (the noise SD from first differences, robust to events), baf_sd (median over chromosomes),
baf_drift, wf (PennCNV's waviness factor), gc_corr_before, gc_corr_after (the correlation of LRR with the 1-Mb GC), y_lrr_raw
(the raw mean LRR over the non-PAR Y, the provisional Y signal), y_call_rate, sex and sex_check (added by `sex`), pass (no
flag), hard_fail (not to be read), flags (comma-joined: excluded, call_rate, lrr_sd, lrr_sd_fail, baf_sd, baf_drift, wave,
lrr_mean, sex_mismatch; later duplicate, unexpected_relative, event_count). Later steps add columns to this table.

`correction.npz`: the GC-wave coefficients per sample (5 x samples), the batch index per sample, the GC standardisation,
the reference per probe and per batch for every chromosome. `summary.json`: counts, every threshold with its source, the
cohort median and robust SD of every metric, the batches and whether a per-batch reference was used.

`figures/qc_distributions`: the metrics as histograms with the thresholds. `figures/gc_wave`: LRR against the 1-Mb GC
before and after correction for a typical sample and the one with the strongest wave.

## DIR/sex/ (`arraykaryo sex`)

`sex.tsv`: sample, reported_sex, qc_pass, x_lrr (median corrected LRR on the X outside the PARs), y_lrr (median LRR over the
Y-specific probes), x_het_rate, x_het_rel (relative to the autosomal rate), y_call_rate, n_x_called, n_x_het, n_y_probes,
x_copies, y_copies (on the cohort's scale), base (the complement the Y implies: XX or XY), complement (the constitutional
reading: XX, XY, X, XXY, XXX, XYY ...), sex (M/F from the base, for the pedigree step), karyotype (46,XY; 47,XXY;
mos 45,X[0.30]/46,XX ...), mosaic (in words), f (the cell fraction read), f_lrr, f_baf (from the level and from the BAF
split), d_x (the band deviation at X heterozygous calls), check (against the reported sex), note (thin evidence, a QC
failure, an X heterozygosity that an XY sample should not have).

`summary.json`: the X modes and the scale per copy, the Y modes, how the Y-specific probes were chosen, counts by
complement and by check. `figures/x_vs_y`: the X level against the Y level, coloured by X heterozygosity, every sample not
plainly XX or XY labelled. `figures/tracks_<sample>`: the X (LRR, BAF) and Y (LRR) tracks of every sample not plainly XX or
XY.
