# Benchmarks

Measured times per version, on the mock unless stated (a laptop: Apple silicon, 10 cores, local SSD). "Example" is the public
1000 Genomes export of `docs/example/README.md` (2,141 samples, 2.42 million probes, 132 GB BCF), which is not in the
repository.

| version | step | scale | threads | wall time | peak RSS |
| --- | --- | --- | --- | --- | --- |
| 0.2.0 | mock | 60 samples x 82k probes (8 chromosomes) | 1 | 10 s | - |
| 0.2.0 | mock | 2,000 x 196k (chr1) | 1 | 12 min | - |
| 0.2.0 | extract, route bcftools (memory-mapped writer) | 2,000 x 196k | 1 | 164 s | 4.0 GB |
| 0.2.0 | extract, route cyvcf2 (memory-mapped writer) | 2,000 x 196k | 1 | 33 s | 3.0 GB |
| 0.2.0 | extract, route cyvcf2 (positioned writes, 8,192-record blocks; the default) | 2,000 x 196k | 1 | 41 s | 0.49 GB |
| 0.2.0 | `bcftools query` alone, the streaming floor of route (a) | 2,000 x 196k | 1 | 41 s | 10 MB |
| 0.2.0 | `bcftools query` alone, example chr21 (32,823 records x 2,141) | - | 1 | 10 s | 20 MB |
| 0.2.0 | extract | 60 x 82k | 4 | 1.3 s | - |
| 0.2.0 | qc (two passes, figures) | 60 x 82k | 4 | 4 s | - |
| 0.2.0 | sex (figures for 6 samples) | 60 x 82k | 1 | 1.7 s | - |
| 0.2.0 | pytest (19 tests, the mock built once) | - | - | 22 s | - |
| 0.3.0 | pedigree (11,370 probes, 60 samples, the plink2 cross-check apart) | 60 x 82k | 1 | 1.3 s | - |
| 0.3.0 | pytest (23 tests) | - | - | 29 s | - |
| 0.4.0 | qc (two passes, figures) | 2,000 x 196k (chr1) | 8 (one chromosome: one used) | 85 s | 3.5 GB resident, file-backed (released block by block from 0.4.0) |
| 0.4.0 | pedigree (4,684 probes on chr1; the HWE exact test over 196k probes is most of it) | 2,000 x 196k | 1 | 19 s | 1.4 GB |
| 0.4.0 | events: window tracks | 2,000 x 196k | 1 (one chromosome) | 65 s | - |
| 0.4.0 | events: calls, before the vectorised scan | 2,000 x 196k | 8 | 188 s (0.73 core-s per sample) | 2.7 GB |
| 0.4.0 | events: one sample's calls on chr1 (1,939 windows), after | - | 1 | 0.30 s | - |
| 0.4.0 | events with figures for 22 samples | 60 x 82k | 4 | 5.7 s | - |
| 0.4.0 | pytest (28 tests) | - | - | 47 s | - |
| 0.5.0 | run (extract to report, figures for 22 samples, 22 sample pages) | 60 x 82k | 4 | 17 s | 0.28 GB |
| 0.5.0 | pytest (30 tests, the end-to-end run included) | - | - | 65 s | - |
| 0.5.2 | extract (route cyvcf2) | example: 2,141 x 2.42M, 24 chromosomes | 8 | 390 s | 1.46 GB (the cache 35 GB) |
| 0.5.2 | qc (two passes, figures) | example | 8 | 373 s | 4.1 GB (the probe table was held in memory; streamed from 0.5.2) |
| 0.5.2 | sex (figures for 51 samples) | example | 1 | about 90 s | - |
| 0.5.2 | pedigree (48,784 probes, 331 trios, the HWE screen) | example | 1 | 134 s | 2.6 GB |
| 0.5.3 | events: window tracks | example | 8 | 546 s | - |
| 0.5.3 | events: calls on 2,124 samples (6.2 core-s per sample over 23 chromosomes) | example | 8 | 1,406 s | 5.0 GB |
