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
