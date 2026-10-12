# Calibration: the detection curve on the mock

Mock cohorts at a probe spacing of 1.3 kb (the example array's), per-sample LRR SD 0.15, BAF SD drawn in 0.03-0.05, 2 replicate(s); every (kind, cell fraction, size) planted in its own sample on chr1, 7, 11, 12, 15 (stretches and arms) or chr21/chr15 (whole chromosomes); 40 null samples per replicate. Detection = an event of any type overlapping the planted span by half or more; the cell-fraction error is |called - planted| on the primary estimate.


## gain

Detection rate (type right / class right), and the median cell-fraction error, by cell fraction (rows) and size (columns, Mb; arm = a whole q arm; whole = a whole chromosome).

| f | 2 | 5 | 10 | 25 | arm | whole |
|---|---|---|---|---|---|---|
| 0.05 | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - |
| 0.10 | 50% (50%/50%) err 0.02 | 0% (0%/0%) err - | 50% (50%/50%) err 0.00 | 50% (50%/50%) err 0.01 | 100% (100%/50%) err 0.02 | 100% (100%/100%) err 0.00 |
| 0.15 | 50% (50%/50%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 |
| 0.20 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/50%) err 0.02 | 100% (100%/100%) err 0.01 |
| 0.30 | 100% (100%/100%) err 0.02 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 |
| 0.50 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/0%) err 0.01 | 100% (100%/100%) err 0.01 |
| 1.00 | 100% (100%/100%) err 0.02 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 |

## loss

Detection rate (type right / class right), and the median cell-fraction error, by cell fraction (rows) and size (columns, Mb; arm = a whole q arm; whole = a whole chromosome).

| f | 2 | 5 | 10 | 25 | arm | whole |
|---|---|---|---|---|---|---|
| 0.05 | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - |
| 0.10 | 0% (0%/0%) err - | 50% (50%/50%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 |
| 0.15 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 |
| 0.20 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 |
| 0.30 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 |
| 0.50 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/0%) err 0.00 | 100% (100%/100%) err 0.00 |
| 1.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 |

## loh

Detection rate (type right / class right), and the median cell-fraction error, by cell fraction (rows) and size (columns, Mb; arm = a whole q arm; whole = a whole chromosome).

| f | 2 | 5 | 10 | 25 | arm | whole |
|---|---|---|---|---|---|---|
| 0.05 | 0% (0%/0%) err - | 0% (0%/0%) err - | 0% (0%/0%) err - | 50% (50%/50%) err 0.00 | 100% (100%/0%) err 0.00 | 100% (100%/50%) err 0.01 |
| 0.10 | 50% (50%/50%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.01 |
| 0.15 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 |
| 0.20 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 |
| 0.30 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 |
| 0.50 | 100% (100%/100%) err 0.01 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/100%) err 0.00 | 100% (100%/0%) err 0.00 | 100% (100%/100%) err 0.00 |

## False positives

0 events on 0 of 80 null samples (the brief asks for no event in 19 of 20): see false_positives.tsv.

