# Calibration

The detection curve of `arraykaryo events` on the mock at the example array's probe spacing (1.3 kb), measured by
`scripts/calibrate.py`: cohorts of 160 samples on eight chromosomes (chr1, 2, 3, 7, 11, 12, 15, 21; 1.1 Gb) with one planted
event per sample over a grid of kinds (gain, loss, copy-neutral LOH), cell fractions (0.05 to 1) and sizes (2, 5, 10, 25 Mb
stretches, a whole q arm, a whole chromosome), per-sample LRR SD 0.15 and BAF SD 0.03-0.05, plus 40 null samples per
replicate; two replicates. No real data: everything is simulated.

- `calibration.md`: the detection rate (type right / class right) and the median cell-fraction error per cell of the grid.
- `calibration.tsv`: one row per planted event (replicate, kind, f, size, chromosome, probes, detected, type and class called,
  the three cell-fraction estimates, the errors, other calls on the sample).
- `false_positives.tsv`: every event called on a null sample (none in the two replicates: the brief asks for no event in 19
  of 20).
- `calibration.png`: the grid as a heat map per kind.

Summary (0.5.2): a whole chromosome or arm is found from a cell fraction of 0.10 for a gain or a loss and 0.05 for a
copy-neutral LOH; a stretch of 5 Mb or more from 0.15 for a gain (0.10 for a loss or LOH at 10 Mb or more); a 2-Mb stretch
from 0.20 for a gain and 0.15 for a loss or LOH; the cell fraction is within 0.02 of the planted one. Below those fractions
nothing is called, and nothing is called on a null sample.

To repeat: `python scripts/calibrate.py --out docs/calibration --work /tmp/calib --threads 8 --replicates 2` (about 20 minutes).
