# The BDF 4 floor of thoracicAorta37: one change at a time

BDF 4 at dt = 1, 0.5, 0.25 ms, 0-0.3 s from the smooth start; differences between successive steps, RMS over the window relative to the RMS of the finest run (to3D.Q, to3D.P), and the order they give; the last column is the finest run against the finest baseline run.

| variant | Q: e(1, 0.5) | e(0.5, 0.25) | order | P: e(1, 0.5) | e(0.5, 0.25) | order | vs baseline (Q) |
|---|---|---|---|---|---|---|---|
| baseline | 2.18e-05 | 1.88e-05 | 0.21 | 6.11e-06 | 5.48e-06 | 0.16 | 0.00e+00 |
| ddtCorrOff | 1.48e-05 | 1.04e-05 | 0.50 | 2.92e-06 | 1.80e-06 | 0.70 | 6.01e-05 |
| tightSolvers | 2.18e-05 | 1.88e-05 | 0.21 | 6.11e-06 | 5.48e-06 | 0.16 | 3.04e-13 |
| cylinder | 2.76e-06 | 1.95e-07 | 3.82 | 1.02e-06 | 6.82e-08 | 3.91 | 7.38e-03 |
| uniformStubs | 2.18e-05 | 1.89e-05 | 0.21 | 6.11e-06 | 5.48e-06 | 0.16 | 8.97e-04 |
| nonOrthCorr | 2.18e-05 | 1.88e-05 | 0.21 | 6.11e-06 | 5.48e-06 | 0.16 | 6.92e-13 |
