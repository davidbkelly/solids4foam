# Time order of womersleySegment (self-convergence, 0-0.6 s, from rest; 1D at CFL 0.1, p residual 1e-10)

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 0.99 | 0.99 | 0.99 | 0.99 |
| 0.002 / 0.001 / 0.0005 | 1.00 | 1.00 | 1.00 | 1.00 |

Finest run order_Euler_dt0.0005: 1200 steps, 24.9 PIMPLE iterations per step, 665 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 1.97 | 1.97 | 1.97 | 1.97 |
| 0.002 / 0.001 / 0.0005 | 1.99 | 1.99 | 1.99 | 1.99 |

Finest run order_backward_dt0.0005: 1200 steps, 19.4 PIMPLE iterations per step, 469 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 2.79 | 2.79 | 2.79 | 2.79 |
| 0.002 / 0.001 / 0.0005 | 2.92 | 2.92 | 2.92 | 2.92 |

Finest run order_BDF3_dt0.0005: 1200 steps, 17.3 PIMPLE iterations per step, 478 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 3.31 | 3.31 | 3.31 | 3.31 |
| 0.002 / 0.001 / 0.0005 | 3.76 | 3.76 | 3.76 | 3.76 |

Finest run order_BDF4_dt0.0005: 1200 steps, 16.0 PIMPLE iterations per step, 448 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 1.36e-04, to3D.P: 1.29e-04, from3D.Q: 1.36e-04, from3D.P: 1.36e-04
- backward: to3D.Q: 4.80e-07, to3D.P: 4.55e-07, from3D.Q: 4.80e-07, from3D.P: 4.80e-07
- BDF3: to3D.Q: 1.41e-08, to3D.P: 1.34e-08, from3D.Q: 1.41e-08, from3D.P: 1.41e-08
