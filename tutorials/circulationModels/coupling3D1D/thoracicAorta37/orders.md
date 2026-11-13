# Time order of thoracicAorta37 (self-convergence, 0-0.3 s, from rest; 1D at CFL 0.1, p residual 1e-10)

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 0.96 | 0.88 | 0.96 | 0.96 | 0.94 | 0.96 |
| 0.001 / 0.0005 / 0.00025 | 0.98 | 0.95 | 0.98 | 0.98 | 0.98 | 0.98 |

Finest run order_Euler_dt0.00025: 1200 steps, 11.1 PIMPLE iterations per step, 450 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.00 | 2.01 | 2.00 | 1.90 | 1.84 | 2.00 |
| 0.001 / 0.0005 / 0.00025 | 1.95 | 1.89 | 1.95 | 1.35 | 1.12 | 1.95 |

Finest run order_backward_dt0.00025: 1200 steps, 11.1 PIMPLE iterations per step, 450 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.34 | 2.38 | 2.34 | 1.41 | 1.06 | 2.34 |
| 0.001 / 0.0005 / 0.00025 | 0.61 | 0.54 | 0.61 | 0.22 | 0.22 | 0.61 |

Finest run order_BDF3_dt0.00025: 1200 steps, 11.1 PIMPLE iterations per step, 448 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 1.08 | 1.30 | 1.08 | 0.37 | 0.23 | 1.08 |
| 0.001 / 0.0005 / 0.00025 | 0.21 | 0.16 | 0.21 | 0.14 | 0.19 | 0.21 |

Finest run order_BDF4_dt0.00025: 1200 steps, 11.1 PIMPLE iterations per step, 449 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 2.28e-03, to3D.P: 3.07e-04, from3D.Q: 2.28e-03, from3D.P: 3.33e-04, mid.P: 2.60e-04, mid.Q: 2.28e-03
- backward: to3D.Q: 1.55e-05, to3D.P: 4.21e-06, from3D.Q: 1.55e-05, from3D.P: 5.47e-06, mid.P: 4.47e-06, mid.Q: 1.55e-05
- BDF3: to3D.Q: 3.34e-06, to3D.P: 9.75e-07, from3D.Q: 3.34e-06, from3D.P: 1.81e-06, mid.P: 1.54e-06, mid.Q: 3.34e-06
