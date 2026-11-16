# Time order of thoracicAorta37 (self-convergence, 0-0.3 s, from rest; 1D at CFL 0.1, p residual 1e-10)

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 0.95 | 0.87 | 0.95 | 0.96 | 0.95 | 0.95 |
| 0.001 / 0.0005 / 0.00025 | 0.98 | 0.94 | 0.98 | 0.99 | 0.98 | 0.98 |

Finest run order_Euler_dt0.00025_cylinder: 1200 steps, 10.0 PIMPLE iterations per step, 461 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 1.99 | 2.01 | 1.99 | 1.98 | 2.00 | 1.99 |
| 0.001 / 0.0005 / 0.00025 | 2.00 | 2.01 | 2.00 | 2.00 | 2.01 | 2.00 |

Finest run order_backward_dt0.00025_cylinder: 1200 steps, 9.6 PIMPLE iterations per step, 454 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.99 | 2.99 | 2.99 | 2.98 | 2.99 | 2.99 |
| 0.001 / 0.0005 / 0.00025 | 3.02 | 3.02 | 3.02 | 3.01 | 3.02 | 3.02 |

Finest run order_BDF3_dt0.00025_cylinder: 1200 steps, 9.6 PIMPLE iterations per step, 423 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P | mid.P | mid.Q |
|---|---|---|---|---|---|---|
| 0.002 / 0.001 / 0.0005 | 3.83 | 3.72 | 3.83 | 3.69 | 3.78 | 3.83 |
| 0.001 / 0.0005 / 0.00025 | 3.82 | 3.91 | 3.82 | 3.82 | 3.91 | 3.82 |

Finest run order_BDF4_dt0.00025_cylinder: 1200 steps, 9.6 PIMPLE iterations per step, 852 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 2.38e-03, to3D.P: 3.15e-04, from3D.Q: 2.38e-03, from3D.P: 3.47e-04, mid.P: 2.67e-04, mid.Q: 2.38e-03
- backward: to3D.Q: 1.50e-05, to3D.P: 3.72e-06, from3D.Q: 1.50e-05, from3D.P: 2.81e-06, mid.P: 1.73e-06, mid.Q: 1.50e-05
- BDF3: to3D.Q: 3.40e-07, to3D.P: 9.73e-08, from3D.Q: 3.40e-07, from3D.P: 6.86e-08, mid.P: 4.15e-08, mid.Q: 3.40e-07
