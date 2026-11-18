# Time order of womersleySegment (self-convergence, 0-0.6 s, from rest; 1D at CFL 0.1, p residual 1e-10, linearised vessels at fourth order (kExact 3, RK4, twice the cells))

Besides the self-convergence, every run is compared with the exact solution of the linear problem (RMS error / RMS of the signal), with the order between successive steps.

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 0.99 | 0.99 | 0.99 | 0.99 |
| 0.002 / 0.001 / 0.0005 | 1.00 | 1.00 | 1.00 | 1.00 |

Against the exact solution (Euler):

| dt | to3D.Q | order | to3D.P | order | from3D.Q | order | from3D.P | order |
|---|---|---|---|---|---|---|---|---|
| 0.004 | 1.08e-03 |  | 1.02e-03 |  | 1.08e-03 |  | 1.08e-03 |  |
| 0.002 | 5.42e-04 | 0.99 | 5.14e-04 | 0.99 | 5.42e-04 | 0.99 | 5.42e-04 | 0.99 |
| 0.001 | 2.72e-04 | 1.00 | 2.57e-04 | 1.00 | 2.72e-04 | 1.00 | 2.72e-04 | 1.00 |
| 0.0005 | 1.36e-04 | 1.00 | 1.29e-04 | 1.00 | 1.36e-04 | 1.00 | 1.36e-04 | 1.00 |

Finest run order_Euler_dt0.0005_linear: 1200 steps, 28.1 PIMPLE iterations per step, 1577 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 1.97 | 1.97 | 1.97 | 1.97 |
| 0.002 / 0.001 / 0.0005 | 1.99 | 1.99 | 1.99 | 1.99 |

Against the exact solution (backward):

| dt | to3D.Q | order | to3D.P | order | from3D.Q | order | from3D.P | order |
|---|---|---|---|---|---|---|---|---|
| 0.004 | 2.62e-05 |  | 2.48e-05 |  | 2.62e-05 |  | 2.62e-05 |  |
| 0.002 | 6.60e-06 | 1.99 | 6.26e-06 | 1.99 | 6.60e-06 | 1.99 | 6.60e-06 | 1.99 |
| 0.001 | 6.15e-06 | 0.10 | 5.83e-06 | 0.10 | 6.15e-06 | 0.10 | 6.15e-06 | 0.10 |
| 0.0005 | 6.83e-06 | -0.15 | 6.47e-06 | -0.15 | 6.83e-06 | -0.15 | 6.83e-06 | -0.15 |

Finest run order_backward_dt0.0005_linear: 1200 steps, 21.9 PIMPLE iterations per step, 1241 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 2.77 | 2.77 | 2.77 | 2.77 |
| 0.002 / 0.001 / 0.0005 | 2.86 | 2.86 | 2.86 | 2.86 |

Against the exact solution (BDF3):

| dt | to3D.Q | order | to3D.P | order | from3D.Q | order | from3D.P | order |
|---|---|---|---|---|---|---|---|---|
| 0.004 | 8.78e-06 |  | 8.33e-06 |  | 8.78e-06 |  | 8.78e-06 |  |
| 0.002 | 7.11e-06 | 0.30 | 6.74e-06 | 0.30 | 7.11e-06 | 0.30 | 7.11e-06 | 0.30 |
| 0.001 | 7.10e-06 | 0.00 | 6.73e-06 | 0.00 | 7.10e-06 | 0.00 | 7.10e-06 | 0.00 |
| 0.0005 | 7.10e-06 | -0.00 | 6.73e-06 | -0.00 | 7.10e-06 | -0.00 | 7.10e-06 | -0.00 |

Finest run order_BDF3_dt0.0005_linear: 1200 steps, 19.3 PIMPLE iterations per step, 1031 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 3.20 | 3.20 | 3.20 | 3.20 |
| 0.002 / 0.001 / 0.0005 | 3.32 | 3.32 | 3.32 | 3.32 |

Against the exact solution (BDF4):

| dt | to3D.Q | order | to3D.P | order | from3D.Q | order | from3D.P | order |
|---|---|---|---|---|---|---|---|---|
| 0.004 | 7.95e-06 |  | 7.54e-06 |  | 7.95e-06 |  | 7.95e-06 |  |
| 0.002 | 7.12e-06 | 0.16 | 6.75e-06 | 0.16 | 7.12e-06 | 0.16 | 7.12e-06 | 0.16 |
| 0.001 | 7.10e-06 | 0.00 | 6.73e-06 | 0.00 | 7.10e-06 | 0.00 | 7.10e-06 | 0.00 |
| 0.0005 | 7.10e-06 | 0.00 | 6.73e-06 | 0.00 | 7.10e-06 | 0.00 | 7.10e-06 | 0.00 |

Finest run order_BDF4_dt0.0005_linear: 1200 steps, 17.9 PIMPLE iterations per step, 956 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 1.35e-04, to3D.P: 1.28e-04, from3D.Q: 1.35e-04, from3D.P: 1.35e-04
- backward: to3D.Q: 4.80e-07, to3D.P: 4.55e-07, from3D.Q: 4.80e-07, from3D.P: 4.80e-07
- BDF3: to3D.Q: 1.54e-08, to3D.P: 1.46e-08, from3D.Q: 1.54e-08, from3D.P: 1.54e-08
