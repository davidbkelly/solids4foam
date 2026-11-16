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
| 0.004 | 1.08e-03 |  | 1.03e-03 |  | 1.08e-03 |  | 1.08e-03 |  |
| 0.002 | 5.48e-04 | 0.99 | 5.19e-04 | 0.99 | 5.48e-04 | 0.99 | 5.48e-04 | 0.99 |
| 0.001 | 2.77e-04 | 0.99 | 2.62e-04 | 0.99 | 2.77e-04 | 0.99 | 2.77e-04 | 0.99 |
| 0.0005 | 1.41e-04 | 0.97 | 1.33e-04 | 0.97 | 1.41e-04 | 0.97 | 1.41e-04 | 0.97 |

Finest run order_Euler_dt0.0005_linear: 1200 steps, 27.5 PIMPLE iterations per step, 1696 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 1.97 | 1.97 | 1.97 | 1.97 |
| 0.002 / 0.001 / 0.0005 | 1.99 | 1.99 | 1.99 | 1.99 |

Against the exact solution (backward):

| dt | to3D.Q | order | to3D.P | order | from3D.Q | order | from3D.P | order |
|---|---|---|---|---|---|---|---|---|
| 0.004 | 3.09e-05 |  | 2.93e-05 |  | 3.09e-05 |  | 3.09e-05 |  |
| 0.002 | 9.86e-06 | 1.65 | 9.34e-06 | 1.65 | 9.86e-06 | 1.65 | 9.86e-06 | 1.65 |
| 0.001 | 6.06e-06 | 0.70 | 5.74e-06 | 0.70 | 6.06e-06 | 0.70 | 6.06e-06 | 0.70 |
| 0.0005 | 5.63e-06 | 0.10 | 5.34e-06 | 0.10 | 5.63e-06 | 0.10 | 5.63e-06 | 0.10 |

Finest run order_backward_dt0.0005_linear: 1200 steps, 21.5 PIMPLE iterations per step, 714 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 2.77 | 2.77 | 2.77 | 2.77 |
| 0.002 / 0.001 / 0.0005 | 2.86 | 2.86 | 2.86 | 2.86 |

Against the exact solution (BDF3):

| dt | to3D.Q | order | to3D.P | order | from3D.Q | order | from3D.P | order |
|---|---|---|---|---|---|---|---|---|
| 0.004 | 7.34e-06 |  | 6.96e-06 |  | 7.34e-06 |  | 7.34e-06 |  |
| 0.002 | 5.57e-06 | 0.40 | 5.28e-06 | 0.40 | 5.57e-06 | 0.40 | 5.57e-06 | 0.40 |
| 0.001 | 5.57e-06 | 0.00 | 5.28e-06 | 0.00 | 5.57e-06 | 0.00 | 5.57e-06 | 0.00 |
| 0.0005 | 5.57e-06 | -0.00 | 5.28e-06 | -0.00 | 5.57e-06 | -0.00 | 5.57e-06 | -0.00 |

Finest run order_BDF3_dt0.0005_linear: 1200 steps, 19.2 PIMPLE iterations per step, 1269 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P | from3D.Q | from3D.P |
|---|---|---|---|---|
| 0.004 / 0.002 / 0.001 | 3.20 | 3.20 | 3.20 | 3.20 |
| 0.002 / 0.001 / 0.0005 | 3.32 | 3.32 | 3.32 | 3.32 |

Against the exact solution (BDF4):

| dt | to3D.Q | order | to3D.P | order | from3D.Q | order | from3D.P | order |
|---|---|---|---|---|---|---|---|---|
| 0.004 | 6.55e-06 |  | 6.21e-06 |  | 6.55e-06 |  | 6.55e-06 |  |
| 0.002 | 5.61e-06 | 0.22 | 5.31e-06 | 0.22 | 5.61e-06 | 0.22 | 5.61e-06 | 0.22 |
| 0.001 | 5.58e-06 | 0.01 | 5.29e-06 | 0.01 | 5.58e-06 | 0.01 | 5.58e-06 | 0.01 |
| 0.0005 | 5.57e-06 | 0.00 | 5.28e-06 | 0.00 | 5.57e-06 | 0.00 | 5.57e-06 | 0.00 |

Finest run order_BDF4_dt0.0005_linear: 1200 steps, 17.7 PIMPLE iterations per step, 983 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 1.36e-04, to3D.P: 1.29e-04, from3D.Q: 1.36e-04, from3D.P: 1.36e-04
- backward: to3D.Q: 4.80e-07, to3D.P: 4.55e-07, from3D.Q: 4.80e-07, from3D.P: 4.80e-07
- BDF3: to3D.Q: 1.54e-08, to3D.P: 1.46e-08, from3D.Q: 1.54e-08, from3D.P: 1.54e-08
