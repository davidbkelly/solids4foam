# Time order of dobroserdova2016 (self-convergence, 0-0.3 s, from rest; 1D at CFL 0.1, p residual 1e-10, linearised vessels at fourth order (kExact 3, RK4, twice the cells))

Besides the self-convergence, every run is compared with the exact solution of the linear problem (RMS error / RMS of the signal), with the order between successive steps.

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 0.98 | 1.01 |
| 0.001 / 0.0005 / 0.00025 | 0.89 | 0.89 |
| 0.0005 / 0.00025 / 0.000125 | 0.98 | 0.97 |

Against the exact solution (Euler):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.58e-01 |  | 1.34e-01 |  |
| 0.001 | 1.56e-01 | 0.02 | 1.45e-01 | -0.11 |
| 0.0005 | 8.19e-02 | 0.93 | 7.70e-02 | 0.91 |
| 0.00025 | 4.14e-02 | 0.98 | 3.90e-02 | 0.98 |
| 0.000125 | 2.09e-02 | 0.99 | 1.97e-02 | 0.99 |

Finest run order_Euler_dt0.000125_smooth_linear: 2400 steps, 7.2 PIMPLE iterations per step, 213 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.06 | 1.75 |
| 0.001 / 0.0005 / 0.00025 | 2.22 | 2.32 |
| 0.0005 / 0.00025 / 0.000125 | 2.01 | 2.01 |

Against the exact solution (backward):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.73e-02 |  | 1.76e-02 |  |
| 0.001 | 4.18e-03 | 2.05 | 4.92e-03 | 1.84 |
| 0.0005 | 9.29e-04 | 2.17 | 1.04e-03 | 2.25 |
| 0.00025 | 2.31e-04 | 2.01 | 2.58e-04 | 2.01 |
| 0.000125 | 5.78e-05 | 2.00 | 6.45e-05 | 2.00 |

Finest run order_backward_dt0.000125_smooth_linear: 2400 steps, 7.1 PIMPLE iterations per step, 216 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 3.03 | 3.02 |
| 0.001 / 0.0005 / 0.00025 | 3.17 | 3.20 |
| 0.0005 / 0.00025 / 0.000125 | 3.00 | 3.00 |

Against the exact solution (BDF3):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.61e-03 |  | 1.53e-03 |  |
| 0.001 | 1.96e-04 | 3.03 | 1.88e-04 | 3.02 |
| 0.0005 | 2.22e-05 | 3.14 | 2.08e-05 | 3.18 |
| 0.00025 | 2.77e-06 | 3.00 | 2.59e-06 | 3.00 |
| 0.000125 | 3.47e-07 | 3.00 | 3.25e-07 | 3.00 |

Finest run order_BDF3_dt0.000125_smooth_linear: 2400 steps, 7.0 PIMPLE iterations per step, 206 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 4.11 | 3.98 |
| 0.001 / 0.0005 / 0.00025 | 4.11 | 4.09 |
| 0.0005 / 0.00025 / 0.000125 | 4.00 | 4.00 |

Against the exact solution (BDF4):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.56e-04 |  | 2.09e-04 |  |
| 0.001 | 9.09e-06 | 4.10 | 1.32e-05 | 3.98 |
| 0.0005 | 5.28e-07 | 4.10 | 7.82e-07 | 4.08 |
| 0.00025 | 3.33e-08 | 3.99 | 4.95e-08 | 3.98 |
| 0.000125 | 2.53e-09 | 3.72 | 3.90e-09 | 3.67 |

Finest run order_BDF4_dt0.000125_smooth_linear: 2400 steps, 7.0 PIMPLE iterations per step, 209 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 2.09e-02, to3D.P: 1.97e-02
- backward: to3D.Q: 5.78e-05, to3D.P: 6.45e-05
- BDF3: to3D.Q: 3.46e-07, to3D.P: 3.24e-07
