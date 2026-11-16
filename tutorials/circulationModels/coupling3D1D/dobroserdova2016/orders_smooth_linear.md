# Time order of dobroserdova2016 (self-convergence, 0-0.3 s, from rest; 1D at CFL 0.1, p residual 1e-10, linearised vessels at fourth order (kExact 3, RK4, twice the cells))

Besides the self-convergence, every run is compared with the exact solution of the linear problem (RMS error / RMS of the signal), with the order between successive steps.

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 0.99 | 1.00 |
| 0.001 / 0.0005 / 0.00025 | 0.90 | 0.89 |
| 0.0005 / 0.00025 / 0.000125 | 0.98 | 0.97 |

Against the exact solution (Euler):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.59e-01 |  | 1.35e-01 |  |
| 0.001 | 1.55e-01 | 0.04 | 1.44e-01 | -0.10 |
| 0.0005 | 8.14e-02 | 0.93 | 7.69e-02 | 0.91 |
| 0.00025 | 4.12e-02 | 0.98 | 3.90e-02 | 0.98 |
| 0.000125 | 2.07e-02 | 0.99 | 1.97e-02 | 0.99 |

Finest run order_Euler_dt0.000125_smooth_linear: 2400 steps, 6.7 PIMPLE iterations per step, 159 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.06 | 1.75 |
| 0.001 / 0.0005 / 0.00025 | 2.22 | 2.32 |
| 0.0005 / 0.00025 / 0.000125 | 2.01 | 2.01 |

Against the exact solution (backward):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.73e-02 |  | 1.75e-02 |  |
| 0.001 | 4.16e-03 | 2.06 | 4.88e-03 | 1.84 |
| 0.0005 | 9.26e-04 | 2.17 | 1.03e-03 | 2.24 |
| 0.00025 | 2.31e-04 | 2.01 | 2.56e-04 | 2.01 |
| 0.000125 | 5.77e-05 | 2.00 | 6.41e-05 | 2.00 |

Finest run order_backward_dt0.000125_smooth_linear: 2400 steps, 6.6 PIMPLE iterations per step, 153 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 3.03 | 3.03 |
| 0.001 / 0.0005 / 0.00025 | 3.16 | 3.20 |
| 0.0005 / 0.00025 / 0.000125 | 3.00 | 3.00 |

Against the exact solution (BDF3):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.59e-03 |  | 1.53e-03 |  |
| 0.001 | 1.94e-04 | 3.03 | 1.87e-04 | 3.03 |
| 0.0005 | 2.20e-05 | 3.14 | 2.07e-05 | 3.17 |
| 0.00025 | 2.74e-06 | 3.00 | 2.58e-06 | 3.00 |
| 0.000125 | 3.44e-07 | 3.00 | 3.23e-07 | 3.00 |

Finest run order_BDF3_dt0.000125_smooth_linear: 2400 steps, 6.5 PIMPLE iterations per step, 157 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 4.11 | 3.98 |
| 0.001 / 0.0005 / 0.00025 | 4.11 | 4.09 |
| 0.0005 / 0.00025 / 0.000125 | 4.00 | 4.00 |

Against the exact solution (BDF4):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.55e-04 |  | 2.08e-04 |  |
| 0.001 | 9.02e-06 | 4.11 | 1.31e-05 | 3.98 |
| 0.0005 | 5.25e-07 | 4.10 | 7.77e-07 | 4.08 |
| 0.00025 | 3.31e-08 | 3.99 | 4.92e-08 | 3.98 |
| 0.000125 | 2.53e-09 | 3.71 | 3.89e-09 | 3.66 |

Finest run order_BDF4_dt0.000125_smooth_linear: 2400 steps, 6.5 PIMPLE iterations per step, 152 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 2.07e-02, to3D.P: 1.97e-02
- backward: to3D.Q: 5.77e-05, to3D.P: 6.41e-05
- BDF3: to3D.Q: 3.43e-07, to3D.P: 3.22e-07
