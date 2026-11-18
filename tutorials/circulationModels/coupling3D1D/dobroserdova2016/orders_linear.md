# Time order of dobroserdova2016 (self-convergence, 0-0.3 s, from rest; 1D at CFL 0.1, p residual 1e-10, linearised vessels at fourth order (kExact 3, RK4, twice the cells))

Besides the self-convergence, every run is compared with the exact solution of the linear problem (RMS error / RMS of the signal), with the order between successive steps.

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 1.85 | 2.00 |
| 0.001 / 0.0005 / 0.00025 | 0.81 | 0.81 |
| 0.0005 / 0.00025 / 0.000125 | 0.95 | 0.93 |

Against the exact solution (Euler):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.36e-01 |  | 2.72e-01 |  |
| 0.001 | 8.58e-02 | 0.66 | 1.59e-01 | 0.78 |
| 0.0005 | 4.67e-02 | 0.88 | 8.71e-02 | 0.87 |
| 0.00025 | 2.40e-02 | 0.96 | 4.50e-02 | 0.95 |
| 0.000125 | 1.22e-02 | 0.98 | 2.30e-02 | 0.97 |

Finest run order_Euler_dt0.000125_linear: 2400 steps, 7.0 PIMPLE iterations per step, 210 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.19 | 1.84 |
| 0.001 / 0.0005 / 0.00025 | 2.18 | 1.99 |
| 0.0005 / 0.00025 / 0.000125 | 2.00 | 1.96 |

Against the exact solution (backward):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.06e-02 |  | 3.33e-02 |  |
| 0.001 | 2.44e-03 | 2.12 | 9.28e-03 | 1.84 |
| 0.0005 | 5.54e-04 | 2.14 | 2.39e-03 | 1.95 |
| 0.00025 | 1.38e-04 | 2.00 | 7.87e-04 | 1.60 |
| 0.000125 | 3.46e-05 | 2.00 | 5.41e-04 | 0.54 |

Finest run order_backward_dt0.000125_linear: 2400 steps, 6.9 PIMPLE iterations per step, 208 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.77 | 2.41 |
| 0.001 / 0.0005 / 0.00025 | 2.94 | 2.45 |
| 0.0005 / 0.00025 / 0.000125 | 2.81 | 2.51 |

Against the exact solution (BDF3):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.50e-03 |  | 1.39e-02 |  |
| 0.001 | 2.21e-04 | 2.76 | 2.65e-03 | 2.39 |
| 0.0005 | 2.94e-05 | 2.91 | 6.78e-04 | 1.96 |
| 0.00025 | 4.53e-06 | 2.70 | 5.25e-04 | 0.37 |
| 0.000125 | 2.23e-06 | 1.02 | 5.29e-04 | -0.01 |

Finest run order_BDF3_dt0.000125_linear: 2400 steps, 6.9 PIMPLE iterations per step, 207 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 3.76 | 2.62 |
| 0.001 / 0.0005 / 0.00025 | 3.01 | 2.54 |
| 0.0005 / 0.00025 / 0.000125 | 2.82 | 2.72 |

Against the exact solution (BDF4):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.30e-03 |  | 1.13e-02 |  |
| 0.001 | 9.73e-05 | 3.74 | 1.92e-03 | 2.56 |
| 0.0005 | 1.23e-05 | 2.99 | 5.96e-04 | 1.69 |
| 0.00025 | 2.17e-06 | 2.50 | 5.35e-04 | 0.15 |
| 0.000125 | 2.28e-06 | -0.07 | 5.34e-04 | 0.00 |

Finest run order_BDF4_dt0.000125_linear: 2400 steps, 6.9 PIMPLE iterations per step, 206 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 1.22e-02, to3D.P: 2.30e-02
- backward: to3D.Q: 3.45e-05, to3D.P: 1.53e-04
- BDF3: to3D.Q: 5.65e-07, to3D.P: 1.23e-05
