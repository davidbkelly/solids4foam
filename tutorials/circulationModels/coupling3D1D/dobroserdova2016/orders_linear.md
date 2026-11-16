# Time order of dobroserdova2016 (self-convergence, 0-0.3 s, from rest; 1D at CFL 0.1, p residual 1e-10, linearised vessels at fourth order (kExact 3, RK4, twice the cells))

Besides the self-convergence, every run is compared with the exact solution of the linear problem (RMS error / RMS of the signal), with the order between successive steps.

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 1.85 | 1.99 |
| 0.001 / 0.0005 / 0.00025 | 0.81 | 0.81 |
| 0.0005 / 0.00025 / 0.000125 | 0.95 | 0.93 |

Against the exact solution (Euler):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.36e-01 |  | 2.72e-01 |  |
| 0.001 | 8.55e-02 | 0.67 | 1.58e-01 | 0.78 |
| 0.0005 | 4.64e-02 | 0.88 | 8.67e-02 | 0.87 |
| 0.00025 | 2.38e-02 | 0.96 | 4.49e-02 | 0.95 |
| 0.000125 | 1.21e-02 | 0.98 | 2.29e-02 | 0.97 |

Finest run order_Euler_dt0.000125_linear: 2400 steps, 6.8 PIMPLE iterations per step, 143 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.20 | 1.84 |
| 0.001 / 0.0005 / 0.00025 | 2.18 | 1.99 |
| 0.0005 / 0.00025 / 0.000125 | 2.00 | 1.96 |

Against the exact solution (backward):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.06e-02 |  | 3.30e-02 |  |
| 0.001 | 2.42e-03 | 2.13 | 9.17e-03 | 1.85 |
| 0.0005 | 5.51e-04 | 2.14 | 2.36e-03 | 1.96 |
| 0.00025 | 1.38e-04 | 2.00 | 7.79e-04 | 1.60 |
| 0.000125 | 3.44e-05 | 2.00 | 5.37e-04 | 0.54 |

Finest run order_backward_dt0.000125_linear: 2400 steps, 6.8 PIMPLE iterations per step, 142 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.77 | 2.41 |
| 0.001 / 0.0005 / 0.00025 | 2.94 | 2.45 |
| 0.0005 / 0.00025 / 0.000125 | 2.81 | 2.51 |

Against the exact solution (BDF3):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.48e-03 |  | 1.37e-02 |  |
| 0.001 | 2.19e-04 | 2.76 | 2.61e-03 | 2.39 |
| 0.0005 | 2.91e-05 | 2.91 | 6.71e-04 | 1.96 |
| 0.00025 | 4.48e-06 | 2.70 | 5.20e-04 | 0.37 |
| 0.000125 | 2.21e-06 | 1.02 | 5.25e-04 | -0.01 |

Finest run order_BDF3_dt0.000125_linear: 2400 steps, 6.8 PIMPLE iterations per step, 142 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 3.76 | 2.62 |
| 0.001 / 0.0005 / 0.00025 | 3.01 | 2.54 |
| 0.0005 / 0.00025 / 0.000125 | 2.82 | 2.72 |

Against the exact solution (BDF4):

| dt | to3D.Q | order | to3D.P | order |
|---|---|---|---|---|
| 0.002 | 1.28e-03 |  | 1.12e-02 |  |
| 0.001 | 9.63e-05 | 3.74 | 1.89e-03 | 2.56 |
| 0.0005 | 1.22e-05 | 2.98 | 5.90e-04 | 1.68 |
| 0.00025 | 2.15e-06 | 2.50 | 5.30e-04 | 0.15 |
| 0.000125 | 2.26e-06 | -0.07 | 5.29e-04 | 0.00 |

Finest run order_BDF4_dt0.000125_linear: 2400 steps, 6.8 PIMPLE iterations per step, 141 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 1.21e-02, to3D.P: 2.29e-02
- backward: to3D.Q: 3.44e-05, to3D.P: 1.51e-04
- BDF3: to3D.Q: 5.60e-07, to3D.P: 1.21e-05
