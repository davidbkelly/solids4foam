# Time order of dobroserdova2016 (self-convergence, 0-0.3 s, from rest; 1D at CFL 0.1, p residual 1e-10)

## Euler (q = 1)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 0.41 | 0.41 |
| 0.001 / 0.0005 / 0.00025 | 0.87 | 0.85 |
| 0.0005 / 0.00025 / 0.000125 | 0.93 | 0.92 |

Finest run order_Euler_dt0.000125: 2400 steps, 7.0 PIMPLE iterations per step, 191 s

## backward (q = 2)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 1.78 | 1.71 |
| 0.001 / 0.0005 / 0.00025 | 1.99 | 1.92 |
| 0.0005 / 0.00025 / 0.000125 | 2.00 | 1.97 |

Finest run order_backward_dt0.000125: 2400 steps, 6.9 PIMPLE iterations per step, 189 s

## BDF3 (q = 3)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 2.87 | 2.45 |
| 0.001 / 0.0005 / 0.00025 | 2.91 | 2.73 |
| 0.0005 / 0.00025 / 0.000125 | 2.97 | 2.94 |

Finest run order_BDF3_dt0.000125: 2400 steps, 6.9 PIMPLE iterations per step, 183 s

## BDF4 (q = 4)

| dt pair | to3D.Q | to3D.P |
|---|---|---|
| 0.002 / 0.001 / 0.0005 | 3.33 | 2.82 |
| 0.001 / 0.0005 / 0.00025 | 3.48 | 3.43 |
| 0.0005 / 0.00025 / 0.000125 | 3.82 | 3.92 |

Finest run order_BDF4_dt0.000125: 2400 steps, 6.9 PIMPLE iterations per step, 185 s

## Finest runs against the finest BDF4 run (RMS of the difference / RMS of the signal)

- Euler: to3D.Q: 1.22e-02, to3D.P: 2.29e-02
- backward: to3D.Q: 3.45e-05, to3D.P: 1.46e-04
- BDF3: to3D.Q: 4.40e-07, to3D.P: 6.07e-06
