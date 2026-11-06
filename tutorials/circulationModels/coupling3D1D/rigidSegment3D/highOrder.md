# Rigid segment 3D-1D at fourth order: self-convergence

Vessels: kExact (degree 3) + RK4, CFL 0.5; pipe: pimpleFoam, BDF 4; full amplitude. dx and dt refined together; reference: implicit 4 at 6400 cells, dt = 0.125 ms. Errors relative to the peaks: 3D flux, pressures at the pipe ends (max of the two), profiles of Q along both vessels at t = 0.35 s.

## implicit order 4

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.76e-03 |  | 1.09e-03 |  | 2.17e-03 |  |
| 800 | 1 | 1.18e-04 | 3.90 | 7.31e-05 | 3.89 | 1.40e-04 | 3.95 |
| 1600 | 0.5 | 7.46e-06 | 3.98 | 4.59e-06 | 3.99 | 8.78e-06 | 4.00 |
| 3200 | 0.25 | 4.41e-07 | 4.08 | 2.70e-07 | 4.09 | 5.16e-07 | 4.09 |

## explicit order 4

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 5.07e-03 |  | 4.17e-03 |  | 3.62e-03 |  |
| 800 | 1 | 3.41e-04 | 3.90 | 2.79e-04 | 3.90 | 2.32e-04 | 3.96 |
| 1600 | 0.5 | 2.18e-05 | 3.97 | 1.76e-05 | 3.98 | 1.45e-05 | 4.00 |
| 3200 | 0.25 | 1.37e-06 | 3.99 | 1.12e-06 | 3.97 | 9.14e-07 | 3.98 |

## implicit order 2

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.76e-03 |  | 1.09e-03 |  | 3.30e-03 |  |
| 800 | 1 | 1.18e-04 | 3.90 | 7.31e-05 | 3.89 | 1.09e-03 | 1.60 |
| 1600 | 0.5 | 7.47e-06 | 3.98 | 4.61e-06 | 3.99 | 2.95e-04 | 1.89 |
| 3200 | 0.25 | 4.43e-07 | 4.08 | 2.73e-07 | 4.08 | 7.51e-05 | 1.97 |

## implicit order 1

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.76e-03 |  | 1.10e-03 |  | 7.73e-02 |  |
| 800 | 1 | 1.18e-04 | 3.89 | 7.94e-05 | 3.79 | 3.92e-02 | 0.98 |
| 1600 | 0.5 | 7.88e-06 | 3.91 | 7.86e-06 | 3.34 | 1.96e-02 | 1.00 |
| 3200 | 0.25 | 1.46e-06 | 2.43 | 2.16e-06 | 1.87 | 9.82e-03 | 1.00 |
