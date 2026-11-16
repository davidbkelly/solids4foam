# Rigid segment 3D-1D at fourth order: against the closed form

Vessels: linearised 1D equations, kExact (degree 3) + RK4, CFL 0.5; pipe: pimpleFoam, BDF 4; full amplitude. dx and dt refined together. Errors against the closed form (tau = 6.5226 ms with the 32-sided patch area), relative to the peaks: 3D flux, pressures at the pipe ends (max of the two), profiles of Q along both vessels at t = 0.35 s (cell averages of the closed form).

## implicit order 4

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.43e-03 |  | 1.10e-03 |  | 1.75e-03 |  |
| 800 | 1 | 9.41e-05 | 3.93 | 7.30e-05 | 3.92 | 1.14e-04 | 3.95 |
| 1600 | 0.5 | 5.99e-06 | 3.97 | 4.62e-06 | 3.98 | 7.14e-06 | 3.99 |
| 3200 | 0.25 | 3.77e-07 | 3.99 | 2.89e-07 | 4.00 | 4.46e-07 | 4.00 |
| 6400 | 0.125 | 2.36e-08 | 4.00 | 1.80e-08 | 4.00 | 2.82e-08 | 3.98 |

## explicit order 4

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.43e-03 |  | 5.62e-03 |  | 1.97e-03 |  |
| 800 | 1 | 9.41e-05 | 3.93 | 3.55e-04 | 3.99 | 1.23e-04 | 4.00 |
| 1600 | 0.5 | 5.99e-06 | 3.97 | 2.22e-05 | 4.00 | 7.66e-06 | 4.01 |
| 3200 | 0.25 | 3.77e-07 | 3.99 | 1.39e-06 | 4.00 | 4.76e-07 | 4.01 |
| 6400 | 0.125 | 2.36e-08 | 4.00 | 8.67e-08 | 4.00 | 2.97e-08 | 4.00 |

## implicit order 2

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.43e-03 |  | 1.10e-03 |  | 2.69e-03 |  |
| 800 | 1 | 9.41e-05 | 3.93 | 7.30e-05 | 3.92 | 8.78e-04 | 1.61 |
| 1600 | 0.5 | 5.99e-06 | 3.97 | 4.62e-06 | 3.98 | 2.38e-04 | 1.89 |
| 3200 | 0.25 | 3.77e-07 | 3.99 | 2.89e-07 | 4.00 | 6.05e-05 | 1.97 |
| 6400 | 0.125 | 2.36e-08 | 4.00 | 1.80e-08 | 4.00 | 1.52e-05 | 1.99 |

## implicit order 1

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.43e-03 |  | 1.10e-03 |  | 6.14e-02 |  |
| 800 | 1 | 9.41e-05 | 3.93 | 7.30e-05 | 3.92 | 3.11e-02 | 0.98 |
| 1600 | 0.5 | 5.99e-06 | 3.97 | 4.62e-06 | 3.98 | 1.56e-02 | 1.00 |
| 3200 | 0.25 | 3.77e-07 | 3.99 | 2.89e-07 | 4.00 | 7.79e-03 | 1.00 |
| 6400 | 0.125 | 2.36e-08 | 4.00 | 1.80e-08 | 4.00 | 3.90e-03 | 1.00 |
