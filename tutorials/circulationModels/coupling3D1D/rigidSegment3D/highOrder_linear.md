# Rigid segment 3D-1D at fourth order: against the closed form

Vessels: linearised 1D equations, kExact (degree 3) + RK4, CFL 0.5; pipe: pimpleFoam, BDF 4; full amplitude. dx and dt refined together. Errors against the closed form (tau = 6.4807 ms with the patch area pi R^2), relative to the peaks: 3D flux, pressures at the pipe ends (max of the two), profiles of Q along both vessels at t = 0.35 s (cell averages of the closed form).

## implicit order 4

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.44e-03 |  | 1.11e-03 |  | 1.76e-03 |  |
| 800 | 1 | 9.43e-05 | 3.93 | 7.33e-05 | 3.92 | 1.14e-04 | 3.95 |
| 1600 | 0.5 | 6.01e-06 | 3.97 | 4.64e-06 | 3.98 | 7.16e-06 | 3.99 |
| 3200 | 0.25 | 3.78e-07 | 3.99 | 2.90e-07 | 4.00 | 4.48e-07 | 4.00 |
| 6400 | 0.125 | 2.36e-08 | 4.00 | 1.81e-08 | 4.00 | 2.84e-08 | 3.98 |

## explicit order 4

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 4.15e-03 |  | 4.18e-03 |  | 2.99e-03 |  |
| 800 | 1 | 2.81e-04 | 3.88 | 2.81e-04 | 3.89 | 1.92e-04 | 3.96 |
| 1600 | 0.5 | 1.80e-05 | 3.97 | 1.77e-05 | 3.99 | 1.20e-05 | 4.00 |
| 3200 | 0.25 | 1.13e-06 | 3.99 | 1.11e-06 | 4.00 | 7.49e-07 | 4.00 |
| 6400 | 0.125 | 7.11e-08 | 3.99 | 6.94e-08 | 4.00 | 4.69e-08 | 4.00 |

## implicit order 2

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.44e-03 |  | 1.11e-03 |  | 2.70e-03 |  |
| 800 | 1 | 9.43e-05 | 3.93 | 7.33e-05 | 3.92 | 8.82e-04 | 1.61 |
| 1600 | 0.5 | 6.01e-06 | 3.97 | 4.64e-06 | 3.98 | 2.39e-04 | 1.89 |
| 3200 | 0.25 | 3.78e-07 | 3.99 | 2.90e-07 | 4.00 | 6.07e-05 | 1.97 |
| 6400 | 0.125 | 2.36e-08 | 4.00 | 1.81e-08 | 4.00 | 1.52e-05 | 1.99 |

## implicit order 1

| cells | dt [ms] | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| 400 | 2 | 1.44e-03 |  | 1.11e-03 |  | 6.16e-02 |  |
| 800 | 1 | 9.43e-05 | 3.93 | 7.33e-05 | 3.92 | 3.12e-02 | 0.98 |
| 1600 | 0.5 | 6.01e-06 | 3.97 | 4.64e-06 | 3.98 | 1.56e-02 | 1.00 |
| 3200 | 0.25 | 3.78e-07 | 3.99 | 2.90e-07 | 4.00 | 7.81e-03 | 1.00 |
| 6400 | 0.125 | 2.36e-08 | 4.00 | 1.81e-08 | 4.00 | 3.91e-03 | 1.00 |
