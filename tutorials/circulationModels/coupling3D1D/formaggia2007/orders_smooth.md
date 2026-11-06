# Time order of the coupled FSI + 1D tube (self-convergence, 0-12 ms, smooth sin^2 pulse)

## Euler (q = 1)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | 0.91 | 1.03 | 1.03 | 0.50 |
| 5e-05 / 2.5e-05 / 1.25e-05 | -0.58 | 0.74 | 0.75 | 0.71 |

## backward (q = 2)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | 3.63 | 1.40 | 1.41 | 2.10 |
| 5e-05 / 2.5e-05 / 1.25e-05 | 2.06 | 2.06 | 2.05 | 2.44 |

## BDF3s2 (q = 3)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | 3.20 | 4.54 | 4.54 | 3.89 |

BDF4s2: 0 finished run(s), need 3
