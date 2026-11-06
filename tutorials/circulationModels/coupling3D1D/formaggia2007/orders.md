# Time order of the coupled FSI + 1D tube (self-convergence, 0-12 ms)

## Euler (q = 1)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | 0.74 | 0.67 | 0.66 | 0.66 |
| 5e-05 / 2.5e-05 / 1.25e-05 | 0.82 | 0.47 | 0.47 | 0.76 |

## backward (q = 2)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | -0.29 | 0.94 | 0.96 | 1.55 |
| 5e-05 / 2.5e-05 / 1.25e-05 | 0.05 | 1.50 | 1.47 | 1.33 |

BDF3: 0 finished run(s), need 3

BDF4: 0 finished run(s), need 3
