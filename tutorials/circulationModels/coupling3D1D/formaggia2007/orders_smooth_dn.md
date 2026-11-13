# Time order of the coupled FSI + 1D tube (self-convergence, 0-12 ms, smooth sin^2 pulse, Dirichlet-Neumann IQN-ILS)

## backward (q = 2)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | -2.56 | 1.78 | 1.78 | 2.33 |

## BDF3s2 (q = 3)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | -0.29 | 1.30 | 1.28 | -0.07 |

## BDF4s2 (q = 4)

| dt pair | P(z=2.5 cm) | P(port) | Q(port) | wall u_r(z=2.5 cm) |
|---|---|---|---|---|
| 0.0001 / 5e-05 / 2.5e-05 | -0.15 | 1.07 | 1.05 | -0.14 |
