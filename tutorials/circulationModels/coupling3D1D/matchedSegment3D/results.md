# B5, the matched segment: refinement against the exact solution

3D cells along the pipe and the 3D step refined together (acoustic Courant 0.49); vessels at 1600 cells, k-exact 3 + RK4, linearised; implicit coupling of order 2 with 20 outer correctors per step, BDF2; errors relative to the incident peak (maximum over the run), the delay relative to L/c0 = 32.404 ms.

## Implicit coupling, Dirichlet interface

| dz [mm] | dt [ms] | P1 | order | Q1 | order | P2 | order | Q2 | order | refl | order | trans | order | delay | exch./step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.4 | 4.16e-03 |  | 4.16e-03 |  | 1.10e-02 |  | 1.10e-02 |  | 3.89e-03 |  | 1.12e-02 |  | +5.3e-03 | 20.0 |
| 2.5 | 0.2 | 1.03e-03 | 2.02 | 1.03e-03 | 2.02 | 2.66e-03 | 2.05 | 2.66e-03 | 2.05 | 9.59e-04 | 2.02 | 2.70e-03 | 2.05 | +1.3e-03 | 20.0 |
| 1.25 | 0.1 | 2.56e-04 | 2.00 | 2.55e-04 | 2.00 | 6.51e-04 | 2.03 | 6.51e-04 | 2.03 | 2.38e-04 | 2.01 | 6.70e-04 | 2.01 | +3.1e-04 | 20.0 |

## Explicit coupling of order 1

| dz [mm] | dt [ms] | P1 | order | Q1 | order | P2 | order | Q2 | order | refl | order | trans | order | delay | exch./step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.4 | 1.78e-02 |  | 1.77e-02 |  | 1.04e-02 |  | 2.63e-02 |  | 1.66e-02 |  | 1.16e-02 |  | +1.1e-02 | 1.0 |
