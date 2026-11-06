# Rigid segment, 3D-1D: pimpleFoam pipe coupled to 1D vessels

c0 = 6.1721 m/s, Zc = 2.0629e+07 Pa s/m3; pipe 8 cm, A3D = 3.121445e-04 m2 (0.9936 pi R^2: the 32-sided patch), tau = 6.5226 ms; incident peak 20.6288 mPa (inflow x 1e-3); 1D dx = 1.25 mm, SSPRK3 (CFL 0.5 with backward, 0.1 with BDF 3 and BDF 4).

## 1. Against the closed form (max error / incident peak)

| run | Q (3D flux) | P1, pipe inlet | P2, pipe outlet |
|---|---|---|---|
| monolithic 1D, CFL 0.5 | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| monolithic 1D, CFL 0.1 | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| backward: explicit 2, dt 4 ms | 9.55e-02 | 7.93e-02 | 7.51e-02 |
| backward: explicit 2, dt 2 ms | 2.88e-02 | 2.46e-02 | 2.04e-02 |
| backward: explicit 2, dt 1 ms | 8.31e-03 | 9.30e-03 | 5.93e-03 |
| backward: explicit 2, dt 0.5 ms | 2.80e-03 | 5.31e-03 | 2.23e-03 |
| backward: explicit 2, dt 0.25 ms | 1.42e-03 | 4.35e-03 | 1.45e-03 |
| backward: implicit 2, dt 4 ms | 3.59e-02 | 3.18e-02 | 3.59e-02 |
| backward: implicit 2, dt 2 ms | 1.03e-02 | 6.86e-03 | 1.03e-02 |
| backward: implicit 2, dt 1 ms | 2.64e-03 | 3.46e-03 | 2.64e-03 |
| backward: implicit 2, dt 0.5 ms | 1.08e-03 | 3.51e-03 | 1.08e-03 |
| backward: implicit 2, dt 0.25 ms | 1.12e-03 | 3.90e-03 | 1.12e-03 |
| backward: implicit 2, dt 0.03125 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| BDF3: explicit 3, dt 4 ms | 5.62e-02 | 5.23e-02 | 5.23e-02 |
| BDF3: explicit 3, dt 2 ms | 9.52e-03 | 9.92e-03 | 1.05e-02 |
| BDF3: explicit 3, dt 1 ms | 1.28e-03 | 3.91e-03 | 2.32e-03 |
| BDF3: explicit 3, dt 0.5 ms | 1.10e-03 | 3.97e-03 | 1.36e-03 |
| BDF3: explicit 3, dt 0.25 ms | 1.20e-03 | 4.04e-03 | 1.24e-03 |
| BDF3: implicit 3, dt 4 ms | 2.04e-02 | 2.04e-02 | 2.04e-02 |
| BDF3: implicit 3, dt 2 ms | 2.46e-03 | 5.18e-03 | 2.46e-03 |
| BDF3: implicit 3, dt 1 ms | 8.50e-04 | 4.19e-03 | 8.50e-04 |
| BDF3: implicit 3, dt 0.5 ms | 1.17e-03 | 4.07e-03 | 1.17e-03 |
| BDF3: implicit 3, dt 0.25 ms | 1.21e-03 | 4.05e-03 | 1.21e-03 |
| BDF3: implicit 3, dt 0.03125 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| BDF4: explicit 4, dt 4 ms | 3.91e-02 | 4.79e-02 | 4.80e-02 |
| BDF4: explicit 4, dt 2 ms | 3.14e-03 | 2.84e-03 | 6.26e-03 |
| BDF4: explicit 4, dt 1 ms | 9.92e-04 | 3.67e-03 | 1.31e-03 |
| BDF4: explicit 4, dt 0.5 ms | 1.20e-03 | 4.02e-03 | 1.22e-03 |
| BDF4: explicit 4, dt 0.25 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 4 ms | 1.62e-02 | 1.91e-02 | 1.62e-02 |
| BDF4: implicit 4, dt 2 ms | 1.44e-03 | 5.41e-03 | 1.44e-03 |
| BDF4: implicit 4, dt 1 ms | 1.22e-03 | 4.14e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 0.5 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 0.25 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 0.03125 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |

The monolithic rows are the 1D discretisation error alone (dx = 1.25 mm); a coupled run cannot do better than them. Only the coupling of the time scheme's own order is listed.

Mass through the rigid pipe (implicit 2, 0.5 ms): max |Qin - Qout| = 1.0e-12 of the peak; outer iterations 11.1 per step on average.

## 2. Convergence in dt, CFD backward (max |Q - Q_fine| / peak)

Q_fine: implicit order 2 at dt = 0.03125 ms.

| dt [ms] | dt/tau | explicit 1 | explicit 2 | explicit 3 | explicit 4 | implicit 2 |
|---|---|---|---|---|---|---|
| 4 | 0.613 | 2.23e-01 | 1.18e-01 | 9.81e-02 | unstable | 4.46e-02 |
| 2 | 0.307 | 1.14e-01 | 3.47e-02 | 2.35e-02 | 1.59e-02 | 1.27e-02 |
| 1 | 0.153 | 5.79e-02 | 9.30e-03 | 4.93e-03 | 3.55e-03 | 3.43e-03 |
| 0.5 | 0.077 | 2.91e-02 | 2.38e-03 | 1.07e-03 | 8.74e-04 | 8.71e-04 |
| 0.25 | 0.038 | 1.46e-02 | 6.02e-04 | 2.43e-04 | 2.18e-04 | 2.18e-04 |

Order between the two finest dt, and least squares over the three finest:

- explicit 1: 1.00 (0.99)
- explicit 2: 1.99 (1.97)
- explicit 3: 2.14 (2.17)
- explicit 4: 2.01 (2.01)
- implicit 2: 2.00 (1.99)

The fine coupled run against the monolithic 1D run (CFL 0.5): 3.42e-06 of the peak.

## 2. Convergence in dt, CFD BDF 3 (max |Q - Q_fine| / peak)

Q_fine: implicit order 3 at dt = 0.03125 ms.

| dt [ms] | dt/tau | explicit 1 | explicit 2 | explicit 3 | explicit 4 | implicit 1 | implicit 3 |
|---|---|---|---|---|---|---|---|
| 4 | 0.613 | 2.81e-01 | 1.01e-01 | 7.00e-02 | 6.42e-02 | 2.65e-02 | 2.65e-02 |
| 2 | 0.307 | 1.26e-01 | 3.03e-02 | 1.19e-02 | 8.73e-03 | 4.28e-03 | 4.28e-03 |
| 1 | 0.153 | 6.13e-02 | 8.24e-03 | 1.65e-03 | 8.83e-04 | 5.47e-04 | 5.47e-04 |
| 0.5 | 0.077 | 3.00e-02 | 2.15e-03 | 2.12e-04 | 9.04e-05 | 6.93e-05 | 6.93e-05 |
| 0.25 | 0.038 | 1.48e-02 | 5.47e-04 | 2.66e-05 | 1.01e-05 | 8.75e-06 | 8.75e-06 |

Order between the two finest dt, and least squares over the three finest:

- explicit 1: 1.02 (1.03)
- explicit 2: 1.97 (1.96)
- explicit 3: 2.99 (2.98)
- explicit 4: 3.17 (3.23)
- implicit 1: 2.99 (2.98)
- implicit 3: 2.98 (2.98)

The fine coupled run against the monolithic 1D run (CFL 0.1): 1.62e-08 of the peak.

## 2. Convergence in dt, CFD BDF 4 (max |Q - Q_fine| / peak)

Q_fine: implicit order 4 at dt = 0.03125 ms.

| dt [ms] | dt/tau | explicit 1 | explicit 2 | explicit 3 | explicit 4 | implicit 1 | implicit 4 |
|---|---|---|---|---|---|---|---|
| 4 | 0.613 | 3.43e-01 | 1.34e-01 | 5.67e-02 | 4.90e-02 | 2.13e-02 | 2.13e-02 |
| 2 | 0.307 | 1.31e-01 | 3.45e-02 | 1.04e-02 | 5.12e-03 | 1.80e-03 | 1.80e-03 |
| 1 | 0.153 | 6.17e-02 | 8.78e-03 | 1.47e-03 | 3.47e-04 | 1.17e-04 | 1.17e-04 |
| 0.5 | 0.077 | 3.00e-02 | 2.22e-03 | 1.93e-04 | 2.23e-05 | 7.53e-06 | 7.53e-06 |
| 0.25 | 0.038 | 1.48e-02 | 5.55e-04 | 2.46e-05 | 1.51e-06 | 5.47e-07 | 5.46e-07 |

Order between the two finest dt, and least squares over the three finest:

- explicit 1: 1.02 (1.03)
- explicit 2: 2.00 (1.99)
- explicit 3: 2.97 (2.95)
- explicit 4: 3.89 (3.92)
- implicit 1: 3.78 (3.87)
- implicit 4: 3.78 (3.87)

The fine coupled run against the monolithic 1D run (CFL 0.1): 4.41e-09 of the peak.

## 3. Explicit stability

Amplification per step of the step-to-step oscillation of the flow after the pulse: predicted (the model's largest root) and measured (a fit over 0.30-0.50 s, the steps fitted in brackets).

| CFD ddt | order | predicted limit dt/tau | run at dt/tau | predicted factor | measured factor | result |
|---|---|---|---|---|---|---|
| backward | 2 | 1.3333 | 1.2000 | 0.8294 | 0.8295 (24) | stable |
| backward | 2 | 1.3333 | 1.4667 | 1.1721 | 1.1721 (19) | unstable |
| backward | 3 | 0.5714 | 0.5143 | 0.8968 | 0.9162 (58) | stable |
| backward | 3 | 0.5714 | 0.6286 | 1.1042 | 1.1088 (47) | unstable |
| backward | 4 | 0.2667 | 0.2400 | 0.9273 | 0.9634 (126) | stable |
| backward | 4 | 0.2667 | 0.2933 | 1.0720 | 1.0705 (103) | unstable |
| BDF 3 | 3 | 0.9524 | 0.8571 | 0.8497 | 0.8517 (34) | stable |
| BDF 3 | 3 | 0.9524 | 1.0476 | 1.1508 | 1.1508 (28) | unstable |
| BDF 3 | 4 | 0.4444 | 0.4000 | 0.9077 | 0.9106 (76) | stable |
| BDF 3 | 4 | 0.4444 | 0.4889 | 1.0923 | 1.0603 (61) | unstable |
| BDF 4 | 1 | 1.6428 | 1.4785 | 0.9733 | 0.9851 (19) | stable |
| BDF 4 | 1 | 1.6428 | 1.8071 | 1.0245 | 1.0352 (15) | unstable |
| BDF 4 | 3 | 1.5238 | 1.3714 | 0.5888 | 0.6147 (21) | stable |
| BDF 4 | 3 | 1.5238 | 1.6762 | 1.2909 | 1.2909 (17) | unstable |
| BDF 4 | 4 | 0.7111 | 0.6400 | 0.8618 | 0.8459 (47) | stable |
| BDF 4 | 4 | 0.7111 | 0.7822 | 1.1371 | 1.1373 (38) | unstable |
