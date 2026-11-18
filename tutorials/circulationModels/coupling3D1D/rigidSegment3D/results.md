# Rigid segment, 3D-1D: pimpleFoam pipe coupled to 1D vessels

c0 = 6.1721 m/s, Zc = 2.0629e+07 Pa s/m3; pipe 8 cm, A3D = 3.141593e-04 m2 (1.0000 pi R^2: the 32-sided patch), tau = 6.4807 ms; incident peak 20.6288 mPa (inflow x 1e-3); 1D dx = 1.25 mm, SSPRK3 (CFL 0.5 with backward, 0.1 with BDF 3 and BDF 4).

## 1. Against the closed form (max error / incident peak)

| run | Q (3D flux) | P1, pipe inlet | P2, pipe outlet |
|---|---|---|---|
| monolithic 1D, CFL 0.5 | 1.22e-03 | 4.04e-03 | 1.22e-03 |
| monolithic 1D, CFL 0.1 | 1.22e-03 | 4.04e-03 | 1.22e-03 |
| backward: explicit 2, dt 4 ms | 9.61e-02 | 7.94e-02 | 7.52e-02 |
| backward: explicit 2, dt 2 ms | 2.90e-02 | 2.46e-02 | 2.04e-02 |
| backward: explicit 2, dt 1 ms | 8.36e-03 | 9.30e-03 | 5.94e-03 |
| backward: explicit 2, dt 0.5 ms | 2.82e-03 | 5.30e-03 | 2.24e-03 |
| backward: explicit 2, dt 0.25 ms | 1.43e-03 | 4.35e-03 | 1.45e-03 |
| backward: implicit 2, dt 4 ms | 3.60e-02 | 3.18e-02 | 3.60e-02 |
| backward: implicit 2, dt 2 ms | 1.03e-02 | 6.89e-03 | 1.03e-02 |
| backward: implicit 2, dt 1 ms | 2.64e-03 | 3.47e-03 | 2.64e-03 |
| backward: implicit 2, dt 0.5 ms | 1.09e-03 | 3.51e-03 | 1.09e-03 |
| backward: implicit 2, dt 0.25 ms | 1.12e-03 | 3.89e-03 | 1.12e-03 |
| backward: implicit 2, dt 0.03125 ms | 1.22e-03 | 4.04e-03 | 1.22e-03 |
| BDF3: explicit 3, dt 4 ms | 5.66e-02 | 5.25e-02 | 5.25e-02 |
| BDF3: explicit 3, dt 2 ms | 9.58e-03 | 9.93e-03 | 1.05e-02 |
| BDF3: explicit 3, dt 1 ms | 1.29e-03 | 3.90e-03 | 2.32e-03 |
| BDF3: explicit 3, dt 0.5 ms | 1.10e-03 | 3.97e-03 | 1.36e-03 |
| BDF3: explicit 3, dt 0.25 ms | 1.21e-03 | 4.03e-03 | 1.24e-03 |
| BDF3: implicit 3, dt 4 ms | 2.05e-02 | 2.05e-02 | 2.05e-02 |
| BDF3: implicit 3, dt 2 ms | 2.46e-03 | 5.19e-03 | 2.46e-03 |
| BDF3: implicit 3, dt 1 ms | 8.54e-04 | 4.19e-03 | 8.54e-04 |
| BDF3: implicit 3, dt 0.5 ms | 1.18e-03 | 4.07e-03 | 1.18e-03 |
| BDF3: implicit 3, dt 0.25 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| BDF3: implicit 3, dt 0.03125 ms | 1.22e-03 | 4.04e-03 | 1.22e-03 |
| BDF4: explicit 4, dt 4 ms | 3.93e-02 | 4.78e-02 | 4.82e-02 |
| BDF4: explicit 4, dt 2 ms | 3.17e-03 | 2.85e-03 | 6.29e-03 |
| BDF4: explicit 4, dt 1 ms | 9.95e-04 | 3.66e-03 | 1.31e-03 |
| BDF4: explicit 4, dt 0.5 ms | 1.21e-03 | 4.02e-03 | 1.22e-03 |
| BDF4: explicit 4, dt 0.25 ms | 1.22e-03 | 4.04e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 4 ms | 1.62e-02 | 1.92e-02 | 1.62e-02 |
| BDF4: implicit 4, dt 2 ms | 1.45e-03 | 5.41e-03 | 1.45e-03 |
| BDF4: implicit 4, dt 1 ms | 1.22e-03 | 4.14e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 0.5 ms | 1.22e-03 | 4.05e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 0.25 ms | 1.22e-03 | 4.04e-03 | 1.22e-03 |
| BDF4: implicit 4, dt 0.03125 ms | 1.22e-03 | 4.04e-03 | 1.22e-03 |

The monolithic rows are the 1D discretisation error alone (dx = 1.25 mm); a coupled run cannot do better than them. Only the coupling of the time scheme's own order is listed.

Mass through the rigid pipe (implicit 2, 0.5 ms): max |Qin - Qout| = 1.0e-12 of the peak; outer iterations 10.4 per step on average.

## 2. Convergence in dt, CFD backward (max |Q - Q_fine| / peak)

Q_fine: implicit order 2 at dt = 0.03125 ms.

| dt [ms] | dt/tau | explicit 1 | explicit 2 | explicit 3 | explicit 4 | implicit 2 |
|---|---|---|---|---|---|---|
| 4 | 0.617 | 2.23e-01 | 1.18e-01 | 9.84e-02 | unstable | 4.46e-02 |
| 2 | 0.309 | 1.15e-01 | 3.48e-02 | 2.35e-02 | 1.59e-02 | 1.27e-02 |
| 1 | 0.154 | 5.80e-02 | 9.34e-03 | 4.94e-03 | 3.55e-03 | 3.42e-03 |
| 0.5 | 0.077 | 2.92e-02 | 2.40e-03 | 1.07e-03 | 8.74e-04 | 8.71e-04 |
| 0.25 | 0.039 | 1.46e-02 | 6.05e-04 | 2.43e-04 | 2.18e-04 | 2.17e-04 |

Order between the two finest dt, and least squares over the three finest:

- explicit 1: 1.00 (0.99)
- explicit 2: 1.99 (1.97)
- explicit 3: 2.14 (2.17)
- explicit 4: 2.01 (2.01)
- implicit 2: 2.00 (1.99)

The fine coupled run against the monolithic 1D run (CFL 0.5): 3.41e-06 of the peak.

## 2. Convergence in dt, CFD BDF 3 (max |Q - Q_fine| / peak)

Q_fine: implicit order 3 at dt = 0.03125 ms.

| dt [ms] | dt/tau | explicit 1 | explicit 2 | explicit 3 | explicit 4 | implicit 1 | implicit 3 |
|---|---|---|---|---|---|---|---|
| 4 | 0.617 | 2.83e-01 | 1.01e-01 | 7.03e-02 | 6.44e-02 | 2.65e-02 | 2.65e-02 |
| 2 | 0.309 | 1.27e-01 | 3.05e-02 | 1.19e-02 | 8.76e-03 | 4.28e-03 | 4.28e-03 |
| 1 | 0.154 | 6.15e-02 | 8.29e-03 | 1.65e-03 | 8.86e-04 | 5.48e-04 | 5.48e-04 |
| 0.5 | 0.077 | 3.01e-02 | 2.16e-03 | 2.13e-04 | 9.07e-05 | 6.94e-05 | 6.94e-05 |
| 0.25 | 0.039 | 1.48e-02 | 5.50e-04 | 2.67e-05 | 1.01e-05 | 8.79e-06 | 8.79e-06 |

Order between the two finest dt, and least squares over the three finest:

- explicit 1: 1.02 (1.03)
- explicit 2: 1.97 (1.96)
- explicit 3: 2.99 (2.98)
- explicit 4: 3.16 (3.23)
- implicit 1: 2.98 (2.98)
- implicit 3: 2.98 (2.98)

The fine coupled run against the monolithic 1D run (CFL 0.1): 1.89e-08 of the peak.

## 2. Convergence in dt, CFD BDF 4 (max |Q - Q_fine| / peak)

Q_fine: implicit order 4 at dt = 0.03125 ms.

| dt [ms] | dt/tau | explicit 1 | explicit 2 | explicit 3 | explicit 4 | implicit 1 | implicit 4 |
|---|---|---|---|---|---|---|---|
| 4 | 0.617 | 3.45e-01 | 1.35e-01 | 5.69e-02 | 4.92e-02 | 2.14e-02 | 2.14e-02 |
| 2 | 0.309 | 1.32e-01 | 3.47e-02 | 1.04e-02 | 5.16e-03 | 1.80e-03 | 1.80e-03 |
| 1 | 0.154 | 6.19e-02 | 8.83e-03 | 1.48e-03 | 3.49e-04 | 1.17e-04 | 1.17e-04 |
| 0.5 | 0.077 | 3.01e-02 | 2.23e-03 | 1.94e-04 | 2.25e-05 | 7.57e-06 | 7.57e-06 |
| 0.25 | 0.039 | 1.48e-02 | 5.58e-04 | 2.47e-05 | 1.55e-06 | 5.82e-07 | 5.81e-07 |

Order between the two finest dt, and least squares over the three finest:

- explicit 1: 1.02 (1.03)
- explicit 2: 2.00 (1.99)
- explicit 3: 2.97 (2.95)
- explicit 4: 3.86 (3.91)
- implicit 1: 3.70 (3.83)
- implicit 4: 3.70 (3.83)

The fine coupled run against the monolithic 1D run (CFL 0.1): 5.30e-09 of the peak.

## 3. Explicit stability

Amplification per step of the step-to-step oscillation of the flow after the pulse: predicted (the model's largest root) and measured (a fit over 0.30-0.50 s, the steps fitted in brackets).

| CFD ddt | order | predicted limit dt/tau | run at dt/tau | predicted factor | measured factor | result |
|---|---|---|---|---|---|---|
| backward | 2 | 1.3333 | 1.2000 | 0.8294 | 0.8297 (24) | stable |
| backward | 2 | 1.3333 | 1.4667 | 1.1721 | 1.1721 (20) | unstable |
| backward | 3 | 0.5714 | 0.5143 | 0.8968 | 0.9156 (58) | stable |
| backward | 3 | 0.5714 | 0.6286 | 1.1042 | 1.1065 (48) | unstable |
| backward | 4 | 0.2667 | 0.2400 | 0.9273 | 0.9634 (127) | stable |
| backward | 4 | 0.2667 | 0.2933 | 1.0720 | 1.0716 (104) | unstable |
| BDF 3 | 3 | 0.9524 | 0.8571 | 0.8497 | 0.8473 (34) | stable |
| BDF 3 | 3 | 0.9524 | 1.0476 | 1.1508 | 1.1508 (28) | unstable |
| BDF 3 | 4 | 0.4444 | 0.4000 | 0.9077 | 0.9204 (76) | stable |
| BDF 3 | 4 | 0.4444 | 0.4889 | 1.0923 | 1.0946 (62) | unstable |
| BDF 4 | 1 | 1.6428 | 1.4785 | 0.9733 | 0.9783 (19) | stable |
| BDF 4 | 1 | 1.6428 | 1.8071 | 1.0245 | 1.0734 (16) | unstable |
| BDF 4 | 3 | 1.5238 | 1.3714 | 0.5888 | 0.6583 (21) | stable |
| BDF 4 | 3 | 1.5238 | 1.6762 | 1.2909 | 1.2909 (17) | unstable |
| BDF 4 | 4 | 0.7111 | 0.6400 | 0.8618 | 0.8701 (47) | stable |
| BDF 4 | 4 | 0.7111 | 0.7822 | 1.1371 | 1.1376 (38) | unstable |
