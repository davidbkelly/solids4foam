# Businaro (2026) pressure Gaussian pulse: 3D-1D coupling against the full FSI tube

Inlet pressure 2000 Pa exp(-1e4 (t - 0.05)^2), R = 1 cm, wall 1 mm, E = 700 kPa, Saint Venant-Kirchhoff, rho = 1060, mu = 4 mPa s; BDF2, window 0-0.1 s. Reference: full3D, 45 cm. RMS and max pressure difference over the window / 2 kPa; RMS flow difference / the reference's peak flow at that section; peak pressure ratio run / reference. z <= 10 cm is 3D in the coupled runs, beyond it the 1D vessel.

| run | z = 5 | z = 9 | z = 10 | z = 11 | z = 15 | z = 20 | z = 25 | z = 30 |
|---|---|---|---|---|---|---|---|---|
| 3D FSI, 30 cm, resistance Rt at the outlet: RMS dP | 9.92e-03 | 1.06e-02 | 1.08e-02 | 1.12e-02 | 9.67e-03 | 6.19e-03 | 5.08e-03 | 1.48e-02 |
| max dP | 2.76e-02 | 2.88e-02 | 2.85e-02 | 2.92e-02 | 2.47e-02 | 1.73e-02 | 1.48e-02 | 6.14e-02 |
| RMS dQ | 7.69e-03 | 7.93e-03 | 6.66e-03 | 4.96e-03 | 3.07e-03 | 5.21e-03 | 7.08e-03 | 1.79e-02 |
| peak P ratio | 1.004 | 1.015 | 1.017 | 1.017 | 1.004 | 0.993 | 0.997 | 1.065 |
| 3D 10 cm + 1D 20 cm, implicit order 2: RMS dP | 8.97e-03 | 1.77e-02 | 2.11e-02 | 2.54e-02 | 4.39e-02 | 6.50e-02 | 7.75e-02 | 1.01e-01 |
| max dP | 2.60e-02 | 5.25e-02 | 5.99e-02 | 6.83e-02 | 1.04e-01 | 1.86e-01 | 2.83e-01 | 3.97e-01 |
| RMS dQ | 1.46e-02 | 2.12e-02 | 2.10e-02 | 2.42e-02 | 4.19e-02 | 6.59e-02 | 7.87e-02 | 7.53e-02 |
| peak P ratio | 1.009 | 1.023 | 1.027 | 1.026 | 1.018 | 1.015 | 1.024 | 1.051 |
| 3D 10 cm + 1D 20 cm, explicit order 1: RMS dP | 1.03e-02 | 1.62e-02 | 1.99e-02 | 2.53e-02 | 4.38e-02 | 6.50e-02 | 7.75e-02 | 1.01e-01 |
| max dP | 3.03e-02 | 4.89e-02 | 5.73e-02 | 6.80e-02 | 1.04e-01 | 1.86e-01 | 2.83e-01 | 3.97e-01 |
| RMS dQ | 1.37e-02 | 2.24e-02 | 2.22e-02 | 2.44e-02 | 4.18e-02 | 6.59e-02 | 7.87e-02 | 7.53e-02 |
| peak P ratio | 1.008 | 1.020 | 1.027 | 1.026 | 1.018 | 1.015 | 1.024 | 1.051 |
| 1D alone, 30 cm: RMS dP | 2.43e-02 | 4.62e-02 | 5.19e-02 | 5.76e-02 | 7.86e-02 | 1.00e-01 | 1.10e-01 | 1.28e-01 |
| max dP | 5.72e-02 | 1.16e-01 | 1.33e-01 | 1.49e-01 | 2.15e-01 | 2.95e-01 | 3.88e-01 | 4.95e-01 |
| RMS dQ | 2.92e-02 | 4.60e-02 | 5.09e-02 | 5.59e-02 | 7.67e-02 | 1.01e-01 | 1.12e-01 | 1.02e-01 |
| peak P ratio | 1.009 | 1.017 | 1.017 | 1.017 | 1.009 | 1.006 | 1.015 | 1.042 |

Wave speed from the half-peak arrivals (least squares):

- full3D, z = 2-28 cm: 5.667 m/s (the 1D vessel's c0 = 6.024 m/s)
- full3D_Rt, z = 2-9 cm: 5.743 m/s
- coupled, z = 12-28 cm: 6.231 m/s
- coupled_explicit, z = 12-28 cm: 6.231 m/s
- only1D, z = 2-28 cm: 6.229 m/s

Interface mismatch (Businaro Eq. 7.12): largest |e_P| and |e_Q| over the window, and relative to the pulse and to the largest interface flow:

| run | max abs e_P [Pa] | / 2 kPa | max abs e_Q [m3/s] | / max Q | max Q at the interface [m3/s] |
|---|---|---|---|---|---|
| 3D 10 cm + 1D 20 cm, implicit order 2 | 5.280e-06 | 2.64e-09 | 2.452e-12 | 2.38e-08 | 1.0283e-04 |
| 3D 10 cm + 1D 20 cm, explicit order 1 | 5.270e-06 | 2.63e-09 | 9.761e-07 | 9.49e-03 | 1.0282e-04 |
| Businaro (2026), explicit (lifex) | 4.452e-01 | 2.23e-04 | 5.757e-06 | 5.57e-02 | 1.0330e-04 |

The full tube's largest flow at z = 10 cm: 1.0673e-04 m3/s, peak pressure 1951.9 Pa.
