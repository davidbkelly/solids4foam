# The step inflow in the resonator: limiters and the unlimited schemes

1D only, linearised, hard coupling (inertance to ground), c0 = 3.5 m/s, step 1e-7 m3/s for 0.1 s, 0-0.5 s, CFL 0.5; the interface flow against the exact delay-equation solution (whose total variation over the window is 19.14 x the inflow).

| cells | scheme | overshoot / inflow | TV / TV exact | RMS error / inflow | max error / inflow |
|---|---|---|---|---|---|
| 100 | minmod | -0.003 | 0.859 | 1.55e-02 | 6.73e-02 |
| 100 | vanLeer | -0.001 | 0.929 | 1.08e-02 | 4.69e-02 |
| 100 | MC | -0.001 | 0.951 | 1.10e-02 | 4.98e-02 |
| 100 | none | -0.000 | 0.968 | 1.17e-02 | 5.08e-02 |
| 100 | kExact2 | +0.000 | 0.969 | 4.00e-03 | 2.29e-02 |
| 100 | kExact3 | +0.000 | 0.999 | 2.37e-03 | 2.00e-02 |
| 400 | minmod | -0.000 | 0.953 | 4.62e-03 | 2.53e-02 |
| 400 | vanLeer | -0.000 | 0.984 | 3.00e-03 | 2.31e-02 |
| 400 | MC | -0.000 | 0.992 | 3.13e-03 | 2.51e-02 |
| 400 | none | -0.000 | 1.003 | 3.48e-03 | 2.74e-02 |
| 400 | kExact2 | +0.000 | 0.996 | 8.86e-04 | 1.05e-02 |
| 400 | kExact3 | +0.000 | 1.002 | 4.69e-04 | 6.54e-03 |
