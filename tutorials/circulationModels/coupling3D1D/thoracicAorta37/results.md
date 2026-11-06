# The thoracic aorta II (s17) of the 37-artery network as a 3D segment

Last cycle (15.713-16.540 s), period 0.827 s.

## AorticArchII

| run | peak P [kPa] | min P [kPa] | peak Q [ml/s] | RMS dP vs original / pulse | RMS dQ vs original / peak Q | RMS dP vs measured | RMS dQ vs measured |
|---|---|---|---|---|---|---|---|
| original | 15.01 | 10.52 | 179.2 | 0.00e+00 | 0.00e+00 | 8.42e-02 | 1.34e-01 |
| stubs1D | 15.01 | 10.52 | 179.2 | 4.26e-07 | 1.65e-06 | 8.42e-02 | 1.34e-01 |
| stubsInertance | 15.38 | 10.47 | 184.1 | 4.42e-02 | 5.04e-02 | 8.64e-02 | 1.23e-01 |
| rigid3D | 15.38 | 10.47 | 182.9 | 4.50e-02 | 5.14e-02 | 8.56e-02 | 1.24e-01 |
| measured | 15.23 | 10.57 | 164.4 | 8.76e-02 | 1.23e-01 | 0.00e+00 | 0.00e+00 |

rigid3D against the inertance limit: RMS dP / pulse 1.99e-03, RMS dQ / peak Q 3.68e-03

## ThoracicAortaII

| run | peak P [kPa] | min P [kPa] | peak Q [ml/s] | RMS dP vs original / pulse | RMS dQ vs original / peak Q | RMS dP vs measured | RMS dQ vs measured |
|---|---|---|---|---|---|---|---|
| original | 15.21 | 10.48 | 127.9 | 0.00e+00 | 0.00e+00 | 1.18e-01 | 2.91e-01 |
| stubs1D | 15.21 | 10.48 | 127.9 | 6.74e-07 | 1.98e-06 | 1.18e-01 | 2.91e-01 |
| stubsInertance | 15.54 | 10.41 | 123.2 | 4.28e-02 | 3.92e-02 | 1.19e-01 | 2.85e-01 |
| rigid3D | 15.60 | 10.45 | 120.8 | 4.37e-02 | 4.21e-02 | 1.22e-01 | 2.85e-01 |
| measured | 15.36 | 10.31 | 79.0 | 1.26e-01 | 1.79e-01 | 0.00e+00 | 0.00e+00 |

rigid3D against the inertance limit: RMS dP / pulse 1.48e-02, RMS dQ / peak Q 6.54e-03

## RIliacFemoralII

| run | peak P [kPa] | min P [kPa] | peak Q [ml/s] | RMS dP vs original / pulse | RMS dQ vs original / peak Q | RMS dP vs measured | RMS dQ vs measured |
|---|---|---|---|---|---|---|---|
| original | 16.29 | 10.09 | 12.4 | 0.00e+00 | 0.00e+00 | 2.99e-01 | 3.78e-01 |
| stubs1D | 16.29 | 10.09 | 12.4 | 1.30e-06 | 3.50e-06 | 2.99e-01 | 3.78e-01 |
| stubsInertance | 16.80 | 9.91 | 13.4 | 4.45e-02 | 5.05e-02 | 3.09e-01 | 4.12e-01 |
| rigid3D | 16.79 | 9.90 | 13.1 | 4.46e-02 | 5.10e-02 | 3.11e-01 | 4.12e-01 |
| measured | 16.66 | 10.54 | 8.2 | 2.97e-01 | 2.50e-01 | 0.00e+00 | 0.00e+00 |

rigid3D against the inertance limit: RMS dP / pulse 4.73e-03, RMS dQ / peak Q 8.09e-03

## Cost of the coupled run

1654 steps, 5.1 outer iterations per step, 10.1 network re-integrations per step, 273 s wall time
