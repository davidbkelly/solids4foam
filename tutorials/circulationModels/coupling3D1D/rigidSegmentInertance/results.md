# Rigid segment (inertance) against the closed form

## Base case (dx = 1.25 mm)

c0 = 6.1721 m/s, Z0 = 2.0629e+07 Pa s/m3, pipe 8 cm: tau = 6.481 ms, sigma = 7.071 ms, incident peak 20.629 mPa (inflow x 0.001).

| quantity | closed form peak | 1D peak | max error / incident peak |
|---|---|---|---|
| pipe flow Q [ul/s] | 0.8043 | 0.8035 | 3.66e-03 |
| P1, pipe inlet [mPa] | 27.9269 | 27.9514 | 1.20e-02 |
| P2 = Pt, pipe outlet [mPa] | 16.5908 | 16.5747 | 3.66e-03 |
| proximal midpoint [mPa] | 20.6288 | 20.6194 | 9.03e-03 |
| distal midpoint [mPa] | 16.5912 | 16.5638 | 5.45e-03 |

Reflected peak 8.4127 mPa against 8.3535 mPa; energy (reflected + transmitted)/incident from the closed form is 1 to round-off (benchmarks3D1D.py rigid).

## Order of accuracy (both vessels at t = 0.6 s)

| cells per vessel | dx [mm] | L1 | order | L2 | order | max | order |
|---|---|---|---|---|---|---|---|
| 200 | 10 | 2.030e-02 |  | 5.314e-02 |  | 2.718e-01 |  |
| 400 | 5 | 6.956e-03 | 1.55 | 2.115e-02 | 1.33 | 1.262e-01 | 1.11 |
| 800 | 2.5 | 1.847e-03 | 1.91 | 5.931e-03 | 1.83 | 3.903e-02 | 1.69 |
| 1600 | 1.25 | 4.668e-04 | 1.98 | 1.511e-03 | 1.97 | 1.022e-02 | 1.93 |
| 3200 | 0.625 | 1.174e-04 | 1.99 | 3.804e-04 | 1.99 | 2.586e-03 | 1.98 |
| 6400 | 0.3125 | 2.978e-05 | 1.98 | 9.672e-05 | 1.98 | 6.580e-04 | 1.97 |

Least-squares order over the 4 finest meshes: L1 1.99, L2 1.98, max 1.97.
