# Order in space of the Womersley segment, alpha = 6

Mesh refined as nr x nz; 0-1.0 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised vessels at fourth order. Errors against the exact transient solution: port pressures (max of the two) and flow, RMS over the window relative to the RMS of the exact signal; mid-length velocity profile, max over eight instants of max |u - u_exact| / max |u_exact|. The exact solution is that of the problem the wedge discretises: the pipe of radius R_c = R cos(theta/2) = 9.99048 mm with every flux times k = 1.000635 (analyse.py).

| nr x nz | P ports | order | Q | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 3.11e-05 |  | 3.11e-05 |  | 6.74e-03 |  | 8.9 | 78 |
| 50 x 100 | 7.80e-06 | 1.99 | 7.80e-06 | 1.99 | 8.56e-04 | 2.98 | 18.6 | 681 |
| 100 x 200 | 1.96e-06 | 2.00 | 1.96e-06 | 2.00 | 1.86e-04 | 2.20 | 45.4 | 12276 |
