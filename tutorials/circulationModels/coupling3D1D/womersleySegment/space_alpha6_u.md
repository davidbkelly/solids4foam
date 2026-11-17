# Order in space of the Womersley segment, alpha = 6 (uniform mesh family)

Mesh refined as nr x nz; 0-1.0 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised vessels at fourth order. Errors against the exact transient solution: port pressures (max of the two) and flow, RMS over the window relative to the RMS of the exact signal; mid-length velocity profile, max over eight instants of max |u - u_exact| / max |u_exact|. The wedge has no wall grading and is scaled by sqrt(theta/sin theta), so that its cross-section is pi R^2. The exact solution is that of the problem the wedge discretises: the pipe of radius R_c = 9.99683 mm (the flat face's distance from the axis) with every flux times k = 1.000635 (analyse.py); the profile is compared for r <= R_c.

| nr x nz | P ports | order | Q | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 1.37e-04 |  | 1.37e-04 |  | 4.64e-02 |  | 5.8 | 69 |
| 50 x 100 | 3.43e-05 | 1.99 | 3.43e-05 | 1.99 | 2.54e-02 | 0.87 | 8.4 | 322 |
| 100 x 200 | 8.59e-06 | 2.00 | 8.59e-06 | 2.00 | 1.48e-02 | 0.78 | 14.8 | 3508 |
| 200 x 400 | 2.15e-06 | 2.00 | 2.15e-06 | 2.00 | 8.54e-03 | 0.79 | 35.9 | 60141 |
