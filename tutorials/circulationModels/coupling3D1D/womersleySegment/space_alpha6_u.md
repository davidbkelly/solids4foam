# Order in space of the Womersley segment, alpha = 6 (uniform mesh family)

Mesh refined as nr x nz; 0-1.0 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised vessels at fourth order. Errors against the exact transient solution: port pressures (max of the two) and flow, RMS over the window relative to the RMS of the exact signal; mid-length velocity profile, max over eight instants of max |u - u_exact| / max |u_exact|, from the axis to the wall cell (inside the wall cell the sampled values converge at first order only). The wedge has no wall grading. With its flat face at the distance R from the axis and the patch flux scaled by pi/tan(theta/2), the discrete problem is the circular pipe of radius R exactly; the exact solution is Womersley's for R, independent of the mesh.

| nr x nz | P ports | order | Q | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 1.37e-04 |  | 1.37e-04 |  | 4.84e-03 |  | 5.8 | 84 |
| 50 x 100 | 3.43e-05 | 2.00 | 3.43e-05 | 2.00 | 1.21e-03 | 2.00 | 8.3 | 394 |
| 100 x 200 | 8.58e-06 | 2.00 | 8.58e-06 | 2.00 | 2.99e-04 | 2.01 | 14.8 | 4398 |
