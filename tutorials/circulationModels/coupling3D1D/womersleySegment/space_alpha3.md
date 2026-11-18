# Order in space of the Womersley segment, alpha = 3

Mesh refined as nr x nz; 0-1.0 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised vessels at fourth order. Errors against the exact transient solution: port pressures (max of the two) and flow, RMS over the window relative to the RMS of the exact signal; mid-length velocity profile, max over eight instants of max |u - u_exact| / max |u_exact|, from the axis to the wall cell (inside the wall cell the sampled values converge at first order only). With its flat face at the distance R from the axis and the patch flux scaled by pi/tan(theta/2), the discrete problem is the circular pipe of radius R exactly; the exact solution is Womersley's for R, independent of the mesh.

| nr x nz | P ports | order | Q | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 5.80e-05 |  | 5.80e-05 |  | 7.69e-03 |  | 16.1 | 237 |
| 50 x 100 | 1.43e-05 | 2.02 | 1.43e-05 | 2.02 | 2.19e-03 | 1.81 | 43.0 | 2145 |
