# Order in space of the Womersley segment, alpha = 12

Mesh refined as nr x nz; 0-1.0 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised vessels at fourth order. Errors against the exact transient solution: port pressures (max of the two) and flow, RMS over the window relative to the RMS of the exact signal; mid-length velocity profile, max over eight instants of max |u - u_exact| / max |u_exact|. The exact solution is that of the circular pipe of the wedge's cross-section (flat wall: 0.5 R^2 sin(theta) per wedge, R_eq = 9.99365 mm); against the pipe of radius R the comparison floors at 6e-5 whatever the mesh.

| nr x nz | P ports | order | Q | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 2.48e-05 |  | 2.48e-05 |  | 1.10e-02 |  | 6.2 | 57 |
| 50 x 100 | 5.51e-06 | 2.17 | 5.51e-06 | 2.17 | 2.50e-03 | 2.14 | 10.2 | 391 |
