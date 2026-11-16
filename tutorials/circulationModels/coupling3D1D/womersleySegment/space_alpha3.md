# Order in space of the Womersley segment, alpha = 3

Mesh refined as nr x nz; 0-1.0 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised vessels at fourth order. Errors against the exact transient solution: port pressures (max of the two) and flow, RMS over the window relative to the RMS of the exact signal; mid-length velocity profile, max over eight instants of max |u - u_exact| / max |u_exact|. The exact solution is that of the circular pipe of the wedge's cross-section (flat wall: 0.5 R^2 sin(theta) per wedge, R_eq = 9.99365 mm); against the pipe of radius R the comparison floors at 6e-5 whatever the mesh.

| nr x nz | P ports | order | Q | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 3.43e-05 |  | 3.43e-05 |  | 7.68e-03 |  | 15.6 | 139 |
| 50 x 100 | 1.35e-05 | 1.35 | 1.35e-05 | 1.35 | 2.18e-03 | 1.82 | 42.1 | 1531 |
