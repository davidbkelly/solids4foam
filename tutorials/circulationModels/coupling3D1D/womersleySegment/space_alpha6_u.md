# Order in space of the Womersley segment, alpha = 6 (uniform mesh family)

Mesh refined as nr x nz; 0-1.0 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised vessels at fourth order. Errors against the exact transient solution: port pressures (max of the two) and flow, RMS over the window relative to the RMS of the exact signal; mid-length velocity profile, max over eight instants of max |u - u_exact| / max |u_exact|. The wedge is scaled by sqrt(theta/sin theta) (its flat wall then spans the cross-section pi R^2) and has no wall grading, so the exact solution is that of the pipe of radius R, independent of the mesh.

| nr x nz | P ports | order | Q | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 1.35e-04 |  | 1.35e-04 |  | 4.49e-02 |  | 5.8 | 69 |
| 50 x 100 | 3.36e-05 | 2.01 | 3.36e-05 | 2.01 | 2.39e-02 | 0.91 | 8.4 | 322 |
| 100 x 200 | 1.05e-05 | 1.68 | 1.05e-05 | 1.68 | 1.33e-02 | 0.84 | 14.8 | 3508 |
