# Womersley segment at physiological amplitude (numerical reference)

alpha = 6, incident wave 5000 Pa, 0-1.0 s from rest, BDF 4 / coupling order 4 at dt = 0.25 ms, nonlinear vessels (kExact 3, RK4, 800 cells). The 100 x 200 run is the reference; the 50 x 100 run against it (RMS / RMS of the signal):

- to3D.P: 5.11e-06
- to3D.Q: 5.68e-06
- from3D.P: 5.75e-06
- mid-length profile, max over 9 instants: 7.62e-04

Against the linear exact solution scaled to 5 kPa (what the linear model would give), RMS / RMS:

- to3D.P: 5.98e-02
- to3D.Q: 1.06e-01
- from3D.P: 5.32e-02
