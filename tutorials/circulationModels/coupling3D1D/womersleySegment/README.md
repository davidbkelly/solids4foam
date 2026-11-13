# The Womersley segment: `womersleySegment`

Prepared by David Kelly

---

## Tutorial Aims

- Verify the 3D-1D coupling against an exact *viscous* solution: a rigid
  no-slip pipe between two 1D vessels, driven by a single-frequency wave, is
  in the linear regime a transmission line with the pipe as the series
  impedance of Womersley (1955), so the reflected and transmitted waves, the
  port pressures and flow, and the velocity profile inside the pipe are all
  in closed form.
- Extend the rigid-segment benchmark (`rigidSegment3D`, inviscid plug flow,
  an inductor) to a 3D side that develops a real velocity profile, where
  the coupler still exchanges only the mean pressure and the flow.

---

## Case Overview

The two 1 m vessels of the rigid-segment benchmark (R = 1 cm, c0 = 6.17 m/s,
inviscid, 400 cells each) joined by a rigid no-slip pipe, R = 1 cm, length
10 cm, an axisymmetric 5-degree wedge in pimpleFoam (50 radial cells graded
10:1 towards the wall, 100 axial cells; `flowScale 72` on both ports),
`backward` at dt = 1 ms, implicit coupling of order 2. The proximal inlet
is a `waveInlet`: a non-reflecting end that sends in
p(t) = 50 Pa sin(2 pi t), ramped over the first period (velocity 7.7 mm/s,
Re <= 35: linear), and absorbs what the pipe reflects; the distal end is
absorbing. At f = 1 Hz the viscosity sets the Womersley number
alpha = R sqrt(w/nu): nu = w R^2/alpha^2 = 7.0e-5, 1.75e-5, 4.36e-6 m2/s for
alpha = 3, 6, 12. Each run lasts three times the slowest viscous decay
R^2/(j01^2 nu) plus three periods (4, 6 and 15 s); the last two periods are
fitted with a sinusoid.

**The exact solution.** For an oscillating axial pressure gradient the
fully developed profile is Womersley's, u(r) ~ 1 - J0(L r/R)/J0(L),
L = i^(3/2) alpha, and the flow is Q = (G A/(i w rho))(1 - F10),
F10 = 2 J1(L)/(L J0(L)). The pipe is therefore the series impedance

    Z(w) = i w rho l / (A (1 - F10)),

the inertance i w rho l/A as alpha -> infinity and Poiseuille's resistance
8 mu l/(pi R^4) as alpha -> 0. Between two lines of impedance Zc = rho c0/A
it reflects Gamma = Z/(Z + 2 Zc) and transmits T = 2 Zc/(Z + 2 Zc), so for
the incident wave P (the inlet's p(t) delayed by L1D/c0): P_in = (1 + Gamma)
P, Q = T P/Zc, P_out = T P. `analyse.py` compares the amplitude and phase of
the fitted port signals with these, tabulates the inviscid inertance's
values alongside, and compares the velocity profile at mid-length, sampled
at eight phases of the last period, with Womersley's profile for the fitted
Q.

```bash
./Allrun            # the three runs together (about an hour), then analyse.py
```

**PIMPLE settings.** The settings of `rigidSegment3D` (momentum predictor,
one corrector, no relaxation) make the outer loop diverge within the first
step on this viscous, wall-graded mesh, and GAMG runs to its iteration cap
on the wedge; `base/system/fvSolution` uses no momentum predictor, two
correctors and PCG, and the outer loop converges slowly (about 15
iterations to a pressure residual of 1e-7; at alpha = 3 it stalls at 4e-7
and runs to its 50-iteration cap, which costs time, not accuracy).

---

## Results

See `results.md`, `plots/ports.png` and `plots/profiles.png`.

| alpha | Z [Pa s/m3] | \|Gamma\| (inviscid) | \|T\| (inviscid) | amplitude error, P_in / Q / P_out | phase error [deg] | profile error |
|---|---|---|---|---|---|---|
| 3 | 1.99e6 + 2.76e6 i | 0.0786 (0.0508) | 0.952 (0.999) | 1.3e-4 / 1.5e-4 / 1.5e-4 | 0.003 / -0.004 / -0.004 | 2.3e-3 |
| 6 | 6.95e5 + 2.56e6 i | 0.0632 (0.0508) | 0.982 (0.999) | 8.5e-5 / 8.8e-5 / 8.9e-5 | 0.004 / -0.004 / -0.004 | 4.6e-3 |
| 12 | 2.94e5 + 2.34e6 i | 0.0568 (0.0508) | 0.991 (0.999) | 7.1e-5 / 7.1e-5 / 7.2e-5 | 0.004 / -0.004 / -0.004 | 8.9e-3 |

The coupled run reproduces the exact periodic state to 1e-4 in amplitude
and 0.004 degrees in phase at every Womersley number, and the mid-length
velocity profile matches Womersley's to 0.2-0.9 % of its peak (the Stokes
layer at alpha = 12 is 1.2 mm thick on a wall cell of 0.05 mm). Viscosity
is not a small effect here: at alpha = 3 the pipe's resistance is as large
as its reactance, it reflects 7.9 % of the wave instead of 5.1 % and
transmits 95 % instead of 99.9 %, and the runs follow the viscous values,
not the inviscid ones. The entrance region at the inlet port leaves no
measurable trace at this length (l/R = 10).

What this adds to `rigidSegment3D`: the coupling is exact for a 3D side
whose velocity profile is not uniform, i.e. the exchange of mean pressure
and flow at the ports loses nothing when the 3D flow is a genuine viscous
solution.

## Order in time

`../orderStudy.py --case womersleySegment` repeats the alpha = 6 run from
rest over the first 0.6 s of the ramped incident wave at dt = 4, 2, 1 and
0.5 ms, with the fluid ddt scheme and the implicit coupling order both at
q = 1 (`Euler`), 2 (`backward`), 3 (`BDF 3`) and 4 (`BDF 4`), the 1D
vessels at CFL 0.1 and the PIMPLE loop to a pressure residual of 1e-10.
Observed orders by self-convergence (`orders.md`, `plots/orders.png`), the
RMS of the port flows and pressures (the four signals give the same orders
to two decimals, the pipe being rigid):

| dt triple [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 4 / 2 / 1 | 0.99 | 1.97 | 2.79 | 3.31 |
| 2 / 1 / 0.5 | 1.00 | 1.99 | 2.92 | 3.76 |

The finest Euler, backward and BDF 3 runs differ from the finest BDF 4 run
by 1.4e-4, 4.8e-7 and 1.4e-8 of the flow: the 3D viscous solution and the
coupling converge together at the scheme's order, down to the 1e-8 level
where the converged runs agree.
