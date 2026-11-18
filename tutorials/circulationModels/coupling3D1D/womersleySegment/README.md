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
10:1 towards the wall, 100 axial cells; its flat outer face at the distance
R from the axis and `flowScale` pi/tan(2.5 deg) = 71.954 on both ports, so
that the coupled problem is exactly the circular pipe of radius R, see
"The wedge's geometry" below),
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
| 3 | 1.99e6 + 2.76e6 i | 0.0786 (0.0508) | 0.952 (0.999) | 1.1e-5 / 1.3e-5 / 1.2e-5 | < 0.0005 | 3.1e-3 |
| 6 | 6.95e5 + 2.56e6 i | 0.0632 (0.0508) | 0.982 (0.999) | 6.9e-6 / 8.0e-6 / 7.1e-6 | < 0.0005 | 3.1e-3 |
| 12 | 2.94e5 + 2.34e6 i | 0.0568 (0.0508) | 0.991 (0.999) | 6.7e-6 / 7.8e-6 / 7.1e-6 | < 0.0005 | 2.5e-3 |

The coupled run reproduces the exact periodic state to 1e-5 in amplitude
and better than 0.0005 degrees in phase at every Womersley number, and the
mid-length velocity profile matches Womersley's to 0.25-0.31 % of its peak
(outside the wall cell: values sampled inside it converge at first order
only, see "Order in space"; the Stokes layer at alpha = 12 is 1.2 mm
thick on a wall cell of 0.05 mm). Viscosity
is not a small effect here: at alpha = 3 the pipe's resistance is as large
as its reactance, it reflects 7.9 % of the wave instead of 5.1 % and
transmits 95 % instead of 99.9 %, and the runs follow the viscous values,
not the inviscid ones. The entrance region at the inlet port leaves no
measurable trace at this length (l/R = 10).

**The wedge's geometry.** A one-cell wedge of angle theta with its flat
outer face at the distance R_c from the axis has, at every distance x,
faces and volumes of width 2 x tan(theta/2) instead of the sector's
x theta: a uniform factor k = 2 tan(theta/2)/theta = 1.000635, whatever
the radial grading. The finite-volume problem is therefore the circular
pipe of radius R_c with every flux times k. The usual wedge, with its edge
vertices at R, is the pipe of radius R cos(theta/2) (and not the pipe of
equal area, which has the same inertance but a viscous resistance k times
smaller); measured against the pipe of radius R its ports floor at 0.7 to
1.5e-4 on every mesh, the geometry and not the coupling. This wedge has
its face at R_c = R (`base/system/blockMeshDict`) and `flowScale`
2 pi/(k theta) = pi/tan(theta/2) on the ports (`makeCases.py`), so the
coupled problem is exactly the pipe of radius R and `analyse.py` compares
with Womersley's solution for R. The same holds for any wedge or polygonal
mesh: the mesh must pose the problem whose exact solution is used (as the
rigid-segment benchmark's 32-gon is scaled to the area pi R^2).

What this adds to `rigidSegment3D`: the coupling is exact for a 3D side
whose velocity profile is not uniform, i.e. the exchange of mean pressure
and flow at the ports loses nothing when the 3D flow is a genuine viscous
solution.

**The exact transient and the reference files.** The linear problem has an
exact solution from rest as well: the incident wave (delayed by L1/c0)
convolved with the impulse responses of 1 + Gamma, T/Zc and T at every
frequency (`analyse.exact_transient`, by FFT, Womersley's impedance at
each frequency with the exponentially scaled Bessel functions), and the
mid-length velocity profile from the flow's spectrum and Womersley's
profile per frequency (`exact_transient_profile`); both agree with the
periodic solution to 1e-12 at late times. `analyse.py --csv` writes them
for alpha = 3, 6 and 12 to `reference/` (ports over the first 6 s, the
periodic complex amplitudes, the profile at eight instants of the first
second), for other codes; `spaceStudy.py` measures the wedge's order in
space against them, and `physiological.py` makes a numerical reference at
5 kPa, where the 1D equations are nonlinear.

**Physiological amplitude (numerical reference).** `physiological.py run`
repeats the alpha = 6 case with an incident wave of 5 kPa (0.12 of
rho c0^2; the ports differ from the linear solution scaled to 5 kPa by
6 % in P_in, 11 % in Q and 5 % in P_out, RMS) over the first second from
rest, on the nonlinear vessels (kExact 3 + RK4, 800 cells), BDF 4 and
coupling order 4 at dt = 0.25 ms, on the 50 x 100 and 100 x 200 wedges.
The finer run is the reference, written to
`reference/alpha6_physiological_ports.csv` (every 0.25 ms) and
`_profiles.csv` (u(r, t) at nine instants), labelled numerical; the
coarser run differs from it by 5.1e-6, 5.7e-6 and 5.7e-6 (RMS / RMS of
P_in, Q, P_out) and 6.1e-4 on the mid-length profile (`physiological.md`;
outside the coarser run's wall cell, and the profile file stops at the
finer run's wall cell),
the same spacing as the linearised study above, so the reference is
converged to that. The two runs take an hour on two cores.

## Order in space

`spaceStudy.py run [--alpha 6]` (`space_alpha6.md`, `plots/space_alpha6.png`;
also alpha 3 and 12 on the two cheaper levels): the wedge at 25 x 50,
50 x 100 and 100 x 200 over the first second from rest, BDF 4 and coupling
order 4 at dt = 0.5 ms (the time error is negligible: dt = 0.25 and
0.125 ms give the same errors to three digits), linearised vessels at
fourth order, against the exact transient of the pipe of radius R; ports
(RMS / RMS of the exact) and the mid-length profile (max over eight
instants):

| family, nr x nz | ports (P, Q) | order | profile | order | PIMPLE it./step | time [s] |
|---|---|---|---|---|---|---|
| graded, alpha 6, 25 x 50 | 3.10e-5 | | 3.05e-3 | | 9.2 | 137 |
| 50 x 100 | 7.78e-6 | 1.99 | 7.54e-4 | 2.02 | 19.1 | 959 |
| 100 x 200 | 1.95e-6 | 2.00 | 1.87e-4 | 2.01 | 46.4 | 13342 |
| uniform, alpha 6, 25 x 50 | 1.37e-4 | | 4.84e-3 | | 5.8 | 84 |
| 50 x 100 | 3.43e-5 | 2.00 | 1.21e-3 | 2.00 | 8.3 | 394 |
| 100 x 200 | 8.58e-6 | 2.00 | 2.99e-4 | 2.01 | 14.8 | 4398 |
| graded, alpha 3, 25 x 50 | 5.80e-5 | | 7.69e-3 | | 16.1 | 237 |
| 50 x 100 | 1.43e-5 | 2.02 | 2.19e-3 | 1.81 | 43.0 | 2145 |
| graded, alpha 12, 25 x 50 | 2.66e-5 | | 3.84e-3 | | 6.5 | 93 |
| 50 x 100 | 6.87e-6 | 1.96 | 9.66e-4 | 1.99 | 10.5 | 507 |

The uniform family is `--family uniform` (`space_alpha6_u.md`), no wall
grading, the face at R and the same `flowScale`. The ports (P and Q agree
to the digits shown) converge at second order on both families and at
every Womersley number, and so does the profile from the axis to the wall
cell. Inside the wall cell the sampled values (`sets`, cellPoint: linear
between the cell centre and the wall's point values) converge at first
order only: on the uniform family, whose wall cell is R/nr, the profile up
to the wall cell's centre gives 1.90e-2, 1.12e-2 and 6.53e-3 (orders 0.76,
0.78), with the maximum in the wall cell; the graded family's wall cell is
ten times thinner and hides it. That is the sampling, not the solution,
and the profile errors stop at the wall cell. The time column is the
serial wall time on xenosim; the PIMPLE outer loop needs more iterations
on the finer wall-resolved meshes, a cost and not an accuracy limit.

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

**Against the exact solution** (`../orderStudy.py --case womersleySegment
--exact`; `orders_linear.md`, `plots/orders_linear.png`): the same runs
with the vessels linearised at fourth order (kExact 3, RK4, 800 cells),
each compared with the exact transient of the pipe of radius R (RMS
error / RMS of the signal, port flow; the four port signals agree to the
digits shown):

| dt [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 4 | 1.08e-3 | 2.62e-5 | 8.78e-6 | 7.95e-6 |
| 2 | 5.42e-4 | 6.60e-6 | 7.11e-6 | 7.12e-6 |
| 1 | 2.72e-4 | 6.15e-6 | 7.10e-6 | 7.10e-6 |
| 0.5 | 1.36e-4 | 6.83e-6 | 7.10e-6 | 7.10e-6 |
| orders | 0.99, 1.00, 1.00 | 1.99, 0.10, -0.15 | 0.30, 0, 0 | 0.16, 0, 0 |

Euler converges at its order to the finest step and backward at 1.99
before the floor; BDF 3 and BDF 4 are at the floor from the coarsest step.
The floor, 7.1e-6, is the same for every family and is the 50 x 100
wedge's spatial error (`spaceStudy.py`: 7.8e-6 for this mesh over 1 s), so
the time orders are those of the self-convergence above, down to the
mesh's own error.
