# Results for the 3D-1D coupling paper

Prepared by David Kelly (with Claude Code), started 2026-11-13. The working
summary of the test-case programme for the paper: for each numbered item of
the plan, the result, the figure, the commands that reproduce it, and
whether it was verified by running the case or inferred. Updated as the
work goes.

## What changed in the paper's claims

- **Dobroserdova case, the 3D run against the exact solution (corrected).**
  The report compared the pimpleFoam hard-coupling run with the delay
  equation using the pipe inductance `rho l/(pi R^2)`; the 3D pipe's
  inductance is `rho l/A3D` with `A3D` the 32-sided patch area (0.9936
  pi R^2), as the rigid-segment benchmark already states. The 0.64 % error
  in L shifts the resonance, and over 21 of its periods the RMS flow error
  read 3.8e-2 of the inflow (sine) and 1.7e-2 (step). With the right area:
  7.7e-3 and 9.3e-3 (RMS), 2.1e-2 and 3.9e-2 (max); the peak flows 1.596
  against 1.596 and 2.263 against 2.265. `dobroserdova2016/analyse.py`
  now takes the pipe area as a parameter of `exact()` and `results.md` is
  regenerated. The conclusion is unchanged (the 3D pipe is the inductor,
  the remaining error is the BDF2 phase error of the resonance), the
  numbers are five times smaller.
- **Reflection sign of the flow inlet (wording).** Γ = +1 in pressure, not
  −1 (item 10).
- **Womersley segment: the agreement is ten times better than reported.**
  The wedge mesh has a flat wall, so its cross-section is
  0.5 R^2 sin(theta) per wedge, 0.99873 of its share of pi R^2, and the 3D
  problem solved is the circular pipe of that area (R_eq = R sqrt(sin
  theta/theta) = 9.99365 mm). Against the pipe of radius R the report's
  amplitude errors (0.7 to 1.5e-4) and phases (0.003 to 0.004 degrees)
  are this geometry, not the coupling or the mesh: against the equal-area
  pipe they are 3e-6 (alpha 6 and 12), 1.2e-5 (alpha 3) and 0.0002 to
  0.0004 degrees. `results.md` now carries both rows, and the spatial
  study (item 3) measures against the equal-area pipe: against pi R^2 it
  floors at 6e-5 whatever the mesh or time step (verified: 25 x 50 and
  50 x 100 at dt = 0.5, 0.25 and 0.125 ms all give 5.9e-5).
- **The aorta floor is the 3D mesh, not open (item 4).** The report said
  the cause was open with the stubs and junctions as candidates; it is the
  non-orthogonal cone mesh in pimpleFoam (the cylinder converges at 3.8 to
  3.9).

## Preliminaries

**Input files.** The report
(`~/OpenFOAM/davidkelly-v2412/circulation3D1D_report/circulation3D1D_report.pdf`,
of 9 November 2026) was used as the baseline. The paper plan
(`paper_plan_3D-1D_coupling.md`) and the literature map
(`3D-1D_coupling_literature_map.md`) were not found anywhere on this
machine (searched the home directory, `~/Papers/3D - 1D Coupling` and the
checkout); the work below does not need them until the paper draft.

**Paths quoted in the plan, checked.** `src/circulationModels/circulation1D`,
`src/circulationModels/circulation1DCoupling` and
`src/circulationModels/docs/benchmarks3D1D.py` exist. The cases are under
`tutorials/circulationModels/coupling3D1D/` as named. Not every case has an
`analyse.py`: `rigidSegment3D` has `sweep.py`, `analyse.py` and
`highOrder.py`; `dobroserdova2016`, `womersleySegment` and `thoracicAorta37`
have `makeCases.py` and `analyse.py`; `formaggia2007` has `makeCases.py`,
`analyse.py` and `orderStudy.py`; the common time-order study is
`coupling3D1D/orderStudy.py`.

**Baseline reproduced (verified).** `rigidSegment3D`, implicit order 2 at
dt = 0.5 ms (`python3 sweep.py --one`), rerun on 2026-11-13: the 3D flux and
both pipe-end pressures agree with the stored run to the last printed digit
(maximum difference 0 at every step), and `analyse.py` regenerates
`results.md` unchanged.

**State of the working tree.** The checkout had uncommitted changes from
another session when this work started (the 0D coupler
`circulationCoupling`, `circulation0D` components, a new
`circulationMeshMotion` library, a `leftVentricle` 0D tutorial, the BDF
scheme of PR #674). They are left untouched and are not part of the commits
made here.

## Phase 1

### 1. Fourth order against the exact solution (rigid segment, B1)

**Linearised 1D mode (implemented).** A uniform vessel with
`linearised true` solves the 1D equations linearised about its initial
state: the tube law becomes its tangent at `Aref` (`linearisedTubeLaw`,
constant `c0`), the convective term, the dynamic pressure at every node and
the nonlinear invariants are replaced by their linear forms
(`U = Q/Aref`, `W = Q/Aref +- c0 (A - Aref)/Aref`), so that every node
imposes its condition on the linear characteristics. The nonlinear path is
unchanged: four tutorial cases rerun with the rebuilt solver reproduce
their stored results to round-off (rigid segment and junction cases to
0; carotid to 1.6e-14 in the physical columns).

**1D only, against the closed form (verified).** `rigidSegmentInertance`,
the pipe as an `inertance` node, linearised vessels at full amplitude
(incident peak 20.6 Pa, inflow x 1), kExact degree 3 + RK4 at CFL 0.5,
errors of P along both vessels at t = 0.6 s against the cell averages of
the closed form (`./Allrun`; `results.md`, `plots/rigidSegment.png`):

| cells per vessel | dx [mm] | L1 | order | L2 | order | max | order |
|---|---|---|---|---|---|---|---|
| 200 | 10 | 1.165e-03 |  | 3.747e-03 |  | 2.440e-02 |  |
| 400 | 5 | 8.173e-05 | 3.83 | 2.797e-04 | 3.74 | 2.013e-03 | 3.60 |
| 800 | 2.5 | 5.110e-06 | 4.00 | 1.757e-05 | 3.99 | 1.286e-04 | 3.97 |
| 1600 | 1.25 | 3.173e-07 | 4.01 | 1.091e-06 | 4.01 | 7.998e-06 | 4.01 |
| 3200 | 0.625 | 1.981e-08 | 4.00 | 6.793e-08 | 4.01 | 4.979e-07 | 4.01 |
| 6400 | 0.3125 | 1.387e-09 | 3.84 | 4.239e-09 | 4.00 | 3.119e-08 | 4.00 |

Fourth order on every norm down to 1e-9 of the peak, with no floor. Two
things learned on the way, both of which a benchmark specification must
state: (i) the comparison must be of cell averages with cell averages
(against point values of the closed form at the cell centres the same runs
show a second-order floor of (dx/sigma_x)^2/24, e.g. 2.6e-5 at 1.25 mm);
(ii) the linearised run at the base resolution (1600 cells, SSPRK2) at full
amplitude gives exactly the errors of the nonlinear run at 1e-3 amplitude
(4.67e-4 / 1.51e-3 / 1.02e-2), as it should.

**3D-1D, against the closed form (verified).** `rigidSegment3D/highOrder.py
run 6 --linear` (20 cases, about 55 minutes with 6 at a time; `highOrder_linear.md`,
`plots/highOrder3D1D_linear.png`): linearised vessels at full amplitude,
kExact 3 + RK4, the pipe with `BDF 4`, dx and dt halved together over five
levels, errors against the closed form with `tau` from the 32-sided patch
area, relative to the peaks; the profiles against the closed form's cell
averages.

| coupling | level (dx, dt) | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| implicit 4 | 2.5 mm, 2 ms | 1.43e-03 | | 1.10e-03 | | 1.75e-03 | |
| | 1.25 mm, 1 ms | 9.41e-05 | 3.93 | 7.30e-05 | 3.92 | 1.14e-04 | 3.95 |
| | 0.625 mm, 0.5 ms | 5.99e-06 | 3.97 | 4.62e-06 | 3.98 | 7.14e-06 | 3.99 |
| | 0.3125 mm, 0.25 ms | 3.77e-07 | 3.99 | 2.89e-07 | 4.00 | 4.46e-07 | 4.00 |
| | 0.15625 mm, 0.125 ms | 2.36e-08 | 4.00 | 1.80e-08 | 4.00 | 2.82e-08 | 3.98 |
| explicit 4 | 2.5 mm, 2 ms | 4.11e-03 | | 4.16e-03 | | 2.96e-03 | |
| | 1.25 mm, 1 ms | 2.78e-04 | 3.88 | 2.80e-04 | 3.89 | 1.90e-04 | 3.96 |
| | 0.625 mm, 0.5 ms | 1.78e-05 | 3.97 | 1.77e-05 | 3.99 | 1.19e-05 | 4.00 |
| | 0.3125 mm, 0.25 ms | 1.12e-06 | 3.99 | 1.11e-06 | 4.00 | 7.42e-07 | 4.00 |
| | 0.15625 mm, 0.125 ms | 7.05e-08 | 3.99 | 6.91e-08 | 4.00 | 4.65e-08 | 4.00 |
| implicit 2 | five levels | as implicit 4 | 3.93 to 4.00 | as implicit 4 | 3.92 to 4.00 | 2.69e-03 to 1.52e-05 | 1.61, 1.89, 1.97, 1.99 |
| implicit 1 | five levels | as implicit 4 | 3.93 to 4.00 | as implicit 4 | 3.92 to 4.00 | 6.14e-02 to 3.90e-03 | 0.98, 1.00, 1.00, 1.00 |

The claim of the report holds against the exact solution, with no floor
and no amplitude caveat: fourth order in space and time with fourth-order
coupling, implicit or explicit, on the 3D flux, the pipe-end pressures and
the waves in the vessels; coupling of order p gives order p on the waves
the ports send into the vessels. One difference from the self-convergence
study of the report (nonlinear vessels at full amplitude): there, implicit
order 1's first-order wave error "leaked" into the 3D flux at the finest
level (order 2.43, 1.46e-6 against 4.4e-7); against the exact solution of
the linear problem the 3D flux of implicit 1 is at the fourth-order error
at every level (3.77e-7 at 0.3125 mm, identical to implicit 4's). The
first-order error is in the waves the ports emit, which in the linear
problem travel away from the pipe and never return within the window; in
the nonlinear problem the erroneous reflected wave interacts with the
incident one at second order in amplitude while they overlap at the pipe,
which is the leak. So the floor of 2e-5 the plan mentions at full
amplitude is the nonlinearity of the 1D equations, not of the coupling.

No order reduction from the time-dependent boundary data at the
Runge-Kutta stages was found: the coupler passes the port flow as a cubic
in time over each step and the 1D RK4 stages evaluate it at the stage
times, and the orders above are clean to 2e-8.

### 2. Time order against the exact solutions

**Dobroserdova resonator (verified).** `orderStudy.py --case dobroserdova2016
--exact` (20 runs, 12 minutes with 3 at a time; `dobroserdova2016/orders_linear.md`,
`plots/orders_linear.png`): the hard-coupling 3D run at dt = 2 to 0.125 ms
with the vessel linearised at fourth order in space and time (kExact 3,
RK4, 200 cells), every run against the delay-equation solution (integrated
with the pipe's patch area, 64000 steps per round trip, accurate to 1e-8).
RMS error / RMS of the signal over the first period, port flow:

| dt [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 2 | 1.36e-1 | 1.06e-2 | 1.48e-3 | 1.28e-3 |
| 1 | 8.55e-2 | 2.42e-3 | 2.19e-4 | 9.63e-5 |
| 0.5 | 4.64e-2 | 5.51e-4 | 2.91e-5 | 1.22e-5 |
| 0.25 | 2.38e-2 | 1.38e-4 | 4.48e-6 | 2.15e-6 |
| 0.125 | 1.21e-2 | 3.44e-5 | 2.21e-6 | 2.26e-6 |
| orders | 0.67, 0.88, 0.96, 0.98 | 2.13, 2.14, 2.00, 2.00 | 2.76, 2.91, 2.70, 1.02 | 3.74, 2.98, 2.50, -0.07 |

Euler and backward reach their orders against the exact solution; BDF 3
and BDF 4 run at 2.9 and 3.7 until a floor of 2.2e-6 in the flow; the
pressure has a floor of 5.3e-4 for every family. Both floors are the 1D
vessel's representation of the exact solution's discontinuities, not the
coupling: the sine inflow starts from rest with a finite slope, so the
pipe pressure `L dQ/dt` jumps at every round-trip arrival of the front
and the flow has a kink; the 1D-only linearised run with the pipe as an
inertance node, same scheme and cells, has exactly these floors (2.33e-6
and 5.29e-4) and they fall at first order in dx for P and second for Q
(verified over 100 to 800 cells). So the exact comparison confirms the
self-convergence orders down to the level the problem's own
non-smoothness allows, and a referee should be told that the resonator
test's exact solution is piecewise smooth.

**Womersley segment (verified).** Exact transient solution implemented
(`womersleySegment/analyse.py: exact_transient`, the Womersley transfer
function at every frequency convolved with the ramped incident wave by
FFT; checked against the periodic Womersley solution to 1e-12 at late
times; its own grid error is 1e-11). The reference is the circular pipe
of the wedge's cross-section (see "what changed"). `orderStudy.py --case
womersleySegment --exact` (16 runs, about 1.5 h with 3 at a time;
`womersleySegment/orders_linear.md`, `plots/orders_linear.png`):
linearised vessels at fourth order, alpha = 6, 0 to 0.6 s. RMS error /
RMS of the signal, port flow (the four port signals agree to the digits
shown):

| dt [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 4 | 1.08e-3 | 3.09e-5 | 7.34e-6 | 6.55e-6 |
| 2 | 5.48e-4 | 9.86e-6 | 5.57e-6 | 5.61e-6 |
| 1 | 2.77e-4 | 6.06e-6 | 5.57e-6 | 5.58e-6 |
| 0.5 | 1.41e-4 | 5.63e-6 | 5.57e-6 | 5.57e-6 |
| orders | 0.99, 0.99, 0.97 | 1.65, 0.70, 0.10 | 0.40, 0, 0 | 0.22, 0.01, 0 |
| self-convergence, finest triple | 1.00 | 1.99 | 2.86 | 3.32 |

Euler converges at its order to the finest step; backward at 1.65 before
the floor; BDF 3 and BDF 4 are at the floor from the coarsest step. The
floor, 5.6e-6, is the 50 x 100 wedge's spatial error (item 3 gives 4.3e-6
for this mesh over a 1 s window) and is the same for every family, so the
exact comparison is consistent with the self-convergence orders down to
it. Against the pipe of radius R instead of the equal-area pipe the floor
is 6e-5 (the earlier, wrongly referenced comparison, kept in the git
history of this file).

### 3. A 3D spatial order (Womersley segment, B2)

`womersleySegment/spaceStudy.py run` (`space_alpha6.md`,
`plots/space_alpha6.png`): the wedge mesh at 25 x 50, 50 x 100 and
100 x 200 (radial x axial, the radial grading 10:1 kept), alpha = 6, 0 to
1 s from rest, BDF 4 and coupling order 4 at dt = 0.5 ms, linearised
vessels at fourth order, against the exact transient of the equal-area
pipe: port pressures and flow (RMS / RMS of the exact), mid-length
velocity profile (max over eight instants of max |u - u_exact| / max
|u_exact|). The time error is negligible: the 25 x 50 level at dt = 0.25
and 0.125 ms gives the same errors to three digits (verified).

Results (verified; `space_alpha{3,6,12}.md`, `plots/space_alpha*.png`):

| alpha | mesh | ports (P, Q) | order | mid-length profile | order | PIMPLE it./step (steps at the cap of 50) | wall time |
|---|---|---|---|---|---|---|---|
| 3 | 25 x 50 | 3.43e-5 | | 7.68e-3 | | 15.6 (4 of 2000) | 139 s |
| 3 | 50 x 100 | 1.35e-5 | 1.35 | 2.18e-3 | 1.82 | 42.1 (240) | 1531 s |
| 6 | 25 x 50 | 2.47e-5 | | 5.32e-3 | | 8.9 (20) | 78 s |
| 6 | 50 x 100 | 4.61e-6 | 2.42 | 1.17e-3 | 2.19 | 18.6 (34) | 681 s |
| 6 | 100 x 200 | 6.49e-6 | -0.49 | 1.36e-3 | -0.23 | 45.4 (1558) | 12276 s |
| 12 | 25 x 50 | 2.48e-5 | | 1.10e-2 | | 6.2 (52) | 57 s |
| 12 | 50 x 100 | 5.51e-6 | 2.17 | 2.50e-3 | 2.14 | 10.2 | 391 s |

Between the 25 x 50 and 50 x 100 meshes the ports converge at 2.2 to 2.4
and the profile at 2.1 to 2.2 at alpha = 6 and 12: second order, as the
schemes are. The third level at alpha = 6 does not continue: its outer
loop runs to the 50-iteration cap on 1558 of its 2000 steps (45
iterations per step; the pressure residual stalls near 1e-7 instead of
reaching 1e-10) and the error rises to 6.5e-6. The same stall caps alpha
= 3 on the 50 x 100 mesh (240 steps at the cap, order 1.35). So the
achievable error of this case is limited by the convergence of the
segregated algorithm without momentum predictor on a wall-graded mesh
(the report's note on alpha = 3), not by the mesh, from about 5e-6 of
the signal on; the second order of the 3D discretisation is established
on the two coarser levels at the two higher Womersley numbers. A
referee should see the stall stated as such; curing it (a momentum
predictor with relaxation, or a coupled solver) is outside this item.
Against the pipe of radius R instead of the equal-area pipe the alpha =
6 port errors are 5.9e-5 and 8.0e-5 on the first two levels, i.e. no
convergence, which is how the wedge-area effect was found.

### 4. The floor in the network case (B6): traced (verified)

`thoracicAorta37/floorStudy.py` (15 BDF 4 runs, 3 at a time, about 50
minutes; `floor.md`, `plots/floor.png`): the BDF 4 runs of the order
study at dt = 1, 0.5 and 0.25 ms with one change each. Differences
between successive steps of the port flow (RMS / RMS of the signal) and
the order they give:

| variant | e(1, 0.5) | e(0.5, 0.25) | order (Q / P) |
|---|---|---|---|
| baseline: cone, tapered stubs | 2.18e-05 | 1.88e-05 | 0.21 / 0.16 |
| linear solvers 1e-15, p residual 1e-12 | 2.18e-05 | 1.88e-05 | 0.21 / 0.16 (runs identical to 3e-13) |
| one non-orthogonal corrector | 2.18e-05 | 1.88e-05 | 0.21 / 0.16 (identical to 7e-13) |
| uniform stubs, cone unchanged | 2.18e-05 | 1.89e-05 | 0.21 / 0.16 |
| `ddtCorr no` | 1.48e-05 | 1.04e-05 | 0.50 / 0.70 |
| cylinder (mean radius) instead of the cone | 2.76e-06 | 1.95e-07 | 3.82 / 3.91 |

The floor is the 3D discretisation on the non-orthogonal cone mesh. The
same network, stubs, junctions, coupling and schemes with the 3D segment
as a cylinder converge at 3.8-3.9, as the uniform-vessel cases do. Not
involved: the tapered stubs' end treatment, the junctions, the solver
and PIMPLE tolerances, output precision (the floor is in coupling.dat's
12-digit columns as much as in the 6-digit function-object output), the
1D sub-stepping (ruled out before). On the cone, pimpleFoam's flux
correction `fvc::ddtCorr` (its coupling coefficient switches on the local
mismatch between the interpolated velocity and the flux, which vanishes on
an orthogonal mesh in plug or developed flow) carries about half; the
rest is the non-orthogonal treatment in the segregated algorithm. Both
sit below the second-order error (`backward` converges cleanly at 1.95)
and above the fourth-order one.

What to claim: the coupling takes the 3D scheme's order on a 37-vessel
network when the 3D discretisation has that order (the cylinder, 3.8-3.9),
and a non-orthogonal pimpleFoam mesh does not go beyond second order in
time whatever the ddt scheme, which is a statement about the 3D solver,
not the coupling. The bound for the cone: 3e-6 of the signal, dt
independent.

### 5. Static against total pressure at the ports (B6)

Implemented in `circulation1DCoupler` as `coupling { pressure total; }`
(default `static`). The patch gets the uniform value
`p = P1D + rho U1D^2/2 - pdyn3D`, with U1D = Q/A1D the velocity of the
vessel's boundary state (the 1D model's own total pressure, as at its
junctions; no dynamic part for a linearised vessel) and pdyn3D the
patch's dynamic pressure from the sampled flux, by `kineticEnergy`:
`mean` (rho U3D^2/2, U3D = Q/A3D; the usual mean total pressure, the
default), `area` (rho <u_n^2>/2, the area average of the local dynamic
pressure, 4/3 of the mean form for Poiseuille flow) or `flux` (the
energy-flux form rho <u_n^3>/(2 <u_n>), 2 times the mean form for
Poiseuille flow, which conserves the mechanical energy flux across the
interface). The 1D limit (the `inertance` node) has its own
`pressure total` with the mean velocities.

**Runs (verified).** The tutorial's coupled run (restart from the periodic
state, two cycles, dt = 1 ms, implicit order 2) with each form, and the
inertance limit with total pressure (`makeCases.py --restart --pressure
total --kineticEnergy mean|area|flux`; `results.md`). Differences of the
total-pressure runs from the static one over the last cycle, relative to
the pulse at the site (4.5 to 6.2 kPa) and the peak flow:

| site | mean: RMS dP / max dP / RMS dQ | area | flux | 1D inertance limit, total vs static |
|---|---|---|---|---|
| aortic arch II | 6.0e-4 / 1.9e-3 / 8.5e-4 | 6.1e-4 / 1.9e-3 / 8.8e-4 | 6.3e-4 / 2.1e-3 / 1.0e-3 | 6.3e-4 / 2.1e-3 / 1.0e-3 |
| thoracic aorta II (3D mid-plane) | 2.9e-3 / 8.9e-3 / 1.6e-3 | 3.3e-3 / 9.4e-3 / 1.6e-3 | 5.7e-3 / 1.0e-1 / 1.8e-3 | 5.5e-4 / 1.6e-3 / 1.8e-3 |
| right iliac-femoral II | 1.3e-3 / 3.9e-3 / 2.0e-3 | 1.3e-3 / 4.0e-3 / 2.1e-3 | 1.5e-3 / 4.3e-3 / 2.4e-3 | 1.4e-3 / 4.7e-3 / 2.4e-3 |

In pascals at the 3D mid-plane: the mean form lowers the pressure by
42 Pa at peak flow (120.8 ml/s, U = 0.38 m/s: rho U^2/2 = 75 Pa at the
stub ends, of which the two ports' corrections leave this net), RMS 14 Pa;
the area form 45 Pa; the flux form 49 Pa at peak flow and the 475 Pa
spike at reversal. At the arch the three forms are within 8 to 9 Pa of
the static run (RMS 2.7 Pa, 0.06 % of the pulse), at the iliac within
24 to 26 Pa. So at physiological amplitude in a vessel of this size the
choice of static or total pressure at the ports moves the network by
0.06 to 0.2 % of the pulse (RMS) and the 3D segment's own pressure by
0.3 % (RMS, 0.9 % at peak flow); the mean and area forms are
interchangeable to 5 % of that, and the energy-flux form is not usable
with reversal. The coupled runs with total pressure sit as close to the
total-pressure inertance limit as the static one to the static limit
(1.0e-3 against 2.0e-3 at the arch, 2.6e-3 against 4.7e-3 at the iliac,
RMS dP / pulse): the limit shifts with the condition and the 3D segment
follows it.

**Finding on the energy-flux form (verified).** It is singular at flow
reversal: `<u_n^3>` stays finite while `<u_n>` passes through zero, so the
imposed pressure spikes. In the aorta run it reaches 475 Pa (10 % of the
pulse) at Q = 0.2 ml/s and 277 Pa at 0.1 ml/s, two samples per cycle;
away from reversal the correction is -49 Pa at peak flow at the 3D
mid-plane (1 % of the pulse) and 9 Pa at the aortic arch (0.2 %). That
form should not be used on a pulsatile case with reversal, which is why
`mean` is the default.

### 6. Cost (verified from the logs)

Serial unless stated; "1D sub-steps" is the number of CFL steps the
network takes per 3D step for one integration (`ceil(dt / (CFL
min(dx/c)))`), and "per step" multiplies it by the network
re-integrations the implicit coupling makes in that step (one per PIMPLE
outer iteration, two per iteration on the Womersley case where the
pressure is evaluated twice). Wall time per step on this machine (Apple
silicon, 14 cores).

| case | run | 3D cells | dt | outer it./step | 1D sub-steps (per step) | wall time | ranks | ms/step |
|---|---|---|---|---|---|---|---|---|
| rigid segment | implicit 2, dt 0.5 ms, 0.35 s | 3072 | 0.5 ms | 11.1 | 5 (55) | 124 s | 1 | 177 |
| rigid segment | implicit 4, 6400 cells, dt 0.125 ms (finest of the fourth-order study) | 3072 | 0.125 ms | 6.8 | 10 (68) | 742 s | 1 | 265 |
| Womersley segment | alpha 6, 6 s, BDF2 | 5000 (wedge) | 1 ms | 21.3 | 5 (213) | 2532 s | 1 | 422 |
| FSI tube (`formaggia2007`) | coupled, Robin-Neumann, 20 ms | 3dTube quarter | 25 us | 3 FSI it./step (40 at dt 100 us) | 1 (the 1D CFL step, 59 us, is longer than dt) | 685 s | 1 | 857 |
| FSI tube | full3D reference (15 cm), 20 ms | 3 x | 25 us | 3 | - | 867 s | 4 | 1084 |
| 37-artery network | rigid3D, two cycles (1.654 s) | 27670 (cone) | 1 ms | 5.1 (10.1 re-integrations) | 14 (142) | 273 s | 1 | 165 |

The 1D network is never the cost: 37 vessels at 2 mm cells take about
0.5 ms per CFL step, so even 142 sub-steps per 3D step are a fifth of the
165 ms. The implicit coupling's price is the outer iterations it needs,
5 on the network at physiological amplitude (the 3D segment is short
against the vessels, dt/tau small), 11 on the rigid segment at dt = 0.5 ms
(dt/tau = 0.077) and 21 on the viscous wedge, where the segregated
algorithm without momentum predictor on the wall-graded mesh sets the
count, not the coupling (the uncoupled pipe needs as many).

## Phase 2

### 8. The compliant FSI tube (B4): in progress

**(a) Elastic support on the cut end of the wall.** Implemented in
`formaggia2007/makeCases.py --springEnd` with the `solidSpringDashpot`
condition of solids4foam on the wall's cut face: axially free (`kNormal`
0) and, in the face's plane (radial and circumferential), the stiffness of
the continuing wall taken as the edge stiffness of a semi-infinite thin
cylindrical shell displaced radially with its edge moment free
(Timoshenko): shear force `Q = 2 D beta^3 u` per unit circumference,
`D = E h^3/(12 (1 - nu^2))`, `beta^4 = 3 (1 - nu^2)/(R^2 h^2)`, divided
by the wall thickness for a traction per unit cut area: `kTangential` =
1.04e7 Pa/m for this wall (`shell_edge_stiffness()`). This is the static
stiffness the continuing tube presents to its cut; it sits between
clamped (infinite) and free (zero). Note on the machine: the installed
`libsolids4FoamModels` (February 2026) has neither this condition nor the
BDF schemes, so the FSI runs of this item use a full build of the
checkout's solids4FoamModels in a second scratch directory; the clamped
case is rerun with it first to show the library change itself does not
move the results (pending).

**Library caveat (verified, needs a decision).** The clamped coupled run
`coupled_c4p92` rerun with the checkout's own solids4FoamModels
(`coupled_c4p92_lib2`, same case files) differs from the run made with
the installed February library by 7.4 % of the peak in the port pressure
(RMS over the 20 ms), 2.6 % in the pressure peak at z = 12.5 mm (1352
against 1388 Pa) and 2.3 % in the peak wall displacement (146.8 against
150.2 um). The case is unchanged; the solids4foam FSI machinery between
the two versions is not (nine months of development on the Robin-Neumann
coupling and the solid solver, plus this checkout's uncommitted BDF
moving-mesh files from another session). The full tube is being rerun
with the same build (`full3D_lib2`) so that every end-condition
comparison of this item is made within one library; the report's FSI
numbers (section on Formaggia) belong to the installed library. Which
version the paper should use is a decision for the author; the
coupled-against-full comparison is what matters and is being redone
consistently.

**(a) Result (verified; all four runs with the checkout build,
`results.md`).** Peaks and RMS differences against `full3D_lib2`:

| z [cm] | full tube [Pa] | clamped | free | spring | RMS(P - full) / pulse: clamped / free / spring |
|---|---|---|---|---|---|
| 2.5 | 1271 | 1303 | 1233 | 1233 | 3.0e-2 / 4.4e-2 / 5.1e-2 |
| 3.75 | 1132 | 1169 | 1108 | 1111 | 4.6e-2 / 2.5e-2 / 2.5e-2 |
| 5 (interface) | 1077 | 1104 | 1093 | 1098 | 4.4e-2 / 1.8e-2 / 2.2e-2 |
| 6.25 | 1045 | 1096 | 1085 | 1089 | 6.9e-2 / 3.6e-2 / 4.6e-2 |

The shell's static edge stiffness (1.04e7 Pa/m) is about a sixth of the
ring's hoop stiffness E h/R^2, so the elastic support behaves nearly as
the free end; the full tube lies between clamped and free/spring at every
section. The static stiffness of the continuing wall is not what holds
the ring in the full tube (its inertia and the pressure on the
continuation are); the condition to build next is the continuity of area
proper, a radial displacement of the cut ring driven by the 1D area,
which needs the coupler to drive a solid boundary (not built; it is
method development of the kind Phase 3 asks about first).

**(b) The wave returning from the 1D side.** `--closedEnd` writes
`full3D_closed` (the 15 cm FSI tube with a wall at its far end) and
`coupled_c4p92_closed` (the 1D vessel's far end a zero-flow inlet, which
reflects with +1 in pressure), so the pulse comes back through the
interface from 1D into 3D; `analyse.py` compares a `_closed` coupled run
with `full3D_closed`. Cases written, runs pending (about 30 minutes: the
coupled run serial, the full tube on 4 ranks).

**(b) Result (verified; 70 ms runs, `results.md`).** The wave returns
through the interface into the 3D part within 1-2 ms of the closed full
tube at every section (peaks at z = 1.25, 2.5, 3.75, 5 cm at 59.0 /
56.4 / 53.5 / 50.1 ms against 58.0 / 56.5 / 54.2 / 51.9 ms) and with the
1D vessel's amplitude: 844-1124 Pa against the full tube's 475-680 Pa,
which the thick wall's inertia has dispersed and attenuated on the way
out and back; the difference is already 20 % on the outgoing wave at
z = 10 cm (1083 against 902 Pa) before anything has returned. RMS
difference over 30-70 ms 0.09-0.12 of the pulse. The coupling passes a
wave from 1D into 3D as it does the other way; the amplitude mismatch is
the 1D model's, counted twice. (The 20 ms window first tried does not
reach the return; the cases now default to that end time only with
`--endTime`.)

**(c) Refinement with the smooth pulse.** Not launched: the estimate is
above the two-hour limit. Per end condition (clamped, free, spring) the
coupled run at the tutorial's mesh and dt = 25 us takes 12 minutes serial
and the next level (half the cells in each direction and half the step)
about 8 times that; the full-tube reference at the fine level about 2
hours on 4 ranks. Two levels for three end conditions plus two reference
levels is 7 to 8 hours of wall time with the machine otherwise idle. A
cheaper version that still answers Malossi's observation (does the
reflection grow or shrink with refinement, and for which end condition):
two levels for the three end conditions at dt = 25 us with the mesh alone
refined (the time error is far below the mesh error at this step, see the
order study), about 3 hours. Awaiting a decision.

### 9. The step input in the resonator: limiters (verified)

`dobroserdova2016/limiterStudy.py` (12 1D runs, 1.5 minutes; `limiters.md`,
`plots/limiters.png`): the step case 1D-only and linearised, MUSCL with
minmod, vanLeer, MC and no limiter under SSPRK2, and the unlimited
k-exact reconstructions of degree 2 (SSPRK3) and 3 (RK4), on 100 and 400
cells at CFL 0.5, against the exact delay-equation solution (a jump at
every round trip, total variation 19.1 x the inflow over 0.5 s).

| cells | minmod | vanLeer | MC | none | kExact 2 | kExact 3 |
|---|---|---|---|---|---|---|
| 100: overshoot / inflow | -0.003 | -0.001 | -0.001 | -0.000 | +0.000 | +0.000 |
| 100: TV / TV exact | 0.859 | 0.929 | 0.951 | 0.968 | 0.969 | 0.999 |
| 100: RMS error / inflow | 1.55e-2 | 1.08e-2 | 1.10e-2 | 1.17e-2 | 4.00e-3 | 2.37e-3 |
| 400: TV / TV exact | 0.953 | 0.984 | 0.992 | 1.003 | 0.996 | 1.002 |
| 400: RMS error / inflow | 4.62e-3 | 3.00e-3 | 3.13e-3 | 3.48e-3 | 8.86e-4 | 4.69e-4 |

No scheme overshoots the exact peak (the inductor integrates what the
vessel sends it). The limiters show as a deficit of total variation (the
smeared jumps: 14 % minmod, 5 % MC at 100 cells); the unlimited schemes
ring by at most 0.3 % of the total variation and have a quarter to a
sixth of the RMS error. The jumps are first order in dx for every scheme.
Documented in the case README; no new limiter.

### 10. Reflection sign of the flow inlet (done)

A prescribed flow reflects an incoming wave with Γ = +1 in pressure (−1 in
flow). The code imposes the flow on the outgoing characteristic, so it has
the right behaviour by construction, and the delay-equation script of
`dobroserdova2016/analyse.py` uses `f(t) = Zc Qin(t - td) + g(t - 2 td)`,
the +1 reflection. Only the report's wording in Section 6.2 ("Γ = −1 in
pressure") was wrong: corrected in the report source, and the sign is now
stated in `flowInlet.H` and the `circulation1D` README. Verified by reading
the code and the script; no run needed.

## Phase 3 (not started; asked for a feasibility note only)

### A 3D segment that carries waves with matched impedance: feasibility

The idea: make the 3D fluid weakly compressible with a sound speed equal
to the 1D wave speed, c0 = 6.17 m/s, so that a rigid-walled 3D segment
has the acoustic impedance rho c0/A of the vessels and transmits the pulse
without the inertance reflection, as a compliant segment of matched wall
would, but without FSI.

What exists in this code base. OpenFOAM.com v2412 has `rhoPimpleFoam`
(pressure-based, compressible, any equation of state through the thermo
library, including a linear `rhoConst`-like law `perfectFluid`
(rho = rho0 + psi p) and `rPolynomial`), and solids4foam's fluid models
include a compressible one (`compressibleFluid`-type solvers wrapping the
same libraries). The 1D coupler is a boundary condition on the kinematic
or static pressure of the patch plus a flux sample; both exist on the
compressible side (`p` is absolute there, the coupler would impose P1D
directly instead of P1D/rho, a one-line switch, and sample `phi` divided
by rho). So the coupling mechanics carry over with small changes.

What the physics needs. A linear acoustic pipe of length l with sound
speed c0 and the same area as the vessels is exactly a matched
transmission line: Gamma = 0 at both ends, delay l/c0, and the exact
solution of the coupled problem is the incident pulse delayed, at any
amplitude in the linear regime. That makes it a cleaner exact benchmark
than the rigid segment (whose exactness rests on plug flow), and one with
a 3D wave inside the segment. Two things must hold for the 3D to be a
transmission line: (i) the pressure wave must be plane, i.e. the segment
must be straight with uniform patches and the wavelength (c0 T = 6 m at
1 Hz, 4 cm for a 7 ms Gaussian's width) much longer than the radius, which
holds; (ii) the 3D time step must resolve the acoustic transit, l/c0 =
13 ms for 8 cm, which it does trivially (the acoustic CFL at c0 = 6 m/s
is benign; this is the virtue of a slow sound speed).

What stands in the way. (1) The equation of state: perfectFluid gives
c^2 = 1/psi at fixed temperature; a sound speed of 6 m/s with rho = 1050
needs psi = 0.028 s^2/m^2, far outside the thermo library's expectations
(it is written for gases and liquids near 1500 m/s) but not forbidden; the
energy equation is a nuisance (set an isothermal, constant-cp fluid and
ignore it, or use the solids4foam compressible fluid model if it allows
dropping the energy equation). (2) The pressure-based algorithm with
`transonic no` treats the compressibility term psi dp/dt implicitly and
is second order with `backward`; whether BDF 3 and 4 can be applied to
the density equation as well depends on `rhoPimpleFoam`'s use of
`fvm::ddt(psi, p)`, which goes through the same ddt scheme interface and
so should work, but the pressure-based compressible algorithm's own
accuracy in time has not been verified here. (3) The outflow boundary at
a matched impedance: the coupler imposes a uniform pressure on the patch,
which is correct for a plane wave; the velocity condition must let the
wave out (`pressureInletOutletVelocity` as now). (4) Nothing of this
exists as a test case; the rigid-segment benchmark's scripts carry over
(the closed form becomes the pure delay).

Estimate: a week of work to have the compressible rigid pipe running with
the 1D coupler and compared with the delay solution; the fourth-order
study would follow the existing highOrder.py. Worth doing if the paper
wants a 3D segment that does not reflect, i.e. a demonstration that the
coupling is exact when the 3D side itself carries the wave; it does not
replace the FSI segment, which is the physical case.

### An impedance (Robin) interface condition

Not touched, as asked. For the record, the ingredient is already in the
code: within one step each vessel end answers the port flow as a resistor
Zc (the coupler's documentation and BENCHMARKS3D1D.md), so a patch
condition `p = P* + Zc (Q - Q*)`, with P* and Q* the network's prediction,
would make the outer iteration converge in one pass for a rigid segment
(the model problem's contraction factor dt/(a0 tau) becomes zero) and
would be the natural place for an "impedance" condition. It is method
development and changes the coupler, so it waits for a decision.
