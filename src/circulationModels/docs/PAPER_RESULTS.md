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
- **The FSI tube's numbers depend on the solids4foam library (item 8).**
  The report's Formaggia results were made with the installed February
  2026 `libsolids4FoamModels`; the same case with the checkout's own build
  differs by 7 % of the peak in the port pressure and 2 to 3 % in the
  pressure and displacement peaks. Every end-condition comparison of item
  8 is made within one library (the checkout's), against its own full-tube
  reference; which library the paper quotes is the author's decision.
- **The rigid segment's "floor" at full amplitude is the 1D nonlinearity
  (item 1).** Fourth order holds against the exact solution of the linear
  problem to 2e-8 with no floor; the 2e-5 floor of the report's
  self-convergence study at full amplitude is the reflected wave of a
  low-order coupling interacting with the incident one at second order in
  amplitude, which the linear problem does not have.

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

## Housekeeping (second round, 2026-11-16)

**Code state.** Every result of the first round was made from the
working tree of `feature/circulation1D` at or before commit `726e5ec61`
(the item 7 commit; the commits of the round are listed in the git log
from `6045b6dcd` on). OpenFOAM: v2412 (`_8dbc61e11c-20241220`, OpenCFD
patch 0) built with Apple clang on Apple silicon (this Mac, 14 cores);
on xenosim the Ubuntu 22.04 package openfoam2412 (build string recorded
below), gcc 11.4.

**The BDF time scheme used by pimpleFoam is not committed.** Every
`BDF 3` and `BDF 4` result of the fluid side (items 1, 2, 3, 4, 7, the
FSI order study) loads `libsolids4FoamModels` through the controlDict's
`libs` for the `bdfDdtScheme`. That scheme is not in any commit of this
branch: it exists in the working tree as the untracked directory
`src/solids4FoamModels/numerics/bdfDdtScheme/` (with one added line in
`Make/files.openfoam`), from the other session's work on PR #674
(`feature/bdfDdtMovingMesh`), and the three files differ from that
branch's committed version (they are newer). Their SHA-256 (the files
the results were made with):

| file | sha256 |
|---|---|
| `bdfDdtScheme.C` | `789e4ec73a35cb855615f73fb9fbc7f2ad6ef9590ea3eb3b110165e5b11d705b` |
| `bdfDdtScheme.H` | `ccdd203d7fd3e558b4531cfeb09f54c7a7c8356e47c48eefae5340d5477e4066` |
| `bdfDdtSchemes.C` | `d2bd25907efe491b982314593890a4c454f675ebb51a33bf5317e3f962a81ee3` |

The library that supplied it on this Mac is the scratch build
`platform/lib/libsolids4FoamModels.dylib` of 5 November (the working
tree's solids4FoamModels, with these files), loaded by pimpleFoam through
`libs`; pimpleFoam itself is OpenFOAM's. For the paper's numbers to come
from a committed state the scheme has to be committed (on the PR #674
branch, or merged here): a decision for the author, since it is the
other session's work.

**The other session's changes are still uncommitted** (the 0D coupler
`circulationCoupling`, `circulation0D` components, `circulationMeshMotion`,
the `leftVentricle` 0D tutorial, `pressureInletOutletDdtCorrVelocity`,
the BDF scheme above, `CHANGELOG.md`, `src/circulationModels/Allwmake`
and README). None of them is in the commits of this work and none is
needed by the 3D-1D cases except the BDF scheme.

**The paper plan and the literature map** were not found on this machine
(the paths in the instructions were placeholders); the work does not
need them yet.

**xenosim.** One node, 192 cores (191 idle on 16 November, one user's
68-core jobs queued on `long`), Slurm partitions `main` (2 days), `long`
(7 days), `verylong` (28 days), Ubuntu 22.04, OpenFOAM v2412 from the
distribution package (`/usr/lib/openfoam/openfoam2412`), PETSc and
Eigen from packages, no solids4foam. The clean HEAD (`git archive
726e5ec61`) plus the three BDF files above is built there in
`~/OpenFOAM/davidkelly-v2412/solids4foam-paper` (`env-paper.sh` sets
the environment, `build.slurm` rebuilt it in 5 minutes on 16 cores);
the job scripts of this round live in that directory. Jobs use at most
80 cores in all.

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

### 7. The benchmark package (done)

**Specifications.** `src/circulationModels/docs/BENCHMARKS3D1D.md` now
carries, for each of the five cases, a table another code can follow
without this one: geometry, material data, inlet and outlet conditions,
initial state, the interface condition, the exact (or numerical)
reference with its formula, the quantities of interest, the error norms
and the refinement path: 1 (rigid segment), 1b (Womersley segment), 1c
(resonator), 1d (FSI tube, numerical reference), 1e (network, numerical
references). Each states the two conventions found in item 1 and in the
"what changed" list: cell averages against cell averages of the closed
form, and the reference taking the area the 3D mesh has (32-gon,
equal-area wedge).

**Exact references as CSV (regenerated in seconds, ignored by git).**

| case | command | files |
|---|---|---|
| rigid segment | `python3 src/circulationModels/docs/benchmarks3D1D.py csv` | `docs/reference/rigidSegment_pipe.csv` (Q, P1, P2 every 0.1 ms), `rigidSegment_{proximal,distal}_t0.35.csv` (P, Q, A every 0.5 mm) |
| Womersley segment | `womersleySegment/analyse.py --csv` | `reference/alpha{3,6,12}_ports.csv` (0 to 6 s every 0.1 ms), `_periodic.csv` (complex amplitudes, Gamma, T), `_profiles.csv` (u(r, t), 101 radii, 8 instants) |
| resonator | `dobroserdova2016/analyse.py --csv` | `reference/{sine,step}_c{350,700,1050}_{hard,soft}_pipe1D.csv`, `..._hard_pipe3D.csv` (every 0.1 ms; the 3D variant with the 32-gon area) |

**Numerical references at physiological amplitude (committed, labelled
numerical in their header).** The nonlinear 1D equations have no closed
form at this amplitude, so the finest run is the reference and the next
coarser run bounds its discretisation error.

| case | reference run | files | coarser run against it |
|---|---|---|---|
| rigid segment (`rigidSegment3D/highOrder.py --reference`) | inflow x 1 (Boileau Eq. 14), implicit 4, 6400 cells per vessel, dt = 0.125 ms, kExact 3 + RK4, BDF 4 pipe | `rigidSegment3D/reference/rigidSegment_physiological_pipe.csv`, `_{proximal,distal}_t0.35.csv` | 4e-7 of the peak (the fourth-order sequence of `highOrder.md`) |
| Womersley segment (`womersleySegment/physiological.py`) | alpha = 6, incident wave 5 kPa (0.12 rho c0^2) over the first second from rest, BDF 4 and coupling 4 at dt = 0.25 ms, kExact 3 + RK4 on 800 cells, 100 x 200 wedge (1 h on 2 cores for the two runs) | `womersleySegment/reference/alpha6_physiological_ports.csv` (P_in, Q, P_out every 0.25 ms), `_profiles.csv` (u(r, t) at 9 instants) | 50 x 100: 5.1e-6 / 5.7e-6 / 5.8e-6 (RMS / RMS, P_in / Q / P_out), 7.6e-4 on the profile (`physiological.md`) |

At 5 kPa the Womersley run differs from the linear exact solution scaled
to that amplitude by 6 % (P_in), 11 % (Q) and 5 % (P_out) RMS: that is
the 1D nonlinearity (at the peak U = P/(rho c0) = 0.77 m/s is 0.12 c0,
and the sqrtElastic wall stiffens by another 3 % in c, so the forward
characteristic runs about 15 % faster than c0 there; inferred from the
equations, not measured), a measure of what the linearised tests leave
out, and the reason the physiological variant needs a numerical
reference.
The two mesh levels agree to 6e-6 at the ports, the same spacing as the
linearised study of item 3, so the reference is converged to that.

## Phase 2

### 8. The compliant FSI tube (B4): (a) and (b) verified, (c) awaiting a decision

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
case was rerun with it first to see whether the library change itself
moves the results (it does: next paragraph).

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

The shell's static edge stiffness (1.04e7 Pa/m, i.e. 2 D beta^3 = (E h/R^2)/(2 beta), the hoop stiffness of only 0.87 mm of tube; an earlier version of this
sentence called it "a sixth of the hoop stiffness", which was wrong) is
weak, so the elastic support behaves nearly as the free end; the full tube lies between clamped and free/spring at every
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

**Decision (16 November): not built.** The paper's limitations state it
in one sentence: the 3D solver is not tested along the direction of wave
travel, since in the rigid and viscous segments the pressure wave
crosses the 3D part instantaneously (incompressible) and in the FSI
segment the wall carries it; a 3D segment that itself propagates the
wave at the vessels' speed (a compressible fluid of matched impedance,
or a long compliant segment resolved in 3D) is left for later work.

### An impedance (Robin) interface condition

Not touched, as asked. For the record, the ingredient is already in the
code: within one step each vessel end answers the port flow as a resistor
Zc (the coupler's documentation and BENCHMARKS3D1D.md), so a patch
condition `p = P* + Zc (Q - Q*)`, with P* and Q* the network's prediction,
would make the outer iteration converge in one pass for a rigid segment
(the model problem's contraction factor dt/(a0 tau) becomes zero) and
would be the natural place for an "impedance" condition. It is method
development and changes the coupler, so it waits for a decision.

## Second round (16 November 2026): the decisions on the first round

Numbered as in the instructions of 16 November. Item 0 is the
"Housekeeping (second round)" section above.

### R3. Resonator with a smoothly ramped inflow (verified)

`dobroserdova2016/makeCases.py` gains the inflow `smooth`, a Gaussian
pulse `1e-7 exp(-(t - 0.105)^2/(2 sigma^2))` m3/s with sigma = 15 ms (7
sigma from the start: 2e-11 of the peak at t = 0), beside the published
sine and step, which are kept unchanged (the tutorial's Allrun now runs
27 1D and 3 pimpleFoam cases; `results.md` gains the smooth rows).
`orderStudy.py --case dobroserdova2016 --exact --inflow smooth` (20 runs,
13 minutes with 3 at a time; `orders_smooth_linear.md`,
`plots/orders_smooth_linear.png`): the hard-coupling 3D run with
linearised fourth-order vessels against the delay-equation solution,
which is now smooth. RMS error / RMS of the signal, port flow (pressure
the same to two digits):

| dt [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 2 | 1.59e-1 | 1.73e-2 | 1.59e-3 | 1.55e-4 |
| 1 | 1.55e-1 | 4.16e-3 | 1.94e-4 | 9.02e-6 |
| 0.5 | 8.14e-2 | 9.26e-4 | 2.20e-5 | 5.25e-7 |
| 0.25 | 4.12e-2 | 2.31e-4 | 2.74e-6 | 3.31e-8 |
| 0.125 | 2.07e-2 | 5.77e-5 | 3.44e-7 | 2.53e-9 |
| orders | 0.04, 0.93, 0.98, 0.99 | 2.06, 2.17, 2.01, 2.00 | 3.03, 3.14, 3.00, 3.00 | 4.11, 4.10, 3.99, 3.71 |

Every family at its nominal order against the exact solution to the
finest step, with no floor above 2.5e-9 (BDF 4's last pair, 3.71, is
where the 1D vessel's error and round-off begin); the self-convergence
orders of the same runs are 0.98, 2.01, 3.00, 4.00 on the finest
triple. This replaces item 2's resonator table as the paper's
time-order evidence: the sine's floors (2.2e-6 in Q, 5.3e-4 in P) were
the problem's jumps, as item 2 said, and the smooth pulse shows the
schemes without them. The published sine and step stay in the paper for
the comparison with Dobroserdova et al.

### R4. Exact solutions independent of the mesh: the cost (awaiting the go-ahead)

**What changes.** The 3D meshes of the rigid segment and the resonator
(the same 32-gon pipe) get their vertex radius scaled by
`sqrt(pi/(16 sin(pi/16)))` = 1.003220, and the wedge of the Womersley
segment by `sqrt(theta/sin theta)` = 1.000635 (theta = 5 degrees), so
that every patch has the area `pi R^2` of the vessels and the closed
forms use `pi R^2` only. `A3D` and `A_WEDGE`/`R_EQ` disappear from
`dobroserdova2016/analyse.py`, `womersleySegment/analyse.py`,
`benchmarks3D1D.py` and `orderStudy.py`; the benchmark tables 1, 1b and
1c lose their "the area the mesh has" rows; the uniform mesh family of
item R2 is already built this way (`spaceStudy.py --family uniform`).
Since the pipe's inertance changes by 0.3 to 0.6 %, every number of
every 3D table changes in its third digit, so all of them have to be
regenerated, not only the exact-solution comparisons.

**Cost of regenerating the headline tables** (from the logged wall
times of the existing runs, serial unless stated):

| table | runs | serial time | notes |
|---|---|---|---|
| rigid segment: `Allrun-sweep` (results.md: accuracy, stability, fine) | about 95 | 1.5 h | 6 at a time: 20 min |
| rigid segment: `highOrder.py` self-convergence and `--linear` exact (highOrder.md, highOrder_linear.md), physiological reference | 41 | 2.2 h | 6 at a time: 40 min |
| resonator: tutorial 3D runs, `orderStudy.py` (self, `--exact`, `--exact --inflow smooth`) | 63 | 1.5 h | 3 at a time: 35 min |
| Womersley: tutorial alpha 3, 6, 12 (6 s each) | 3 | 2.5 h | 3 at once: 1 h |
| Womersley: `orderStudy.py` (self and `--exact`) | 32 | 6 h | 3 at a time: 2 h |
| Womersley: spatial order (R2's new family, counted there) | - | - | - |
| Womersley: physiological reference (50 x 100 and 100 x 200, 1 s) | 2 | 4 h | the 100 x 200 run is the long pole (3.7 h with the graded mesh; about 1 h with R2's uniform mesh and settings) |
| **total** | about 240 | **18 h** | **on this Mac 6 to 8 at a time: 4 to 5 h; on xenosim all at once (about 40 cores): about 4 h, set by the Womersley physiological run** |

The FSI tube and the network are not affected (their meshes are not
compared with a closed form). The first-round text of items 1, 2, 3, 7
and the "what changed" bullets on the areas would be rewritten: the
areas stop being a finding and become a rule of the specification
("the 3D mesh must have the vessels' area"). Not started.

### R9. The network's time order on the cylinder segment (verified)

`orderStudy.py --case thoracicAorta37 --cylinder` (16 runs, 50 minutes
with 2 at a time; `thoracicAorta37/orders_cylinder.md`,
`plots/orders_cylinder.png`): the first-round order study (item 4's
network, smooth start, CFL 0.1, p residual 1e-10) with the cone replaced
by the cylinder of its geometric-mean radius (10.107 mm; `floorStudy.py`'s
variant), all four families. Self-convergence orders, port flow / port
pressure:

| dt triple [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 2 / 1 / 0.5 | 0.95 / 0.87 | 1.99 / 2.01 | 2.99 / 2.99 | 3.83 / 3.72 |
| 1 / 0.5 / 0.25 | 0.98 / 0.94 | 2.00 / 2.01 | 3.02 / 3.02 | 3.82 / 3.91 |

(The mid-plane flow and pressure of the 3D segment give the same orders.)
For the paper this replaces the cone's table as the network's time-order
result: the coupling takes the 3D scheme's order in the 37-artery network
when the 3D mesh is orthogonal; the cone's BDF 3/4 floor is pimpleFoam's
on a non-orthogonal mesh (item 4, `floor.md`), and is stated as such.

### R1. The 7 % between the two solids4foam builds: diagnosed (in progress)

**Set-up.** The solids4foam `3dTube` tutorial (no 1D side: the
Robin-Neumann FSI pulse tube, 20 ms at dt = 25 us, outerCorrTolerance
1e-6), run with three builds on this Mac: the installed February 2026
`libsolids4FoamModels` + `solids4Foam` (what the report's FSI numbers
used), the checkout's working tree (`platform2`, the item 8 runs), and a
clean `git archive 726e5ec61` without the uncommitted files
(`platform3`). The wall's radial displacement at the tutorial's point
(z = 2.5 cm), 800 steps:

| build | peak [um] | at [ms] | FSI it/step | RMS / peak vs checkout | max / peak |
|---|---|---|---|---|---|
| installed (Feb 2026) | 149.86 | 7.22 | 3.03 | 2.0e-2 | 6.7e-2 |
| checkout, working tree | 146.51 | 7.12 | 3.36 | 0 | 0 |
| checkout, clean HEAD | 146.51 | 7.12 | 3.36 | 0 (bit-identical) | 0 |
| installed, tolerance 1e-8 | 149.87 | 7.22 | 4.12 | 2.0e-2 | 6.7e-2 |
| checkout, tolerance 1e-8 | 146.51 | 7.12 | 6.13 | 1.0e-7 | 6.8e-7 |

So: (i) the difference is in solids4foam itself, not in the 1D
coupling (the tube alone shows it); (ii) the uncommitted files (BDF
scheme, 0D coupler, mesh motion) are not involved (clean HEAD is
bit-identical to the working-tree build); (iii) the FSI tolerance is not
the cause (tightening it by 100 moves either build by less than 3e-5 of
the peak; the two builds stay 2 % apart). The same case on xenosim with
the clean HEAD (gcc, Ubuntu package OpenFOAM) agrees with the Mac's
checkout build to 3e-7 of the peak.

**Which commit.** The candidate found first, `28751407c` of 6 October
("FSI: fix two defects in the wall force and the interface point
transfer": a tmp-aliasing error in the wall viscous force at gcc -O3, and
degenerate barycentric weights in the interface point transfer), changes
the 3dTube by 5e-4 of the peak (built and run on xenosim at that commit
and at its parent), so it is not the 2 %. Of the 288 commits touching the
solid, fluid and FSI code between the two builds the suspects are the
solid-model rewrite of September 2026 (the `mechanicalConstitutiveLaw`
framework made the solid models' only path, `25c7a7758`, with the
hydrostatic smoothing `902c51d2c`), the FSI predictor made the default
(`7709d76a4`), and the elasticWallPressure changes of 22-24 September. An
automated bisect is running on xenosim (the clean 3dTube case of HEAD
run for 8 ms on 4 ranks at each commit, 12 minutes a step; the Mac's
shared clone `scratchpad/bisclone`, script `scratchpad/bis/step.sh`,
results in `scratchpad/bis/results.txt`); its first step checks that the
development commit of 4 February (`4b254fa52`) reproduces the installed
build.

**Which build is right.** Not decided until the bisect names the
commit; what can be said now: the checkout build is the one that
corresponds to a committed state of solids4foam, it passes the
tutorial's own CI regression tests, and the February binary cannot be
tied to a commit. If the bisect lands on a fix (as `28751407c` was, for
its 0.1 %), the February numbers were wrong and the paper uses the
checkout; if it lands on a modelling change (e.g. the constitutive-law
rewrite), that is a change in solids4foam to raise with its developers,
with the 3dTube case as the evidence.

### R7. The Robin (impedance) interface condition (in progress)

**Implemented** as an option of the coupler, `coupling { interface
robin; }` (default `dirichlet`, the coupling as before; the Dirichlet
path reproduces the stored rigid-segment baseline to the last digit,
1e-16). The patch condition becomes

    p = P*/rho + zu (u_n - u_n*),   zu = Zc flowScale A3D / rho,

per face: P* the vessel end's pressure for the flow Q* the network was
integrated with (the iterate; for explicit coupling the extrapolated
flow), u_n* the face velocity of that flow, Zc = rho c/A the vessel
end's characteristic impedance at its current state (`network1D::
portImpedance`). It is a mixed condition, p = f R + (1 - f) p_P with f =
1/(1 + k), k = zu rAU_f deltaCoeffs, R = P*/rho + zu (phiHbyA - phi*)/
|Sf|, so the pressure equation sees the 1D side's response to the flux it
is computing; rAU_f is recovered exactly from the solver's
`constrainPressure` call (the condition implements `updateableSnGrad`)
and `phiHbyA`. At the fixed point u_n = u_n* on every face and p = P*
uniform: the converged solution is the Dirichlet one (verified: the
rigid segment at dt = 0.5 ms agrees with the stored Dirichlet run to
6e-9).

**The scalar model, extended.** Within a step a vessel end answers the
port flow as P = P0 - Zc Q (inlet port; + at the outlet), P0 the
incoming wave. With the Dirichlet interface the pipe step sees P(Q*) of
the previous iterate: the fixed point contracts by dt/(a0 tau) and the
explicit scheme is stable below the dt/tau limits of the first round.
With the Robin interface the step sees P(Q*) -/+ Zc (Q - Q*) = P0 -/+
Zc Q: the iterate drops out, so the implicit iteration converges in one
exchange and the explicit scheme has no dt/tau limit at all for
linearised vessels (what it extrapolates is only P0, which within a
step does not depend on Q). For nonlinear vessels the residual loop gain
is (1 - Zc_imposed/Zc_tangent) dt/tau, a few per cent of the Dirichlet
one.

**The 1D surrogate (verified; `Test-circulation1DPorts` gains
`interface robin`).** BDF2 pipe of tau = 6.48 ms, the tutorial's ports:

- implicit order 2 at dt = 0.5 ms (dt/tau = 0.077): with ONE exchange
  per step the Robin run agrees with the converged Dirichlet run to
  3.2e-9 of the peak (the Dirichlet run with two exchanges: 1.6e-3);
  the Dirichlet iteration contracts by 0.051 per pass as the model says;
- explicit, orders 1 to 4, dt/tau = 1, 2, 4, 8, 16 (the Dirichlet limits
  are 4, 1.33, 0.57, 0.27): every Robin run is stable (the peak flow
  stays within a factor 3 of the resolved pulse's at dt/tau = 16, where
  the step is 100 ms and the pulse is unresolved); every Dirichlet run
  beyond its limit diverges.

**In 3D (rigid segment, implicit order 2, dt = 0.5 ms; in progress).**
The interface pressure P* converges in 2 to 3 exchanges (to 1e-10
relative) where the Dirichlet needs 9; but the PIMPLE loop then needs 17
outer iterations instead of 9 to bring the pressure residual to 1e-10,
because the mixed condition couples the boundary pressure to the
momentum predictor's flux (the p residual then decays by 0.4 per outer
iteration, a mode the Dirichlet condition does not have). So at this
small dt/tau the Robin interface converges the coupling faster but the
step is not cheaper.

**In 3D, measured (verified; `scratchpad/item7`).** PIMPLE outer
iterations per step to a p residual of 1e-10 (the stored cases' own
criterion), the converged solutions agreeing to 1e-8 or better:

| case | dt | dt/tau | Dirichlet | Robin (mixed) |
|---|---|---|---|---|
| rigid segment, implicit 2 | 0.5 ms | 0.077 | 11.1 | 15.3 (coupling converged after 2-3) |
| rigid segment, implicit 2 | 2 ms | 0.31 | 19.6 | 28.4 |
| rigid segment, implicit 2 | 4 ms | 0.62 | 22.6 | 29.0 |
| rigid segment, implicit 2 | 8 ms | 1.23 | 30.2 | 30.4 |
| resonator (sine, 3D), implicit 2 | 1 ms | - | 11.0 | 18.5 |
| 37-artery network, implicit 2 | 1 ms | - | 5.0 | 15.0 |

At dt = 2 to 8 ms the Dirichlet count is no longer the coupling's
(dt/(a0 tau) = 0.2 to 0.8 would need 14 to 100 passes): it is the 3D
solver's own, and the Robin count is the same or higher. So the mixed
condition does not save outer iterations in pimpleFoam; what it does is
remove the explicit limit:

| rigid segment, explicit order 2 | dt | dt/tau | Dirichlet | Robin |
|---|---|---|---|---|
| | 16 ms | 2.5 | diverges (2300 x the peak) | stable, max error 0.36 of the peak (16 ms steps on a 50 ms pulse) |
| | 32 ms | 4.9 | - | stable |

(The Dirichlet limit for order 2 with BDF2 is dt/tau = 1.33.)

**The relaxed Dirichlet form (`interface aitken`).** Since the Robin
fixed point is the Dirichlet one under-relaxed by omega = 1/(1 + g), g
the loop gain (dt/(a0 tau) in the model), the same response can be
applied without the mixed condition: Aitken's secant on the port
pressures over the exchanges of a step measures omega (first exchange of
a step at the previous step's omega; omega kept in (0, 1]); the boundary
condition stays Dirichlet, so the 3D solver's iteration is untouched.
Measured (verified, same cases and criterion; the converged solutions
agree with the Dirichlet ones to 7e-9 or better):

| case | dt/tau | Dirichlet | Robin (mixed) | Aitken-relaxed Dirichlet |
|---|---|---|---|---|
| rigid segment, 0.5 ms | 0.077 | 11.1 | 15.3 | 9.2 |
| rigid segment, 4 ms | 0.62 | 22.6 | 29.0 | 19.5 |
| rigid segment, 8 ms | 1.23 | 30.2 | 30.4 | 23.6 |
| resonator, 1 ms | - | 11.0 | 18.5 | 9.5 |
| network, 1 ms | - | 5.0 | 15.0 | 5.0 |

The relaxed form saves 15 to 20 % of the outer iterations where the
coupling's contraction is a visible part of the count, and nothing on
the network, where five passes are the 3D solver's own. The conclusion
for the paper: in a segregated pressure-based 3D solver the outer
iterations are mostly the solver's, not the coupling's, once dt/tau is
below about 0.1 (as it is in the physiological cases: 0.02 on the
network); the Robin interface's gain is the explicit scheme's stability,
not the implicit scheme's cost.

**Order study with the Robin interface (verified; `rigidSegment3D/
highOrder.py run 6 --linear --robin`, 20 runs, 1.5 h with 6 at a time;
`highOrder_linear_robin.md`, `plots/highOrder3D1D_linear_robin.png`).**
Against the closed form, implicit order 4, 2 and 1 give the same errors
as with the Dirichlet interface to three digits at every level (the
same fixed point, orders 3.93 to 4.00 on the flux and pressures, 2 and 1
on the waves for orders 2 and 1); explicit order 4 keeps fourth order
at every level with the Robin interface, its flux equal to the
implicit's to three digits and its pipe-end pressure error 1.35 times
the Dirichlet-explicit one (5.62e-3 to 8.67e-8 over the five levels,
orders 3.99 to 4.00). So the Robin interface changes neither the
converged solution nor the order of either scheme.

**Status of item 7: done.** Files: `circulation1DCoupler.{H,C}`
(`interface dirichlet | robin | aitken`, `portImpedance`, the iterate
fluxes), `circulation1DPressureFvPatchScalarField.{H,C}` (now a mixed
condition with `updateableSnGrad`; Dirichlet mode bit-identical to the
old fixed value), `coupledFlow::characteristicImpedance`,
`network1D::portImpedance`, `Test-circulation1DPorts` (`interface
robin`), `highOrder.py --robin`, BENCHMARKS3D1D.md (the model extended).

### R2. Womersley spatial order: the stall cured (in progress)

**Diagnosis (verified, short runs of 0.25 to 0.3 s at alpha = 6).** On
the tutorial's wall-graded wedge the outer loop's slow convergence is the
segregated algorithm's, and it worsens with refinement (18.6 passes per
step at 50 x 100, 45 at 100 x 200, to a p residual of 1e-10). Tried:

| mesh | predictor | p solver | passes per step (wave present) | notes |
|---|---|---|---|---|
| 50 x 100 graded (tutorial) | off | PCG | 18.6 (1 s window) | the first round |
| 50 x 100 uniform | off | PCG | 8.0 | converges to 1e-10 every step |
| 50 x 100 uniform | on | PCG | 6.4 (nCorrectors 1), 5.3 (3) | 6 steps near the wave's arrival do not reach 1e-10: with the predictor on, the second corrector's p residual sits at 5e-9 whatever the outer iteration, because the port condition is evaluated on the predicted flux in the first corrector and on the corrected one in the second |
| 50 x 100 uniform | off | GAMG | 8.0 | no faster than PCG (57 s against 65 s) |
| 100 x 200 uniform | off | PCG | 14.2 | 1 s in about 50 minutes |
| 100 x 200 uniform | off | GAMG | 14 | slower than PCG |
| 100 x 200 graded | on | PCG | 9.8 | same non-convergence near the arrival |
| 200 x 400 uniform | off | PCG | - | 14 h per second of simulation: not used |

So the family for the paper is the **uniform mesh with the predictor off
and PCG** (the only combination that converges to 1e-10 at every step),
at 25 x 50, 50 x 100 and 100 x 200; the wedge is scaled by sqrt(theta/sin
theta) so that its cross-section is pi R^2 and the exact solution is that
of the pipe of radius R, independent of the mesh (`spaceStudy.py
--family uniform`, runs `space_alpha6_u_x*`). At alpha = 6 the Stokes
layer (2.4 mm) has 6, 12 and 24 uniform cells on the three levels. The
three runs are on xenosim (`space_u.slurm`; 25 x 50 in 2 minutes,
50 x 100 in 20, 100 x 200 in about 2 hours).

**Result (verified; `space_alpha6_u.md`, `plots/space_alpha6_u.png`),**
alpha = 6, 0 to 1 s from rest, BDF 4 and coupling order 4 at dt = 0.5
ms, linearised vessels at fourth order, against the exact transient of
the pipe of radius R:

| nr x nz | P ports | order | Q | order | profile (max over r) | order | profile, r < 0.97 R | order | passes/step | time [s] |
|---|---|---|---|---|---|---|---|---|---|---|
| 25 x 50 | 1.35e-04 | | 1.35e-04 | | 4.49e-02 | | 7.41e-03 | | 5.8 | 69 |
| 50 x 100 | 3.36e-05 | 2.01 | 3.36e-05 | 2.01 | 2.39e-02 | 0.91 | 2.33e-03 | 1.67 | 8.4 | 322 |
| 100 x 200 | 1.05e-05 | 1.68 | 1.05e-05 | 1.68 | 1.33e-02 | 0.84 | 1.51e-03 | 0.63 | 14.8 | 3508 |

Every step converges to the 1e-10 p residual except 34 to 62 steps per
run at the start, at t = 0.52 and at t = 0.99 s, where the normalised
residual cannot be met because the pressure field passes through zero
(the normalisation vanishes); they are harmless. Two things are not
clean yet and are being checked:

- the profile's max-norm error sits at the last sample point, 5 um from
  the wall (r = 0.9995 R, the `sets` function object's cellPoint
  interpolation between the wall value and the first cell), i.e. it is
  the sampling's first-order interpolation, not the solution's; away
  from the wall (r < 0.97 R) the error is 3 to 30 times smaller but its
  order is not 2 either on the finest pair;
- the port error's order drops from 2.01 to 1.68 on the finest pair,
  and the error grows with time within each run (level 2: 1.9e-6,
  4.7e-6, 1.3e-5, 1.6e-5 of the signal over the four quarters of the
  second; level 1: 8.3e-6, 1.8e-5, 5.1e-5, 3.8e-5), with a ratio between
  levels that falls to 2.4 in the last quarter: a non-spatial
  contribution appearing after about 0.5 s, i.e. after the pipe's
  reflection has travelled back to the vessels' ends (2 x 0.162 s). Two
  checks run on xenosim (`space_u2_checks.slurm`): the 100 x 200 level
  with the 1D side refined (800 cells, CFL 0.25) and with the time step
  halved; a fourth level, 200 x 400, runs on the long partition (about
  a day). Whichever moves the late-time error names the floor.

### R6. The FSI tube with pulses of increasing length (verified)

`formaggia2007/pulseSweep.py run` on xenosim (8 runs, 13 h on 20 cores;
`pulse.md`, `plots/pulse.png`): smooth sin^2 pulses of 1333 Pa and
duration T = 3, 6, 12, 24 ms (R/(c0 T) = 0.32, 0.16, 0.079, 0.040) on the
full FSI tube (15, 20, 20, 25 cm, p = 0 at the far end, long enough that
no fluid wave returns within the window T + 30 ms) and on the coupled
5 cm tube + 1D vessel (clamped wall end, c0 = 5.25 m/s, the static
plane-strain match), dt = 25 us, checkout build of solids4foam.

**The plain RMS difference does not fall cleanly with pulse length**
(max over the 3D part 3.1, 4.0, 4.0, 3.1 % of the pulse; at z = 12.5 cm
11.4, 14.0, 10.4, 6.9 %). Split into its parts, it does where it should
and two other parts remain, both identified:

| T [ms] | R/(c0 T) | full tube c, z 1.25-5 / 5-12.5 cm [m/s] | peak ratio coupled/full at 12.5 cm | shift that aligns at 12.5 cm [ms] | RMS after the shift, 12.5 cm | max RMS, 3D part (z <= 5) |
|---|---|---|---|---|---|---|
| 3 | 0.317 | 3.77 (dispersed) | 1.377 | 1.65 | 3.8e-2 | 3.1e-2 |
| 6 | 0.159 | 5.09 / 5.17 | 1.323 | 1.55 | 5.8e-2 | 4.0e-2 |
| 12 | 0.079 | 5.02 / 5.11 | 1.129 | 1.15 | 4.4e-2 | 4.0e-2 |
| 24 | 0.040 | 5.25 / 5.06 | 1.021 | 1.17 | 1.6e-2 | 3.1e-2 |

1. **Dispersion (falls, as expected).** The thick wall's inertia
   disperses short pulses, which the 1D model cannot: the coupled run's
   peak at z = 12.5 cm is 38 % above the full tube's for 3 ms and 2 %
   for 24 ms, and the shape difference left after aligning the arrival
   falls to 1.6 %. This is the long-wave mismatch, and it goes as the
   pulse lengthens.
2. **Wave speed of the static match (does not fall).** In the 1D part
   the coupled pulse arrives 1.2 ms early at z = 12.5 cm for every
   pulse of 12 ms or more. Away from the clamped inlet the long-pulse
   wave speed of the full tube is 5.06-5.11 m/s, the static match of a
   thick cylinder with free ends (generalised plane stress, 5.07 m/s);
   near the clamped inlet it is 5.25 m/s, the plane-strain match used
   for the 1D. So for long pulses the 1D vessel is 3.7 % too fast,
   which is a choice of the wall match (plane strain), not a coupling
   error and not dispersion; matching plane stress would remove it.
   (Inferred from the measured speeds and the two static formulas; not
   rerun.)
3. **The cut end of the wall (does not fall).** In the 3D part (z <= 5
   cm) the difference is 2-4 % of the pulse at every pulse length, with
   no arrival shift: it is the clamped condition on the wall's cut end
   (items 8a and R5), which does not depend on the wavelength.

So the expectation holds for the part it concerns: the model mismatch of
the 3 ms pulse (R/(c0 T) = 0.3) is mostly dispersion and falls by a
factor of about 20 in peak and 3 in shape by R/(c0 T) = 0.04. What
remains at long pulses is the wall match (fixable by matching plane
stress) and the structural end condition. Recommendation for the paper:
quote the long-pulse case as the one in which the coupled-against-full
comparison measures the coupling, and match the 1D wall to the long
tube's long-wave speed.

### R2 (concluded). The Womersley floor was the reference, not the solver (verified)

The finest pair of the uniform family gave 1.68 on the ports, and the
fourth level (200 x 400, 17 h on xenosim) 0.42. Checks at 100 x 200: the
1D side refined (800 cells, CFL 0.25) changes the ports by 1e-8 and the
time step halved by 1e-9, so neither is the cause. The error was linear
(identical relative error at 50, 5 and 0.5 Pa, so not the convection
term), and doubling the pipe length doubled it (so not an end effect at
the ports): a per-length error in the segment's impedance, independent
of the mesh.

**Cause.** A wedge of angle theta with one cell across and a flat outer
face at the distance R_c = R cos(theta/2) from the axis has, at every
distance x from the axis, faces and volumes of width 2 x tan(theta/2)
instead of the sector's x theta, a uniform factor k = 2 tan(theta/2)/theta
= 1.000635. The axisymmetric finite-volume problem it discretises is
therefore the circular pipe of radius R_c with every flux multiplied by
k: Womersley's impedance of radius R_c divided by k (area k pi R_c^2 =
A_WEDGE). The "equal-area pipe" used since the first round (radius
sqrt(A_WEDGE/pi)) has the same inertance but a viscous resistance k
times smaller, and the comparison floors at about 1e-5 against it. The
profile sets also sampled 5 um beyond the flat face (r = 0.9995 R >
R_c), which made the profile's maximum error sit at the wall.

**With the correct reference** (`analyse.py: R_WEDGE, K_WEDGE`; used by
`spaceStudy.py`, `orderStudy.py --exact` and `results.md`):

| family | mesh | ports | order | profile (r <= R_c) | order |
|---|---|---|---|---|---|
| uniform, alpha 6 | 25 x 50 | 1.37e-04 | | 4.64e-02 | |
| | 50 x 100 | 3.43e-05 | 1.99 | 2.54e-02 | 0.87 |
| | 100 x 200 | 8.59e-06 | 2.00 | 1.48e-02 | 0.78 |
| | 200 x 400 | 2.15e-06 | 2.00 | 8.54e-03 | 0.79 |
| graded 10:1 (tutorial), alpha 6 | 25 x 50 | 3.11e-05 | | 6.74e-03 | |
| | 50 x 100 | 7.80e-06 | 1.99 | 8.56e-04 | 2.98 |
| | 100 x 200 | 1.96e-06 | 2.00 | 1.86e-04 | 2.20 |
| graded, alpha 3 | 25 x 50 / 50 x 100 | 5.82e-05 / 1.43e-05 | 2.02 | 7.69e-03 / 2.19e-03 | 1.81 |
| graded, alpha 12 | 25 x 50 / 50 x 100 | 2.67e-05 / 6.88e-06 | 1.95 | 1.40e-02 / 1.79e-03 | 2.97 |

The ports converge at second order on both families at every alpha, to
2e-6 of the signal. So the first-round statement that the 100 x 200
graded mesh "is limited by the outer loop ... about 5e-6, not the mesh's"
was wrong: that was the reference's floor too (the PIMPLE stall, 45
passes per step, is a cost, not an accuracy limit). The profile
converges at second order or better on the graded mesh; on the uniform
mesh its maximum sits in the first cell at the wall, where the sampling
interpolates linearly between the wall value and the first cell centre
(first order), and should be measured on cell values instead.

Periodic state of the tutorial runs against the correct reference:
amplitudes 6.7e-6 to 8.0e-6 (alpha 6 and 12) and 1.1-1.3e-5 (alpha 3),
phases 0.0002-0.0004 degrees (essentially the first round's numbers,
since the periodic state is dominated by the inertance). Time order
against the exact transient (alpha 6, 50 x 100 graded): Euler 0.99-1.00
to the finest step; backward 1.99 before the floor; BDF 3/4 at the floor
from the coarsest step; the floor, 7.1e-6, is the spatial error of this
mesh (7.8e-6 in the space study over the 1 s window).

### R1 (continued). Locating the build difference

The first automated bisect treated a crash as "old" behaviour and found
`df5551e9f` ("Fix PETSc Jacobian preallocation", 1 August), which only
fixes a segfault with the older PETSc of the Ubuntu packages (it does not
change results with newer PETSc): every commit before it crashed on
xenosim for that reason, not because it behaves like the February build.
The difference between the installed February build and the checkout
(2.0 % RMS of the 3dTube wall displacement) splits into:

| range | part of the 2.0 % | commit | what it is |
|---|---|---|---|
| February - 1 August | 1.47 % | `739fa4afa` "Fix for least squares boundary face usage (#276)", 8 June | a bug fix in the solid's least-squares gradient (which 3dTube's solid uses): boundary-face values were used where they should not be (new `useBoundaryFaceValues_DD`/`_sigmaHyd` flags in `solidModel`); the commit itself moved 3dTube's regression reference (max displacement 7.83e-7 to 2.24e-7, mean force 0.0489 to 0.0321) |
| 10 - 24 August | 0.62 % | `53b669539` "Synchronise cached interface geometry and build mappings from reference geometry (#375)", 19 August | a bug fix: the FSI interface applied the fluid pressure with the face normals of the initial configuration for the whole run |
| October | 0.05 % | `28751407c` "FSI: fix two defects in the wall force and the interface point transfer" | a bug fix (tmp aliasing in the wall viscous force at gcc -O3; degenerate barycentric weights) |

Range A (bisected with the PETSc fix patched into each commit,
`scratchpad/bis/driveA.sh`; the February-era commits reproduce the
installed build to 3.5e-4 once patched) ended on four candidates, three of
which do not build on xenosim and were skipped: `55a8fa68f` (#248, block-coupled
incompressible: only new files and build lists, no change to code 3dTube
runs), `1f56ccd34` (#272, tutorial READMEs only), `739fa4afa` (#276) and
`a2e43c33e` (#281, a test application only). Only #276 changes code that
3dTube runs, and it updated 3dTube's own regression values, so it is the
1.47 %.

All three parts are fixes of defects in the February code, so the
checkout build is the correct one. The paper's FSI numbers
should come from the checkout build; the report's (installed-library)
Formaggia numbers carry those defects.

### R4 (go-ahead given 17 November). Meshes whose exact solutions are the plain pipe: done

Set-up committed in `bebb98128`: the rigid pipe's O-grid scaled radially
by 1.003222 (32-gon patch of area pi R^2; A3D = pi R^2 everywhere), the
Womersley wedge with its face at distance R and flowScale
pi/tan(theta/2) = 71.954 (the coupled problem is exactly the pipe of
radius R). Validated before the reruns: B1 baseline 1.09e-3 against the
pi R^2 closed form (the 1D error, as before); 25 x 50 uniform wedge
1.37e-4 against the radius-R solution.

Regeneration on xenosim from `bebb98128` + the BDF files (sha256 as in
the housekeeping section), tree `~/OpenFOAM/davidkelly-v2412/solids4foam-r4`,
35 cores alongside the running FSI refinement and bisect (73 of the 80):

| job | what | cores |
|---|---|---|
| r4_B1 | rigidSegment3D: Allrun-sweep (accuracy, stability), highOrder.py (nonlinear, --linear, --linear --robin), --reference | 8 |
| r4_B3 | dobroserdova2016: Allrun (incl. smooth), orderStudy self / --exact / --exact --inflow smooth, analyse.py --csv | 4 |
| r4_B2tut | womersleySegment: Allrun (alpha 3, 6, 12 periodic state) | 3 |
| r4_B2time | orderStudy womersleySegment self / --exact | 8 |
| r4_B2space | spaceStudy graded alpha 6 (3 levels), uniform alpha 6 (3 levels), alpha 3 and 12 (2 levels) | 6 |
| r4_B2phys | physiological.py (numerical reference, 2 levels) | 2 |
| r4_B4ps | pulseSweep.py --planeStress: the long-pulse coupled runs with the 1D wall matched in plane stress (5.07 m/s), against the existing full tubes (FSI code unchanged between 726e5ec61 and bebb98128) | 4 |

Expected wall time 7-8 h, set by the 100 x 200 graded Womersley run and
the 100 x 200 physiological reference (about 7 h each on xenosim). The
results will replace the tables of items 1, 2, 3, 7, R2, R3 and the
benchmark specifications 1, 1b, 1c.

**Result (18 November; all seven jobs ended cleanly, tables fetched into
the tutorials).** With the meshes posing the plain pipe, every result of
the first round stands: the orders are unchanged and the numbers move in
their third digit at most, which confirms the earlier per-mesh references
(32-gon area; wedge = pipe of R_c with fluxes times k) were right. What
changes is that no reference needs a mesh correction any more.

| table | before (mesh-corrected reference) | now (pi R^2 / radius R) |
|---|---|---|
| B1 tau | 6.5226 ms (32-gon) | 6.4807 ms (= 1D inertance's) |
| B1 closed form, implicit 4, five levels (3D flux) | 1.43e-3 ... 2.36e-8, orders 3.93-4.00 | 1.44e-3 ... 2.36e-8, orders 3.93-4.00 |
| B1 explicit 4, finest | 7.05e-8 | 7.11e-8 |
| B1 time orders (finest pair) BDF3 expl. 4 / BDF4 expl. 4 / BDF4 impl. 4 | 3.17 / 3.89 / 3.78 | 3.16 / 3.86 / 3.70 |
| B1 stability, BDF4 p = 4 (measured vs model 1.1371) | 1.1373 | 1.1376 |
| B3 3D pipe peaks: sine / step / smooth | 1.596 (1.596) / 2.263 (2.265) / 1.698 (1.698) | 1.594 (1.594) / 2.265 (2.267) / 1.700 (1.700): now equal to the lumped pipe's |
| B3 smooth, BDF4 finest | 2.53e-9 | 2.53e-9 |
| B2 periodic ports, alpha 3 / 6 / 12 (amplitude) | 1.1-1.3e-5 / 6.9-8.0e-6 / 6.7-7.8e-6 (vs the wedge's pipe) | the same (vs the pipe of R); phase < 0.0005 deg |
| B2 time floor (50 x 100) | 7.1e-6 | 7.1e-6 |
| B2 physiological, 50 x 100 vs 100 x 200 ports | 5.1-5.8e-6 | 5.1-5.7e-6 |

**A sampling artefact found on the way.** With the wedge's face now at R,
the profile set's last point (r = 0.9995 R) lies inside the wall cell
(before, it was outside the mesh, beyond R cos(theta/2), and was dropped),
so the max-norm profile errors jumped (graded alpha 6 space study 5.9e-3 /
3.1e-3 at order 0.93; tutorial up to 1.4e-2; physiological 2.7e-3). Inside
the wall cell `sets` with cellPoint interpolates between the cell value and
the wall's point values, and the sampled values converge at first order
only: a property of the sampling, not of the solution (the same artefact
as the uniform family's 0.8 order in R2). A first fix, comparing up to the
wall cell's centre, cured the graded family but left the uniform one at
1.90e-2 / 1.12e-2 / 6.53e-3 (orders 0.76, 0.78), the maximum in the wall
cell itself; outside the wall cell the uniform family is at order 2.00
(the graded family's wall cell is ten times thinner and hides it).
`spaceStudy.py`, `analyse.py` and `physiological.py` therefore compare
from the axis to the inner face of the wall cell (r <= R - wall_cell(level),
from the mesh's grading), and the physiological profile reference stops at
the finer run's wall cell. Re-analysed on xenosim, no reruns:

| profile error | all samples | excluding the wall cell |
|---|---|---|
| space study, graded alpha 6, 25 x 50 / 50 x 100 / 100 x 200 | 1.13e-2 / 5.93e-3 / 3.12e-3 (0.93, 0.93) | 3.05e-3 / 7.54e-4 / 1.87e-4 (2.02, 2.01) |
| space study, uniform alpha 6 (first column: to the wall cell's centre) | 1.90e-2 / 1.12e-2 / 6.53e-3 (0.76, 0.78) | 4.84e-3 / 1.21e-3 / 2.99e-4 (2.00, 2.01) |
| tutorial periodic, alpha 3 / 6 / 12 | 3.3e-3 / 6.6e-3 / 1.4e-2 | 3.1e-3 / 3.1e-3 / 2.5e-3 |
| physiological, 50 x 100 vs 100 x 200 | 2.68e-3 | 6.1e-4 |

The space study in full (`womersleySegment/space_alpha*.md`; the ports are
unchanged by the profile norm, P and Q agree to the digits shown):

| family | nr x nz | ports | order | profile | order |
|---|---|---|---|---|---|
| graded, alpha 6 | 25 x 50 / 50 x 100 / 100 x 200 | 3.10e-5 / 7.78e-6 / 1.95e-6 | 1.99, 2.00 | 3.05e-3 / 7.54e-4 / 1.87e-4 | 2.02, 2.01 |
| uniform, alpha 6 | 25 x 50 / 50 x 100 / 100 x 200 | 1.37e-4 / 3.43e-5 / 8.58e-6 | 2.00, 2.00 | 4.84e-3 / 1.21e-3 / 2.99e-4 | 2.00, 2.01 |
| graded, alpha 3 | 25 x 50 / 50 x 100 | 5.80e-5 / 1.43e-5 | 2.02 | 7.69e-3 / 2.19e-3 | 1.81 |
| graded, alpha 12 | 25 x 50 / 50 x 100 | 2.66e-5 / 6.87e-6 | 1.96 | 3.84e-3 / 9.66e-4 | 1.99 |

The uniform 200 x 400 run of R2 (old geometry) was not repeated: the three
levels give order 2.00 already. The physiological reference is unchanged
(5.10e-6 / 5.67e-6 / 5.74e-6 and 6.10e-4).

**FSI wall match in plane stress (r4_B4ps, `formaggia2007/pulse_planeStress.md`).**
The long-pulse coupled runs with the 1D wall at c0 = 5.07 m/s (plane
stress) instead of 5.25 (plane strain), against the same full tubes. In the
1D part (z = 6.25-12.5 cm) at T = 24 ms the time shift that best aligns
coupled and full grows by 0.32 ms (0.35 to 0.67) instead of 0.72 ms (0.45
to 1.17), the peak ratios are 0.99-1.00 instead of 1.01-1.02, and the raw
RMS difference at z = 12.5 cm falls from 6.9e-2 to 4.2e-2 of the pulse; in
the 3D part nothing changes (as it must), and the 0.35 ms shift at the port
(z = 5 cm) is the same in both, so it comes from the 3D part's cut end. The
full tube's speed measured over z = 1.25-5 cm is 5.17 m/s, between the two
models: the thick tube with clamped ends is neither in plane strain nor in
plane stress. Plane stress halves the lag that builds up along the 1D part
but does not remove it; the remaining split is unchanged (dispersion falls
with pulse length, a 0.35 ms port shift from the cut end, a 2-4 % residual
after alignment).

Updated with these numbers: `rigidSegment3D` (README, results, highOrder
tables, references), `dobroserdova2016` (README, results, orders, CSVs),
`womersleySegment` (README, results, space and order tables, physiological
reference), `BENCHMARKS3D1D.md` specifications 1, 1b (wedge geometry and the
profile norm) and 1c, the paper (`sections/exact.tex`, `analysis.tex`) and
the thesis (chapters 7-10, 13).

