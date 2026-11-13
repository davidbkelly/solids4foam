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

### 10. Reflection sign of the flow inlet (done)

A prescribed flow reflects an incoming wave with Γ = +1 in pressure (−1 in
flow). The code imposes the flow on the outgoing characteristic, so it has
the right behaviour by construction, and the delay-equation script of
`dobroserdova2016/analyse.py` uses `f(t) = Zc Qin(t - td) + g(t - 2 td)`,
the +1 reflection. Only the report's wording in Section 6.2 ("Γ = −1 in
pressure") was wrong: corrected in the report source, and the sign is now
stated in `flowInlet.H` and the `circulation1D` README. Verified by reading
the code and the script; no run needed.
