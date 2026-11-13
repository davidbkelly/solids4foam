# 3D-1D coupling benchmarks

Reference problems for verifying the 3D-1D coupling (circulation1D phase 4).
Section 1 is implemented and measured (the `inertance` and `coupledFlow`
nodes, `circulation1DCoupling`, the tutorials `rigidSegmentInertance` and
`rigidSegment3D` in `tutorials/circulationModels/coupling3D1D`; results below); Section 2 is not run yet.
`benchmarks3D1D.py` evaluates and checks every formula below:

```bash
python3 benchmarks3D1D.py rigid    # closed form vs ODE, energy, |Gamma|
python3 benchmarks3D1D.py mms      # symbolic check of the 2016 MMS
python3 benchmarks3D1D.py csv      # the rigid segment's reference as CSV
```

The literature has no physical analytical benchmark for 3D-1D coupling:
the coupling is checked against a longer 3D or 1D-only run, or against a
manufactured solution (Section 2). Section 1 proposes one.


## 1. Rigid segment in a compliant line (proposed)

A straight rigid 3D pipe, length `l`, area `A`, joins two 1D vessels of the
same area and wave speed `c`; the far end is absorbing and a pulse arrives
from the left.

**Why it is exact.** The 3D fluid is incompressible and the wall rigid, so
the flow rate `Q(t)` is the same at both ends. With slip walls (or no-slip
and a plug inflow profile in the inviscid limit) the 3D Navier-Stokes
solution is uniform velocity `Q/A` and linear pressure, so

    p_in - p_out = (rho l / A) dQ/dt                      (exact, nonlinear)

The 3D region is an inertance `L = rho l/A` in series. In the linear 1D
lines, `Z_c = rho c/A`, and the load seen from the left is `Z_c + i w L`:

    Gamma(w) = i w tau / (1 + i w tau),   T(w) = 1/(1 + i w tau),
    tau = l/(2c)

In time, the transmitted pressure is a first-order low pass of the
incident one and the reflected pressure is what it removes:

    tau dp_t/dt + p_t = p_inc,     p_r = p_inc - p_t

For a Gaussian `p_inc = P exp(-T^2/(2 sigma^2))` (T from its peak at the
interface):

    p_t = P (sigma/tau) sqrt(pi/2) exp(sigma^2/(2 tau^2) - T/tau)
            erfc((sigma^2/tau - T)/(sqrt(2) sigma))

and `|Gamma|^2 + |T|^2 = 1`: incident energy = reflected + transmitted.

| 5 cm rigid segment, c = 5 m/s | 1 Hz | 2 Hz | 5 Hz | 10 Hz | 20 Hz |
|---|---|---|---|---|---|
| \|Gamma\| | 0.031 | 0.063 | 0.155 | 0.300 | 0.532 |

Everything depends on `w l/c` only: a rigid 3D region is acceptable when
`w l / c << 1`.

### Specification (for another code)

Everything SI. `benchmarks3D1D.py csv` writes the reference below from the
class `RigidSegment` of that script.

| | |
|---|---|
| 1D vessels | two straight vessels, proximal `L1` = 1 m and distal `L2` = 1 m, radius `R` = 10 mm at the diastolic pressure `Pd` = 0, wall thickness `h` = 1.5 mm, Young's modulus `E` = 400 kPa, external pressure 0, inviscid (`mu` = 0), blood density `rho` = 1050 kg/m3 (the single-pulse vessel of Boileau et al. 2015, Table II) |
| 1D model | the linearised 1D equations about the state at rest: `dA/dt + dQ/dx = 0`, `dQ/dt + c0^2 dA/dx = 0`, `P = (rho c0^2/A0) (A - A0)`, with `A0 = pi R^2` and `c0` the wave speed of the sqrtElastic law at `A0`: `c0^2 = 2 E h/(3 rho R)` = (6.1721 m/s)^2; `Zc = rho c0/A0` = 2.0629e7 Pa s/m3. A code with only the nonlinear equations should run the inflow scaled by 1e-3 (then the nonlinear error is 2e-5 of the peak at most) |
| 3D pipe | a straight rigid pipe of length `l` = 80 mm and radius 10 mm between the proximal vessel's downstream end and the distal vessel's upstream end; incompressible inviscid flow (or slip walls), so the exact 3D solution is plug flow; `rho` as above |
| pipe cross-section `A3D` | the exact solution's inertance is `L = rho l/A3D` with `A3D` the area of the 3D mesh's inlet patch, which the 3D solution has exactly. The tutorial's O-grid has a 32-sided wall polygon, `A3D = 16 R^2 sin(pi/16)` = 0.99359 pi R^2 = 3.121445e-4 m2, `tau = L/(2 Zc)` = 6.5226 ms. A mesh with another patch area uses its own `A3D` (`RigidSegment(A3D=...)`); a 1D reduction of the pipe (an inertance element) uses `pi R^2` (`tau` = 6.4807 ms) |
| inflow | prescribed flow at the proximal vessel's upstream end: `Q(t) = 1e-6 exp(-1e4 (t - 0.05)^2)` m3/s (a Gaussian of standard deviation `sigma = 1/sqrt(2e4)` = 7.071 ms peaking at `t` = 0.05 s; the incident pressure peak is `Zc` x 1e-6 = 20.629 Pa). The inlet reflects returning waves with +1 in pressure |
| outlet | non-reflecting end of the distal vessel |
| initial state | everything at rest, `P` = 0, `Q` = 0, `A = A0` |
| interface | continuity of flow and of (static) pressure at both pipe ends |
| window | 0 to 0.35 s: the incident peak reaches the pipe at `t_arrive = 0.05 + L1/c0` = 0.212 s; the reflected pulse has returned to the inlet and re-entered the first centimetres of the proximal vessel by 0.35 s (in the reference) |
| quantities of interest | (i) the flow through the pipe `Q(t)` and the pressures at its inlet `P1(t)` and outlet `P2(t)`, at the coupling times; (ii) the profiles of `P`, `Q` and `A` along both vessels at `t` = 0.35 s |
| reference | `Q = p_t/Zc`, `P1 = 2 p_inc - p_t`, `P2 = p_t` with `p_inc(T) = Zc 1e-6 exp(-T^2/(2 sigma^2))`, `T = t - t_arrive`, and `p_t` the closed form above; along the vessels the incident, reflected (`p_inc - p_t`, backward) and re-reflected (at the inlet, forward) waves, `P = f + g`, `Q = (f - g)/Zc`, `A = A0 (1 + P/(rho c0^2))`. `rigidSegment_pipe.csv` (every 0.1 ms) and `rigidSegment_{proximal,distal}_t0.35.csv` (every 0.5 mm) |
| error norm | relative to the incident peak (`Zc` x 1e-6 for `P`, 1e-6 for `Q`): the maximum over the window for (i), the maximum over both vessels for (ii). A finite-volume code compares its cell averages with the cell averages of the reference (`RigidSegment.profile_avg`, five-point Gauss quadrature per cell): against point values the comparison is second order, `(dx/(c0 sigma))^2/24`, whatever the scheme |
| refinement | `dx` and `dt` halved together from 2.5 mm and 2 ms (the 1D cells; the 3D axial cells may stay at 5 mm, plug flow being exact on any mesh) |
| physiological amplitude | the same with the inflow x 1 on the nonlinear 1D equations: no closed form; reference values from the finest runs of `rigidSegment3D/highOrder.py` (implicit order 4, 0.15625 mm, 0.125 ms, kExact 3 + RK4 vessels, BDF 4 pipe), labelled numerical, see the results below |

**What it tests.** With slip walls, finite volumes represent uniform
velocity and linear pressure exactly, so the 3D discretisation error is
about zero and any mismatch is the coupling's (or the 1D solver's, whose
order is known). It gives:
- the convergence order in dx and in dt; partitioned explicit coupling
  usually drops to first order in time;
- an energy check;
- a number for the "spurious" reflection of a rigid 3D region, which is
  shown only qualitatively in the literature.

**Coupling stability.** The 3D side computes `Q`, the 1D sides return the
pressures. Within one 3D step each vessel end answers as a resistor `Z_c` (its
reflections take a round trip to return), so the pipe sees
`Pin - Pout = 2 P_inc - 2 Z_c Q` and the coupling number is
`z = 2 Z_c dt/L = dt/tau`: the benchmark's own time constant sets the
coupling's stability. A model problem (the 3D side backward Euler or BDF2, the
pressure from the flow extrapolated to `t1` by order `p`) gives explicit
stability for

| explicit order p | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| Euler: `dt/tau <` | 2 | 0.667 | 0.286 | 0.133 |
| BDF2: `dt/tau <` | 4 | 1.333 | 0.571 | 0.267 |
| BDF3: `dt/tau <` | 4.5 | 2.222 | 0.952 | 0.444 |
| BDF4: `dt/tau <` | 1.643 | 2.25 | 1.524 | 0.711 |

(`2/(2^p - 1)` and `4/(2^p - 1)` for Euler and BDF2; BDF4's order 1 goes
unstable through a complex pair of roots, the others at -1, the step-to-step
mode). The plain fixed-point iteration of implicit coupling contracts by
`dt/(a0 tau)` per iteration, `a0` the leading BDF coefficient (1, 3/2, 11/6,
25/12), so it converges for `dt/tau < a0`.
`applications/test/circulation1DPorts` checks this with a rigid-pipe
surrogate for the CFD. For a general 3D domain, `tau` becomes `L/(2 Z_c)`
with `L` the domain's inertance between its coupled ends.

Measured (8 cm pipe, `tau` = 6.48 ms, the 1D ports driven as the coupler
does; milestone M1):

- explicit stability, Euler to BDF4, orders 1-4, at 0.9, 1.1 and 1.25 times
  each limit: every run decays below the limit and grows above it, except
  BDF3 order 1 at 0.9 times (dt = 26 ms, 17 steps after the pulse: a fit of
  1.0005 per step where the model decays by 0.981); where the growth can be
  fitted it is the model's largest root (1.1209 against 1.1210, 1.1508
  against 1.1508). (A first test that compared two windows of a run put the
  limits of orders 3-4 3-4% too high: it missed the slow growth just past
  them.)
- implicit contraction rate: `dt/(a0 tau)` to four digits for BDF1 to BDF4;
- convergence in `dt` against the same coupling at a step 8 times finer:
  explicit order 1 at order 1.00, every other scheme at 2.0 (the BDF2 pipe
  caps it; implicit orders 1-4 coincide); the coupled flow tends to the
  monolithic `inertance` run (2.2e-4 of the peak at `dt/tau` = 0.04).
- the `inertance` node itself against the closed form: order 2 in dx
  (tutorial `rigidSegmentInertance`).

With the CFD (milestone M3, tutorial `rigidSegment3D`: pimpleFoam, slip
walls, `circulation1DCoupling` at both ends, `tau` = 6.52 ms with the 3D
patch area):

- the 3D flux and pipe pressures agree with the closed form to the 1D
  discretisation error (1.1e-3 and 3.5e-3 of the peak); the finest coupled
  run differs from the 1D run with the inertance node by 3.4e-6;
- order in `dt` with BDF2 in the CFD: explicit 1 1.00, explicit 2-4 1.99 to
  2.14, implicit 2 2.00; the errors match the surrogate's to about 1%;
- with `BDF 3` and `BDF 4` (bdfDdtScheme, solids4FoamModels), measured on the
  3D flux: explicit order `p` at `min(p, q)`; implicit coupling of orders 1
  and `q` at `q` (2.98; 3.78 where the error, 5e-7, nears the vessels' time
  error), since within a step the pipe sees only the end-of-step pressure, a
  function of the end-of-step flow. The waves the ports send into the vessels
  carry the shape of the forcing over each step, though: there implicit order
  `p` gives `min(p, q)` (order 1 for I1), and asymptotically that error leaks
  back into the 3D flux too (see below);
- explicit stability: 16 runs at 0.9 and 1.1 times the limits, all as
  predicted, with the model's amplification per step (e.g. 1.1508 against
  1.1508).

Fourth order in space and time (`rigidSegment3D/highOrder.py`): vessels with
the k-exact reconstruction of degree 3 and RK4, the pipe with `BDF 4`, full
amplitude, dx and dt refined together (2.5 mm and 2 ms to 0.3125 mm and
0.25 ms) against a run at half of both:

| coupling | 3D flux | pipe-end pressures | waves in the vessels |
|---|---|---|---|
| implicit 4 | 3.90, 3.98, 4.08 | 3.89, 3.99, 4.09 | 3.95, 4.00, 4.09 |
| explicit 4 | 3.90, 3.97, 3.99 | 3.90, 3.98, 3.97 | 3.96, 4.00, 3.98 |
| implicit 2 | 3.90, 3.98, 4.08 | 3.89, 3.99, 4.08 | 1.60, 1.89, 1.97 |
| implicit 1 | 3.89, 3.91, 2.43 | 3.79, 3.34, 1.87 | 0.98, 1.00, 1.00 |

The coupled 3D-1D problem is fourth order in space and time with fourth-order
coupling (error 4.4e-7 of the peak at the finest level with I4, three times
that with E4). Every quantity converges at `min(p, q)` in the end: a
low-order coupling spoils the waves in the vessels at once, and the 3D side
later, with a much smaller constant (I1's first-order wave error reaches the
3D flux at the finest level; I2's not yet).

Against the closed form itself, with the vessels on the linearised 1D
equations (`linearised true`, see the circulation1D README), of which it
is the exact solution at any amplitude (`highOrder.py --linear`, one level
further, the profiles against the closed form's cell averages):

| coupling | 3D flux | pipe-end pressures | waves in the vessels | error at 0.15625 mm, 0.125 ms |
|---|---|---|---|---|
| implicit 4 | 3.93, 3.97, 3.99, 4.00 | 3.92, 3.98, 4.00, 4.00 | 3.95, 3.99, 4.00, 3.98 | 2.4e-8, 1.8e-8, 2.8e-8 |
| explicit 4 | 3.88, 3.97, 3.99, 3.99 | 3.89, 3.99, 4.00, 4.00 | 3.96, 4.00, 4.00, 4.00 | 7.1e-8, 6.9e-8, 4.7e-8 |
| implicit 2 | as I4 | as I4 | 1.61, 1.89, 1.97, 1.99 | 1.5e-5 (waves) |
| implicit 1 | as I4 | as I4 | 0.98, 1.00, 1.00, 1.00 | 3.9e-3 (waves) |

Fourth order to 2e-8 of the peak with no floor; and the 1D-only reduction
(`rigidSegmentInertance`, the pipe as an `inertance` node, kExact 3 + RK4)
at 3.95-4.01 down to 1.4e-9 (L1) over 200 to 6400 cells per vessel. In the
linear problem the 3D side of I1 and I2 stays at the fourth-order error at
every level: the waves a low-order coupling sends into the vessels never
return to the pipe, and the leak seen in the self-convergence study is the
nonlinear interaction of I1's erroneous reflected wave with the incident
one while they overlap at the pipe. Two rules for anyone comparing with the
closed form: cell averages against cell averages (against point values the
comparison is second order, `(dx/(c0 sigma))^2/24`), and `tau` with the
mesh's own patch area.

**Ladder.**
1. Linear, inviscid, slip walls: the closed form above.
2. No-slip with a Womersley inflow profile: `L` becomes the Womersley
   rigid-tube impedance `i w rho l / (A (1 - F10))`; the 1D lines use their
   viscous characteristic impedance. Closed form in frequency.
3. Nonlinear: the inertance relation stays exact, so the reference is a
   fine 1D run with the pipe as an `inertance` node (`pressure total` for
   total-pressure continuity). At finite amplitude the 1D area differs from the 3D area
   at the interface, which separates static- from total-pressure
   continuity.
4. Compliant 3D (solids4foam FSI): a segment of wave speed `c3D` and
   impedance `Z3D` is a finite transmission line,
   `Z_in = Z3D (Z_c + i Z3D tan(k3D l)) / (Z3D + i Z_c tan(k3D l))`,
   `Gamma = (Z_in - Z_c)/(Z_in + Z_c)`; rigid is the `c3D -> inf` limit and
   matched walls give `Gamma = 0`. The 1D sqrtElastic law has
   `c^2 = 2Eh/(3 rho R)` (the 4/3 is `1/(1 - nu^2)`, nu = 0.5), so the 3D
   wall must match it.
5. Any rigid 3D geometry (e.g. a bifurcation) is, at linear order, an
   inertance matrix `L_ij` from a potential-flow solve.

The formula is textbook acoustics (series inertance in a transmission
line); its use as an exact 3D-1D verification problem is the new part.
Malossi et al. (2013), Passerini et al. (2009) and Formaggia et al. (2001)
still need a full read to confirm nobody has compared with it.


## 2. Dobroserdova, Olshanskii & Simakov (2016), manufactured solution

T.K. Dobroserdova, M.A. Olshanskii, S.S. Simakov, Multiscale coupling of
compliant and rigid walls blood flow models, Int. J. Numer. Meth. Fluids
(2016), doi:10.1002/fld.4241. Transcribed from the preprint (University of
Houston NASC preprint #46, January 2016), Section 4.1; the journal version
was not checked and may differ.

It verifies a 1D-0D-3D chain: 1D vessel, an absorbing 0D "balloon", rigid
3D pipe. It is a manufactured solution with source terms and negative,
time-dependent resistances: it checks an implementation, not the physics.

### Models (paper units: cm, g, s)

1D, `x in [0, b]`, unknowns `S`, `u` (mean velocity), `pbar`:

    dS/dt + d(S u)/dx = phi
    du/dt + d(u^2/2 + pbar/rho)/dx = psi                   (friction off)
    pbar = rho c0^2 f(S/Shat),  f(s) = exp(s - 1) - 1 (s > 1),  ln s (s <= 1)

0D balloon (eq. 7-10), with `Q1D = S u` at `x = b` and `Q3D` the flux into
the 3D domain through `Gamma`:

    I d2V/dt2 + R0 dV/dt + (V - V0)/C = p0D
    dV/dt = Q1D - Q3D
    pbar - p0D = R1D0D Q1D           at x = b
    p0D  - p   = R0D3D Q3D           on Gamma

3D, Navier-Stokes with body force `f`, rigid no-slip wall, traction on the
open ends. Their discretisation enforces `(p0D - p - R0D3D Q3D) n =
-nu du/dn` on `Gamma`; for this solution `du/dn = 0` there.

### Set-up

- 3D: `x in (-1, 1)`, `y^2 + z^2 < 1`; `Gamma` at `x = -1`, outlet at
  `x = 1`.
- 1D: length 5, coupled at its outlet `x = b` to `Gamma`; inlet data from
  the exact solution.
- `Shat = pi` cm^2 (the 3D section), `rho = 1`, `c0 = 350` cm/s,
  `nu = 0.04` cm^2/s, `alpha = 0.11`, `beta = 10` (period 0.1 s).
- 0D: `I = R0 = 0`, `V0 = 0`, `C = 1/c0^2`.
- Final time: not stated.
- Their meshes: 1272, 8403 and 63384 P2-P1 tetrahedra with
  `dt3D = 0.008, 0.004, 0.002`; 1D first order with
  `dx = 5/16, 5/64, 5/256` and `dt1D = 1e-5, 2.5e-6, 6.125e-7`.
- They report second-order 3D velocity (L2) and pressure errors, Tables 1-2.

### Exact solution (corrected)

With `c = cos(2 pi beta t)`, `s = sin(2 pi beta t)`, `E = exp(1 - c)`:

    S    = Shat exp(alpha (c - 1))
    u    = alpha exp(-alpha (c - 1)) (1 - exp(c - 1))
    pbar = c0^2 alpha (c - 1)                       (the ln branch, S <= Shat)

    V    = alpha (c - 1) + (alpha/c0^2)(E - 1)      (V means V - V0)
    dV/dt = 2 pi beta s alpha (E/c0^2 - 1)
    R1D0D = -R0D3D = -E/Shat

    u3D  = ( (2/pi) [Shat alpha (1 - exp(c - 1)) - dV/dt] (1 - y^2 - z^2), 0, 0 )
    p3D  = c0^2 alpha (c - 1) + (2 pi beta/Shat) alpha s E (E/c0^2 - 1)
           + 10 (-1 - x)

`Q1D = Shat alpha (1 - exp(c - 1))`, and the 3D flux is the bracket in
`u3D` (the profile integrates to 1 over the unit disc), so `dV/dt = Q1D -
Q3D` holds.

### Errata in the printed version

`benchmarks3D1D.py mms` checks each relation symbolically. All hold except
eq. 10 with the printed 3D pressure.

1. **3D pressure.** The paper prints the second term as
   `- (2 pi beta/Shat) alpha s alpha E (E/c0^2 - 1)`: wrong sign and an
   extra `alpha`. Eq. 10 then fails by
   `(2 pi beta/Shat) alpha (1 + alpha) s E (E/c0^2 - 1)`, under every sign
   convention for `R0D3D` and `Q3D`. The corrected term is above; it is the
   unique `p` that satisfies eq. 10. The authors' code probably applied the
   correct sources, since their errors converge.
2. **`u3D`** prints `cos(beta pi t)` inside one exponential; eq. 8 needs
   `cos(2 beta pi t)`.
3. **`R1D0D`** prints `cos(2 pi t)`; it needs `cos(2 beta pi t)`.
4. `V` is negative: it is `V - V0` (their Remark 2.1).

### Sources

    phi  = dS/dt = -2 pi beta alpha Shat s exp(alpha (c - 1))
    psi  = du/dt
    f_x  = rho dU/dt (1 - y^2 - z^2) + 4 nu U - 10,
           U = (2/pi) [Shat alpha (1 - exp(c - 1)) - dV/dt]

Outlet traction at `x = 1`: `-p3D(1) n`. The script prints them in full.

The 1D solver here is conservative in `(A, Q)`. There the momentum source
is `psi_Q = S psi + u phi`, since `dQ/dt + d(Q^2/A)/dx + (A/rho) dP/dx =
S psi + u phi`.

### What it covers, and what it doesn't

- The 1D solution is uniform in `x`, so the 1D spatial discretisation is
  not tested; the 1D side is effectively an ODE.
- The 3D velocity is Poiseuille, exact for P2 elements; finite volumes
  will show a spatial error.
- It tests the coupling algebra, the sub-iterations and the time
  integration, with the 0D element in between. It does not test wave
  transmission or reflection.

**Needed to run it here.**
- 1D: the log tube law, or the solution redone with sqrtElastic, and mass
  and momentum sources in `vessel1D`.
- 0D: a compliance between two time-dependent, negative resistances.
  circulation0D resistors are constant and must be positive, so this needs
  a dedicated node or an expression resistor.
- 3D: a body force (fvOptions) and the outlet traction.
