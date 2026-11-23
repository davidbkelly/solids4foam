# The pulse tube continued by a 1D vessel: `formaggia2007`

Prepared by David Kelly

---

## Tutorial Aims

- Reproduce the 3D-1D test of Formaggia, Moura and Nobile, *On the stability
  of the coupling of 3D and 1D fluid-structure interaction models for blood
  flow simulations*, M2AN 41 (2007) 743-769, section 5, first posed by
  Formaggia, Gerbeau, Nobile and Quarteroni (CMAME 2001): the pressure-pulse
  tube of the `3dTube` tutorial (Greenshields and Weller 2005) continued by a
  1D vessel, against the same FSI tube made longer.
- Couple the 1D network to a solids4foam FSI fluid (the `circulation1DPressure`
  condition on the outlet of a Robin-Neumann FSI case) and match the 1D wall to
  the 3D one.

---

## Case Overview

The tube: radius 5 mm, wall 1 mm, length 5 cm, water (rho = 1000 kg/m3,
nu = 3e-6 m2/s), wall linear elastic (rho = 1200 kg/m3, E = 300 kPa,
nu = 0.3) with both ends clamped; a pressure of 1333.2 Pa at the inlet for
3 ms, then 0. Quarter tube with two symmetry planes, 80 axial cells, 20 ms of
simulation: the `3dTube` tutorial's Robin-Neumann set-up, unchanged, at its
coarser time step of 100 us (`makeCases.py --dt 2.5e-5` for its finer one,
four times the cost). `makeCases.py` writes

| run | what |
| --- | --- |
| `full3D` | their "exact solution": the same tube three times as long (15 cm; theirs was 10 cm, with a reflecting end the pulse does not reach), p = 0 at the far end |
| `coupled` | the 5 cm tube with a 10 cm 1D vessel on its outlet, `coupledFlow` port with `flowScale 4` (quarter tube), absorbing end, implicit coupling of order 2 |
| `only1D` | the pulse (a `pressureInlet` node) in a 15 cm 1D vessel |

**Matching the walls.** The 1D wall is given the static compliance of the 3D
one: a thick cylinder in plane strain (clamped ends), inner radius a = 5 mm,
outer b = 6 mm, has `u(a)/p = a (1 + nu) [(1 - 2 nu) a^2 + b^2] / (E (b^2 -
a^2))`, so `dA/dp = 2 pi a u/p` and `c0 = sqrt(A / (rho dA/dp))` = 5.25 m/s
(the thin-wall Moens-Korteweg value with nu = 0.3 is 5.74 m/s, and 6.32 m/s
with the 1D law's own 4/3 factor, i.e. nu = 1/2). The sqrtElastic law has
`c0^2 = 2 E1D h / (3 rho R)`, so `E1D = 3 rho R c0^2 / (2 h)` = 207 kPa.
`makeCases.py --c0 <m/s>` sets another speed, e.g. the one measured in
`full3D`.

**Sections.** Face zones at z = 1.25, 2.5, 3.75 (in the 3D part of both
runs), 5 (the interface: a zone in `full3D`, the port in `coupled`), 6.25,
7.5, 10 and 12.5 cm (3D in `full3D`, 1D probes in `coupled`) record the mean
pressure and the flow; `solidPointDisplacement` the radial wall displacement
at z = 2.5 cm. A coded function object, `sectionFields`, gives the zones
what they measure: the solver leaves `phi` relative to the mesh and the
sections move axially with the wall, so the flow is the sum of
`phiAbsolute = phi + meshPhi`; on a processor boundary `p` holds the
neighbour cell's value, so the pressure is the area average of `pSection`,
the mean of the two cells of every face. Without them the flows at the
sections were 10-20 % off, and a section on a processor boundary read the
pressure half a cell away (3.5 % of the pulse).

```bash
./Allrun            # full3D (~10 min), coupled (~4 min) and only1D together; results.md, plots/
```

The dictionary of the 1D network is `constant/fluid/circulation1DProperties`
(the coupler registers on the fluid region's mesh). Needs solids4foam with
PETSc; the function objects use `libfieldFunctionObjects` and
`libutilityFunctionObjects`.

---

**Orders used.** As the tutorial: `Euler` in time for the fluid and the
solid (first order), second order in space; Robin-Neumann FSI converged to
1e-6 each step; 1D: MUSCL/MC/HLL in space, SSPRK3 in time; 3D-1D coupling
implicit, order 2 requested, which the first-order CFD scheme limits to
order 1 (the coupler says so at start-up). The published results also used
Euler. The time error is the same in `full3D` and the coupled runs and
cancels in their comparison; for any accuracy statement switch both
`fvSchemes` to `backward` (the commented line) or to `BDF` 3/4.

**Time-order study (`orderStudy.py`, `orders.md`).** The coupled run
repeated at dt = 100, 50, 25 and 12.5 us for the scheme families Euler
(q = 1), backward + `BDF 2` d2dt2 (q = 2), `BDF 3` and `BDF 4` on both
regions, with the coupling order set to q; observed orders from three
successive steps, 0-12 ms. Findings on 2026-11-05: the `BDF` ddt scheme
refuses a moving mesh (`checkStaticMesh`), so BDF 3/4 cannot run on the FSI
fluid as it stands; Euler gives observed orders 0.5-0.8 and backward 1-1.5
(port, wall) with no convergence at all in the mid-tube pressure, because the
inlet pulse is a step with a 0.1 ms ramp (one time step at 100 us) and the
FSI iterations stop at 1e-6, which is the size of the differences being
measured at the finer steps. An order study on this case needs a smooth
pulse and tighter FSI and pressure tolerances (`orderStudy.py --smooth
--tight`, `orders_smooth.md`): with them Euler gives about 1 and backward
2.05-2.44 on every quantity, i.e. the coupling keeps the order of the
solvers. The `BDF` ddt scheme now supports moving meshes (solids4foam PR
#674), but on this case orders 3 and 4 meet two limits of the FSI
coupling, not of the schemes: `BDF 3`/`BDF 4` d2dt2 on the wall amplify
its undamped ring mode and the Robin-Neumann loop diverges (as the d2dt2
scheme's documentation warns); and with the wall kept at `BDF 2`, the
Robin-Neumann loop still diverges on the first step at which the fluid's
start-up ramp reaches order 3 or 4 (its pressure Robin coefficient is
derived for Euler/backward), whereas the same case with the
Dirichlet-Neumann IQN-ILS coupling runs through.

The IQN-ILS study (`orderStudy.py --smooth --tight --dn`,
`orders_smooth_dn.md`; `newMovingWallVelocity`/`zeroGradient` on the wall,
`IQNILS` with an initial relaxation of 0.05) settles what the FSI case can
show. Its interface residual levels off near 1e-7 at the finer steps (the
least-squares modes become collinear as the increments shrink), so that is
its tolerance. With it, fluid `BDF 3`/`BDF 4` with the coupling at order
3/4 and the wall at `BDF 2` (`BDF3s2`, `BDF4s2`) run through the whole
pulse at 17-21 FSI iterations per step at dt = 100, 50 and 25 us, which the
Robin-Neumann loop could not; the full `BDF 3`/`BDF 4` families (wall
d2dt2 at order 3/4) diverge in the FSI loop mid-pulse under IQN-ILS too,
confirming the ring-mode limit of the d2dt2 scheme. The observed orders
with IQN-ILS are 1-2 for every family (backward 1.78 at the port, 2.33 on
the wall; `BDF3s2` 1.3, `BDF4s2` 1.05 at the port): the differences between
successive steps are about 1e-2 of the signal and the same size whatever
the fluid scheme, i.e. they belong to the second- and first-order parts of
the FSI machinery (the interface velocity from the mesh displacement, the
mesh motion, the `BDF 2` wall), not to the fluid scheme or the 3D-1D
coupling. The conclusion for this case: the coupling preserves the order of
the FSI solver, which is 2; orders 3 and 4 of the coupling itself are
established on the pimpleFoam cases (`rigidSegment3D`, `dobroserdova2016`,
`womersleySegment`, `thoracicAorta37`, see `../orderStudy.py`).

## The wall's end at the interface: clamped, free, elastic support

`makeCases.py --c0 4.92 --springEnd` puts the cut end of the wall on an
elastic support (`solidSpringDashpot`): free axially and, in the plane of
the cut, the stiffness the continuing wall presents to its edge, taken as
that of a semi-infinite thin cylindrical shell displaced radially with its
edge moment free (Timoshenko), `Q = 2 D beta^3 u` per unit circumference
with `D = E h^3/(12 (1 - nu^2))` and `beta^4 = 3 (1 - nu^2)/(R^2 h^2)`,
divided by h: 1.04e7 Pa/m here (`shell_edge_stiffness()`). The three end
conditions against the full tube, all four runs made with the same build
of solids4foam (the installed library has no `solidSpringDashpot`; see
`results.md`, the `_lib2` runs and `full3D_lib2`):

| z [cm] | full tube peak [Pa] | clamped | free | spring | RMS(P - full) / pulse: clamped / free / spring |
|---|---|---|---|---|---|
| 2.5 | 1271 | 1303 | 1233 | 1233 | 3.0e-2 / 4.4e-2 / 5.1e-2 |
| 3.75 | 1133 | 1169 | 1108 | 1111 | 4.7e-2 / 2.4e-2 / 2.5e-2 |
| 5 (interface) | 1077 | 1104 | 1093 | 1098 | 4.4e-2 / 1.8e-2 / 2.2e-2 |
| 6.25 | 1045 | 1096 | 1085 | 1089 | 6.9e-2 / 3.6e-2 / 4.6e-2 |

(Wall displacement at 2.5 cm: RMS difference 3.6 / 3.2 / 4.1 um.) The
shell's edge stiffness is about a sixth of the ring's own hoop stiffness
`E h/R^2`, so the elastic support behaves nearly as the free end: the
overshoot at 3.75 cm and the interface is gone, the 3D part gains the
free end's oscillation. The full tube's peaks lie between the clamped and
the free/spring values at every section, so the condition that matches it
is stiffer than the continuing shell's static edge stiffness: the cut ring
in the full tube is held by the wall's inertia and the pressure acting on
the continuation, not by its static edge stiffness alone. The condition
to try next is Formaggia's continuity of area itself: the ring's radial
displacement prescribed from the 1D vessel's area at the interface (a
coupled Dirichlet condition on the solid), which needs the coupler to
drive a solid boundary and is not built.

## The wave returning from the 1D side

`makeCases.py --c0 4.92 --closedEnd --endTime 0.07` closes the far end:
the full tube gets a wall at 15 cm, the coupled run's 1D vessel a zero-flow
end (which reflects with +1 in pressure), and both run for 70 ms so that
the pulse, after leaving the 3D part at 10 ms, returns through the
interface at about 50 ms and reaches z = 1.25 cm at 58 ms. Against the
closed full tube (`coupled_c4p92_closed` against `full3D_closed`, same
build; the outgoing peak is the largest pressure before 30 ms, the
returning one after it: at 12.5 cm the full tube's largest pressure
overall, 974 Pa at 36.6 ms, is the wave reflected from its wall at 15 cm):

| z [cm] | outgoing peak: full / coupled [Pa] | returning peak: full / coupled [Pa], at [ms] | RMS(P - full) / pulse, 30-70 ms |
|---|---|---|---|
| 1.25 | 1352 / 1352 | 475 / 844 at 58.0 / 59.0 | 0.12 |
| 2.5 | 1271 / 1303 | 656 / 1054 at 56.5 / 56.4 | 0.10 |
| 3.75 | 1133 / 1169 | 680 / 1086 at 54.2 / 53.5 | 0.09 |
| 5 (interface) | 1077 / 1104 | 678 / 1124 at 51.9 / 50.1 | 0.12 |
| 10 | 902 / 1083 | | |
| 12.5 | 880 / 1077 | | |

The returning wave comes back through the interface into the 3D part at
the right time (within 1-2 ms of the full tube at every section) and
with the 1D vessel's amplitude: the full tube disperses and attenuates
the pulse on its way out and back through the thick wall (475-680 Pa on
return), the 1D vessel carries it almost undispersed (844-1124 Pa), and
the difference is already there on the outgoing wave at 10 and 12.5 cm
(1083 against 902, 1077 against 880 Pa) before anything has returned. So
the coupling passes a wave from 1D into 3D as it does from 3D into 1D;
what the comparison measures beyond that is the 1D model's lack of wall
inertia and dispersion, now counted twice.

## Results

See `results.md` and `plots/sections_P.png`, `sections_Q.png`,
`displacement.png`. Timings on 2026-11-05: `full3D` 14.5 min on 4 ranks,
each coupled run 11.5 min serial (3 FSI iterations per step, as the
tutorial; at dt = 100 us it needs 40 and is no cheaper).

- **The 3D part.** Until the wave reflected at the interface returns, the
  coupled tube is the full tube: the radial displacement at z = 2.5 cm agrees
  to 3.3-3.5 um RMS of its 150 um peak over the 20 ms, the mean pressure to
  2.5-4 % of the pulse RMS at z = 1.25-5 cm.
- **The wall match.** The pulse's half-peak travels the full 3D tube at
  4.92 m/s, below the 5.25 m/s of the static plane-strain match: the thick
  wall's inertia slows it and disperses it (the peak pressure falls from
  1388 Pa at z = 1.25 cm to 1093 Pa at 5 cm; the 1D vessel alone keeps
  1330 Pa). The 1D vessel at 5.25 m/s is the stiffer line and reflects about
  5 % of the pulse: an overshoot of +70 Pa at the interface and a return of
  about +80 Pa at z = 1.25 cm near 18 ms. Matching the 1D to the measured
  speed (`coupled_c4p92`) halves the overshoot (1124 against 1093 Pa) and
  brings the arrival times beyond the interface to within 0.2 ms of the full
  tube; the RMS differences in the 1D part drop from 9-10 % to 5-6 %.
- **What remains.** Beyond the interface the 1D vessel carries the pulse
  undispersed, so at z = 6.25-7.5 cm it is sharper and earlier than in the
  full tube: a difference between the models, as Formaggia et al. found. In
  the 3D part the overshoot at z = 3.75 cm (1171 against 1133 Pa) is the
  same for both 1D wave speeds, so it is not the line's impedance but the
  structural end condition at the interface: the wall is clamped there in
  the coupled tube and free to breathe in the full one, which Formaggia et
  al. address through their "continuity of area" condition.
  `makeCases.py --c0 4.92 --freeEnd` (`coupled_c4p92_free`) frees that end
  (zero traction): the overshoots at z = 3.75 and 5 cm turn into small
  undershoots (1107 and 1110 against 1133 and 1093 Pa) and the RMS
  differences beyond the interface halve again, to 3.4-3.7 %, but a free
  ring is softer than the continuing wall and the 3D part gains an
  oscillation of its own (+200 Pa at z = 2.5 cm near 13 ms; 4.7 % RMS at
  z = 1.25 cm). Clamped and free bracket the full tube; the condition in
  between is an elastic support with the stiffness of the continuing wall,
  which is what a "continuity of area" condition supplies.
