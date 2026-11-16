# A 3D segment in the 37-artery network: `thoracicAorta37`

Prepared by David Kelly

---

## Tutorial Aims

- Embed a 3D CFD segment in a validated physiological 1D network at
  physiological amplitude: the thoracic aorta II of the 37-artery in vitro
  model of Matthys et al. (Boileau et al. 2015, section 4.5), computed by
  pimpleFoam as a rigid no-slip cone, with the other 36 vessels, the
  junctions and the Windkessel terminals in 1D.
- Measure what the finite-element studies of this configuration (Blanco et
  al. 2013, Malossi et al. 2013, Passerini et al. 2009) do not report: the
  effect of the rigid segment on the whole network against its exact
  inviscid limit (an inertance node), the coupling's cost at physiological
  amplitude, and the comparison with the in vitro measurements.

---

## Case Overview

Segment s17 (72 mm, radius 11 to 9.26 mm, the `ThoracicAortaII` measurement
site at its midpoint, junctions at both ends) is split into an 8 mm 1D stub
`s17a`, a 56 mm middle and an 8 mm stub `s17b`; `makeCases.py` writes the
network of the Boileau case with that split and:

| run | the middle |
| --- | --- |
| `stubs1D` | a 1D vessel `s17m`: must equal the original network; 18 cycles from rest (64 s) and the periodic state at 14.886 s |
| `stubsInertance` | an `inertance` node, `L = rho l/(pi R0 R1)` (a cone's inertance): the exact inviscid rigid limit; restarted from that state, two cycles (7 s) |
| `rigid3D` | the pimpleFoam cone (O-grid, 8 x 8 core and four outer blocks of 8 cells graded 4:1 to the wall, 28 axial cells; no-slip, mu = 2.5 mPa s; `backward`, dt = 1 ms; implicit coupling of order 2); restarted from that state with the 3D as plug flow at the state's flow, two cycles (4.5 min) |

The original network's results and the in vitro measurements are read from
the Boileau case. Face zones at the 3D mid-plane (the measurement site) and
2 mm from each port record the mean pressure and the flow; the coupler's
log holds every 1D probe of the network.

**Restart.** `circulation1DCoupler` reads `circulation1DState` at the start
time, as `circulation1DFoam` does, so the coupled run resumes the periodic
state of `stubs1D` (vessels of the same names and meshes; the vessel the 3D
domain replaces is ignored, and a node the stored network did not have, the
inertance, keeps its initial state). The volume-drift line of the restarted
`stubsInertance` run reports the volume of the vessel it no longer has, not
a drift. The coupling's history is seeded with the state's flow through the
segment (`seed`), and the 3D velocity is initialised as plug flow at that
flow.

```bash
./Allrun                  # stubs1D (1 min), then stubsInertance and rigid3D (5 min), then analyse.py
```

---

## Results

See `results.md` and `plots/sites.png`.

| site | run | peak P [kPa] | peak Q [ml/s] | RMS dP vs original / pulse | RMS dQ vs original / peak Q | vs in vitro: RMS dP / dQ |
|---|---|---|---|---|---|---|
| aortic arch II (upstream) | original 1D | 15.01 | 179 | - | - | 0.084 / 0.134 |
| | stubs1D | 15.01 | 179 | 4e-7 | 2e-6 | 0.084 / 0.134 |
| | inertance | 15.38 | 184 | 4.4e-2 | 5.0e-2 | 0.086 / 0.123 |
| | rigid3D | 15.38 | 183 | 4.5e-2 | 5.1e-2 | 0.086 / 0.124 |
| thoracic aorta II (3D mid-plane) | original 1D | 15.21 | 128 | - | - | 0.118 / 0.291 |
| | inertance | 15.54 | 123 | 4.3e-2 | 3.9e-2 | 0.119 / 0.285 |
| | rigid3D | 15.60 | 121 | 4.4e-2 | 4.2e-2 | 0.122 / 0.285 |
| right iliac-femoral II (downstream) | original 1D | 16.29 | 12.4 | - | - | 0.299 / 0.378 |
| | inertance | 16.80 | 13.4 | 4.5e-2 | 5.1e-2 | 0.309 / 0.412 |
| | rigid3D | 16.79 | 13.1 | 4.5e-2 | 5.1e-2 | 0.311 / 0.412 |

rigid3D against the inertance limit (RMS, of the pulse and of the peak
flow): 2.0e-3 / 3.7e-3 at the arch, 1.5e-2 / 6.5e-3 at the mid-plane,
4.7e-3 / 8.1e-3 at the femoral.

- **Splitting the vessel** changes nothing (1e-6).
- **Making 56 mm of the aorta rigid** changes the pressure everywhere in the
  network by 4.4 % of the pulse (RMS) and the flow by 4-5 % of the peak:
  +0.4 kPa on the systolic peak at the arch and the thoracic aorta, +0.5 kPa
  at the femoral, with the ripples of the reflected harmonics (the segment's
  tau = L/(2 Zc) is 5.6 ms, so it reflects 4 % of the first harmonic and 18 %
  of the fifth). That is the physics of a rigid 3D segment in a compliant
  network, and it is the same in the exact inertance limit and in the
  coupled 3D run.
- **The coupled 3D run** differs from that limit by 0.2-0.5 % of the pulse
  upstream and downstream and by 1.5 % at its own mid-plane, where the
  viscous pressure loss and the taper of the cone (the inertance node uses
  the equivalent radius) act: the coupling adds nothing visible to the
  physics of the segment.
- **Against the in vitro data** the rigid segment moves the Boileau errors
  by less than 0.01, in either direction: the 1D model's own error (8-30 %
  RMS at these sites in Boileau's Table X) dominates.
- **Cost**: 5.1 outer iterations per step (10 network re-integrations), 273 s
  for two cycles of 827 steps on one core.

Planned on this case: static against total pressure continuity at the
ports, the order in dt at physiological amplitude, and the compliant segment
with solids4foam FSI (silicone wall, E = 1.2 MPa, h = 0.34 mm, as the in
vitro tube).

## Static against total pressure at the ports

The coupler imposes the static pressure of the vessel end on the patch by
default; `coupling { pressure total; kineticEnergy mean; }` makes the total
pressure continuous instead, `p = P1D + rho U1D^2/2 - pdyn3D` with U1D the
velocity of the vessel's boundary state and pdyn3D the patch's dynamic
pressure from the sampled flux: `mean` (rho U3D^2/2 of the mean normal
velocity), `area` (the area average of the local dynamic pressure) or
`flux` (the energy-flux form rho <u^3>/(2 <u>), which conserves the
mechanical energy flux but is singular at flow reversal, where <u^3> stays
finite as <u> passes through zero: in this case it spikes to 475 Pa, 10 %
of the pulse, twice per cycle). The `inertance` node has the same choice
with the mean velocities (`pressure total`). `makeCases.py --restart
--pressure total --kineticEnergy mean|area|flux` writes
`rigid3D_total_<form>` and `stubsInertance_total`; `results.md` compares
them with the static runs over the last cycle:

| site | total (mean) vs static: RMS dP / pulse, max dP / pulse, RMS dQ / peak Q | area | 1D limit, total vs static |
|---|---|---|---|
| AorticArchII | 6.0e-4 / 1.9e-3 / 8.5e-4 | 6.1e-4 / 1.9e-3 / 8.8e-4 | 6.3e-4 / 2.1e-3 / 1.0e-3 |
| ThoracicAortaII (3D mid-plane) | 2.9e-3 / 8.9e-3 / 1.6e-3 | 3.3e-3 / 9.4e-3 / 1.6e-3 | 5.5e-4 / 1.6e-3 / 1.8e-3 |
| RIliacFemoralII | 1.3e-3 / 3.9e-3 / 2.0e-3 | 1.3e-3 / 4.0e-3 / 2.1e-3 | 1.4e-3 / 4.7e-3 / 2.4e-3 |

At physiological amplitude in a vessel of this size the choice moves the
network by 0.06-0.2 % of the pulse (RMS) and the 3D segment's own pressure
by 0.3 % (42 Pa at peak flow, where rho U^2/2 = 75 Pa at the stub ends);
the mean and area forms agree to 5 % of that. The coupled run with total
pressure follows the total-pressure inertance limit as closely as the
static one follows the static limit.

## Order in time

`../orderStudy.py --case thoracicAorta37` runs the coupled case from rest
over the first systole (0-0.3 s) at dt = 2, 1, 0.5 and 0.25 ms, with the
fluid ddt scheme and the implicit coupling order both at q = 1 (`Euler`), 2
(`backward`), 3 (`BDF 3`) and 4 (`BDF 4`), the network at CFL 0.1 and the
PIMPLE loop to a pressure residual of 1e-10. For this study the network
starts in equilibrium at the terminal pressure (432.6 Pa) and the inflow
rises under a sin^2 ramp of 50 ms: the tutorial's own start (P = 0 against
the terminals, the inflow's Fourier series at its t = 0 value) launches
fronts one 1D cell wide, 0.3 ms at the port, i.e. features at the scale of
dt (`orders_impulsiveStart.md` has that version; the orders below are the
same either way). Observed orders by self-convergence (`orders.md`,
`plots/orders.png`), port flow / port pressure:

| dt triple [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 2 / 1 / 0.5 | 0.96 / 0.88 | 2.00 / 2.01 | 2.34 / 2.38 | 1.08 / 1.30 |
| 1 / 0.5 / 0.25 | 0.98 / 0.95 | 1.95 / 1.89 | 0.61 / 0.54 | 0.21 / 0.16 |

Euler and backward converge at their orders. BDF 3 and BDF 4 do not reach
theirs on this case: their differences between successive steps stop
falling at about 3e-6 of the signal (RMS; 3e-5 of the peak flow at most),
spread uniformly over the active part of the window and of the same size
for both schemes, i.e. a floor rather than an order. The finest BDF 3 and
BDF 4 runs agree to 3.3e-6 of the flow and with the finest backward run to
1.6e-5. The floor is not the 1D network: the 1D limiter off and SSPRK3
(`orders_none_SSPRK3` variant, since removed) and the network at CFL 0.02
instead of 0.1 (`orders_cfl0.02.md`) change the runs by less than 1e-7. It
is not the start (the smooth start above leaves it unchanged), not the
PIMPLE convergence (the pressure residual reaches 3e-12 every step), and
not the coupling order (the coupler runs at order 4).

**The floor traced** (`floorStudy.py`, `floor.md`, `plots/floor.png`): the
BDF 4 runs at dt = 1, 0.5 and 0.25 ms repeated with one change at a time,
the difference between successive steps (RMS, relative to the signal) of
the port flow and the order it gives:

| variant | Q: e(1, 0.5) | e(0.5, 0.25) | order | P order |
|---|---|---|---|---|
| baseline (cone, tapered stubs) | 2.18e-05 | 1.88e-05 | 0.21 | 0.16 |
| linear solvers to 1e-15, p residual 1e-12 | 2.18e-05 | 1.88e-05 | 0.21 | 0.16 |
| one non-orthogonal corrector | 2.18e-05 | 1.88e-05 | 0.21 | 0.16 |
| uniform stubs (the cone unchanged) | 2.18e-05 | 1.89e-05 | 0.21 | 0.16 |
| `ddtCorr no` (flux correction off) | 1.48e-05 | 1.04e-05 | 0.50 | 0.70 |
| cylinder of the mean radius instead of the cone | 2.76e-06 | 1.95e-07 | 3.82 | 3.91 |

The floor is the 3D side's, on the cone: with the same network, stubs,
coupling and schemes, the 3D segment as a cylinder (an orthogonal mesh)
converges at 3.8-3.9, as the uniform-vessel cases do. The tapered stubs,
the junctions, the solver tolerances and the non-orthogonal correction
are not involved. On the cone, pimpleFoam's flux correction
(`fvc::ddtCorr`, whose coupling coefficient switches on the local
mismatch between the interpolated velocity and the flux, which is zero on
an orthogonal mesh in plug or developed flow and not on a cone) carries
about half of it; the rest is in the non-orthogonal mesh's own treatment
in the segregated algorithm. Both are below the second-order error, which
is why `backward` converges cleanly, and above the fourth-order one. For
the physics measured here (differences of 4 % of the pulse) a floor at
3e-6 is of no consequence; the claim to make is that the coupling takes
the 3D scheme's order on the network when the 3D discretisation has that
order, which the cylinder shows, and that a non-orthogonal 3D mesh in
pimpleFoam does not go beyond second order in time, whatever the ddt
scheme.
