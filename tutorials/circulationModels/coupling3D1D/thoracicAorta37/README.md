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
