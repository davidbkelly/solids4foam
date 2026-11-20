# A cross-code case: the pressure Gaussian pulse of Businaro (2026): `businaro2026`

Prepared by David Kelly

---

## Tutorial Aims

- Repeat, with the solids4foam FSI solver and the circulation1D coupler, the
  one coupled 3D-1D test of Businaro's thesis (Politecnico di Milano, 2026,
  Section 7.4 and Table F.1), computed there with lifex (deal.II): a
  compliant FSI tube continued by a 1D vessel, driven by a Gaussian pressure
  pulse. It is the only case of this kind with a second, independent code.
- Compare the coupled runs with a full FSI tube, with the 1D vessel alone,
  and with the interface mismatch Businaro reports for his explicit coupling.

---

## Case Overview

Businaro's set-up (Table F.1): a tube of radius R = 1 cm and wall thickness
1 mm, the wall Saint Venant-Kirchhoff (rho_s = 1000 kg/m3, E = 700 kPa,
nu = 0.3), blood rho = 1060 kg/m3, mu = 4 mPa s; the inlet pressure
P(t) = 2 kPa exp(-1e4 (t - 0.05)^2) (his Eq. 7.11); the wall's ends fixed
axially and free radially, at the inlet and at the interface (his Eq. 5.12;
here `solidSymmetry` end faces); BDF2 at dt = 0.1 ms over 0.1 s; 1 mm axial
cells. The Womersley number of the pulse is 62 (his 62.0), so the fluid
mesh is graded towards the wall (wall cell 0.125 mm against a Stokes layer
of 0.23 mm). The 3D cases are quarter tubes.

The coupled tube is 10 cm of FSI tube and a 20 cm 1D vessel with an
absorbing end. The 1D vessel is Businaro's: his tube law (Eq. 2.3),
beta = sqrt(pi) h E/(1 - nu^2), so c0 = 6.02 m/s (the sqrtElastic law with
E1D = 3E/(4(1 - nu^2))), and the friction of his blunt profile, KR =
2 pi (gamma + 2) nu with gamma = 9.

The runs:

- `full3D`, the reference: the FSI tube 45 cm long, p = 0 at its far end;
  nothing reflected there reaches z <= 30 cm by t = 0.1 s.
- `full3D_Rt` (`--withRt`), Businaro's own reference: the tube 30 cm long
  with his terminal resistance Rt = 2.0494e7 Pa s/m3 at the outlet (through
  a 1 cm 1D stub ending in a `resistance` node).
- `coupled`: 10 cm FSI + 20 cm 1D, implicit coupling of order 2.
- `coupled_explicit`: the same with explicit coupling of order 1, the 1D
  vessel advanced with the 3D flow of the previous step and the 3D tube
  with the 1D pressure at the end of the step, as Businaro's algorithm (his
  Eq. 5.21).
- `only1D`: the 30 cm tube as a 1D vessel alone.

Differences from Businaro's set-up, all deliberate:

- **The reference tube is 45 cm long** with p = 0 at its end, instead of
  30 cm with the resistance Rt: Businaro notes that his resistance outlet,
  applied with the flow of the previous step, sends spurious reflections
  back into the tube. His version is `full3D_Rt`, with the resistance
  applied implicitly.
- **alpha = 1** in the 1D vessel (his alpha = 11/10 for gamma = 9; the
  circulation1D solver needs alpha = 1), with his friction coefficient.
- **The interface conditions**: the 3D outlet takes the 1D pressure as a
  uniform static pressure and the vessel the 3D flow (Businaro: the 1D
  pressure as the normal stress, the vessel the 3D state through its
  incoming characteristic). His explicit scheme is reproduced by
  `coupled_explicit`; `coupled` iterates the coupling to convergence.
- **Two solver settings for a flow at rest.** The pulse is 1e-8 Pa at
  t = 0, so relative residuals have no scale: the wall's Newton iteration
  has an absolute tolerance `snes_atol 1e-12` (at rest its residual is
  round-off, 2e-15; under load the iteration has already reduced it by 1e-9
  when it stops), and the Robin FSI residuals are scaled by
  `robinPressureReference` (the 2 kPa pulse, as kinematic pressure) and
  `robinFluxReferenceVelocity` (0.3 m/s, the peak flow over the area) in
  `constant/fsiProperties`. The tolerances themselves are the defaults.

**Measuring at the sections.** The pressure and flow at the sections inside
the 3D tube come from a coded function object, `sectionFields`, in each
run's `system/controlDict`. The solver leaves `phi` relative to the mesh,
and the sections move axially with the wall, so the flow is the sum of
`phiAbsolute = phi + meshPhi`. On a processor boundary `p` holds the
neighbour cell's value, so the pressure is the area average of `pSection`,
the mean of the two cells of every face. Without these the flows at the
sections were 10 % off and a section on a processor boundary read the
pressure half a cell away.

---

## Running the Case

```bash
NPROC=16 ./Allrun --withRt  # about 5 h; the 30 cm Rt tube is serial
python3 analyse.py          # results.md and plots/ again
```

`Allrun` writes the runs with `makeCases.py`, meshes them, runs the 1D
vessel alone, the coupled runs and `full3D_Rt` serially alongside the 45 cm
reference on `NPROC` ranks, and calls `analyse.py`. A finished run is not
run again. `./Allclean` removes the runs, the plots and the results.

---

## Expected Results

See `results.md`, `plots/snapshots.png` (pressure and flow along the tube at
t = 0.06, 0.08 and 0.10 s, as Businaro's Figure 7.12), `plots/sections.png`
and `plots/interface.png` (the interface mismatch, his Figure 7.9).

**Interface mismatch** (Businaro's Eq. 7.12, the largest value over the
window; max Q is the largest flow through the interface):

| run | e_P / 2 kPa | e_Q / max Q | max Q [m3/s] |
| --- | --- | --- | --- |
| `coupled` | 2.6e-9 | 2.4e-8 | 1.0283e-4 |
| `coupled_explicit` | 2.6e-9 | 9.5e-3 | 1.0282e-4 |
| Businaro (2026), explicit | 2.2e-4 | 5.6e-2 | 1.0330e-4 |

- e_P is zero by construction here: the 3D outlet takes the 1D pressure as
  its boundary value, and what remains is the rounding of the logged
  values. In Businaro's scheme the pressure is a normal stress and e_P is
  a measured difference.
- The implicit coupling converges the flow to 2e-8. The explicit one lags
  the flow by a step, 0.95 % of the peak against Businaro's 5.6 %.
- The peak flow through the interface agrees with lifex to 0.5 %
  (1.028e-4 against 1.033e-4 m3/s).

**Against the full FSI tube.** RMS pressure difference over the window,
relative to the 2 kPa pulse, at z in cm; the interface is at z = 10 cm
(`results.md` also has the maximum, the flow and the peak ratios):

| run | z = 5 | z = 10 | z = 15 | z = 20 | z = 30 |
| --- | --- | --- | --- | --- | --- |
| `full3D_Rt` | 0.99 % | 1.1 % | 0.97 % | 0.62 % | 1.5 % |
| `coupled` | 0.90 % | 2.1 % | 4.4 % | 6.5 % | 10.1 % |
| `coupled_explicit` | 1.0 % | 2.0 % | 4.4 % | 6.5 % | 10.1 % |
| `only1D` | 2.4 % | 5.2 % | 7.9 % | 10.0 % | 12.8 % |

- In the 3D part the coupled runs are within 0.9-2.1 % of the full tube,
  less than half the difference of the 1D vessel alone.
- Beyond the interface the difference grows with distance, because the
  1D vessel carries the pulse faster than the FSI tube: the half-peak
  arrivals give 6.23 m/s in the 1D vessel (c0 = 6.02 m/s from Businaro's
  tube law, plus the pulse's own amplitude) and 5.67 m/s in the full tube.
  Over 20 cm of vessel the coupled pulse arrives 3.2 ms early, against a
  pulse 17 ms wide at half its peak (`plots/snapshots.png`). The 1D vessel
  alone has the same speed, so this is the 1D model, not the coupling.
- The same mismatch of impedance shows at the interface: the coupled
  interface pressure peaks 2.7 % above the full tube's at z = 10 cm and
  the flow 3.6 % below (1.028e-4 against 1.067e-4 m3/s).
- Beyond the interface explicit and implicit coupling differ from the
  full tube by the same amount, to 0.01 % of the pulse.
- Businaro's own reference, `full3D_Rt`, is within 1.1 % of the 45 cm
  tube up to z = 25 cm. At its outlet the pressure peaks 6.5 % higher:
  Rt is close to the 1D vessel's characteristic impedance, not the 3D
  tube's, so the outlet reflects part of the wave.
