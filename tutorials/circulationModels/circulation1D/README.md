# 1D arterial flow benchmarks: `circulation1D`

---

Prepared by David Kelly

---

## Tutorial Aims

- Verify and validate `circulation1DFoam` against the benchmark of Boileau et
  al. (2015), *A benchmark study of numerical schemes for one-dimensional
  arterial blood flow modelling*, Int. J. Numer. Meth. Biomed. Engng. 31,
  e02732, with its published data for six schemes and a 3D model.
- Measure the order of accuracy of the scheme.

---

## Cases

| Case | Boileau | Reference | What it checks |
|---|---|---|---|
| `boileau2015/singlePulseInviscid` | 4.1, Table II | linear theory, six schemes | a 20 Pa pulse travels 10 m at `c0` with < 1% peak decay |
| `boileau2015/singlePulseViscous` | 4.1, Eq. 15 | Eq. 15, six schemes | the peak damps as `exp(-(zeta + 2) pi mu x/(rho c0 A0))` |
| `singlePulseOrder` | (none) | linear solution | order 2 in space and time |
| `boileau2015/commonCarotid` | 4.2, Tables III-IV | 3D, six schemes | midpoint P, Q, dr and inlet-outlet dP |
| `boileau2015/upperThoracicAorta` | 4.3, Tables V-VI | 3D, six schemes | the same |
| `junctionReflection` | (none) | linear theory | reflection and transmission at a 1 -> 2 junction |
| `boileau2015/aorticBifurcation` | 4.4, Tables VII-VIII | 3D, six schemes | a 3-vessel junction; P, Q, dr at three sites (Fig. 5) |
| `boileau2015/thirtySevenArtery` | 4.5, Tables IX-X | in vitro, six schemes | a 37-segment network, resistance terminals (Figs. 6-9) |
| `fluidAtRest` | (none) | rest | a tapered network at uniform pressure stays at rest |
| `taperedOrder` | (none) | finer meshes | order 2 in a tapered vessel |
| `boileau2015/ADAN56` | 4.6, Tables XI-XIII | six schemes | a 77-vessel network, tapering, Eq. 17 walls, Windkessels (Figs. 11-16) |
| `singlePulseHighOrder` | (none) | finest run | the k-exact reconstruction and RK4: order 2, 3, 4 in space and time |

`common/` holds the shared analysis: `boileauMetrics.py` (Eqs. 10-13, the
published tables and the figures; `python3 boileauMetrics.py selftest`
recomputes Tables IV, VI, VIII and X from the six schemes' data) and
`singlePulse.py`.

Each case has `Allrun`, `Allclean`, `constant/circulation1DProperties`,
`system/controlDict`, `referenceData/` (from the paper's data supplement:
the 3D data unchanged, the six schemes trimmed to the compared cycle or
thinned along x) and `plot.py`, which writes `metrics.md` or `peaks.md` and the
figures in `plots/`.

---

## Running

The cases need OpenFOAM.com, the `circulationModels` module and
`circulation1DFoam` built, and a Python with numpy and matplotlib:

```bash
./Allrun                              # every case, about 8 minutes in all
PYTHON=~/venv/bin/python ./Allrun     # another Python
cd boileau2015/commonCarotid && ./Allrun     # one case
./Allclean
```

To see a case's network, draw it with `circulation1DTree` (in
`src/circulationModels/tools`, needs matplotlib), for example
`circulation1DTree boileau2015/ADAN56`, or `circulation1DTree -o tree.png
boileau2015/ADAN56` for a
file.

---

## Set-up notes

- The physiological inflows are the Fourier series of the data supplement,
  written as `expression` `Function1`s with `arg()` the time, from `t = 0` as
  in the paper. The vessels start from rest at `P = 0` (Tables III and V), and
  run 9 (carotid) and 20 (aorta) cycles.
- The Windkessels are `lumped0D` nodes with the `windkessel3` preset of
  circulation0D: `Rch = R1`, `C`, `R = R2`, `Pout = 0`.
- **Time alignment of the 3D data.** The supplement does not say how the 3D
  cycle lines up with the 1D runs, and its README windows do not reproduce the
  tables. The windows that do are `3D(tau) <-> 1D(k T + 0.056 + tau)` for the
  carotid and `1D(k T + 0.002 + tau)` for the aorta (`selftest` reproduces the
  LCG, FEM, FVM and McC columns to 0.01%). The metrics use the last whole
  cycle of the run with that phase, and report the cycle-to-cycle change.
- **Resolution.** The scheme is second order, against fourth and fifth order
  for the paper's FVM and DCG at `dx = 1-2 cm`. The single pulse (standard
  deviation 4.4 cm along x) needs `dx = 1.25 mm` and the unlimited
  reconstruction to keep its peak decay below 1%. Its peak decay at x = 9 m:

  | dx | `none`, CFL 0.9 | `none`, CFL 0.4 | `MC`, CFL 0.4 |
  |---|---|---|---|
  | 10 mm | 34.3% | | |
  | 5 mm | 14.6% | | |
  | 2.5 mm | 4.6% | | 6.4% |
  | 1.25 mm | 0.93% | 0.58% | 2.0% |
  | 0.625 mm | 0.21% | | 0.94% |

  The carotid and aorta waveforms are long (`c T` of several metres), and
  their metrics do not change, to two decimals, from `dx = 1 mm` to 0.5 mm or
  from `MC` to `none`.
- **CFL.** A limited reconstruction is TVD for `CFL <= 0.5` only. At 0.9 with
  `MC` the single pulse decays *more* on a finer mesh (3.4% at 8000 cells,
  11.8% at 16000), as the limiter's clipping accumulates over more steps; the
  solver warns.
- **37-artery data.** The supplement's `mesh.txt` labels the terminal
  resistances Pa s/m3, but they are in 1e9 Pa s/m3 (as labelled, the mean
  pressure would be 0.4 kPa; scaled, 12.1 kPa against 12.7 kPa measured). It
  names no segments: each of the eight measured sites is the segment whose
  midpoint waveforms match the published FVM ones (relative RMS 0.2-0.4%,
  against 4% or more for every other segment). The in vitro windows are
  `1D(tini - 0.001 + tau)`, with `tini` of `time_shifts.txt`.
- **37-artery tapering.** 22 segments taper and the paper's model has
  `A0(x)`, the radius linear along a segment. `makeNetwork.py` writes them as
  tapered vessels. With `--taper pieces` it writes them instead as uniform
  pieces with radius steps of at most 2% joined by two-vessel junctions (the
  approach of phase 2): the metrics agree within 0.06. They change by at most
  0.11 with `dx = 1 mm`.
- **ADAN56 data.** `mesh.txt` is in cgs (converted on writing the dictionary)
  and has no wall thickness: it follows Boileau Eq. 17, whose constants b and d
  are negative (the minus signs are easy to lose from the paper's text). The
  inflow is a table with period 1 s, linearly interpolated. The data
  supplement lists the times of the pressure along the aorta as 9.1, 9.2, 9.3,
  9.5, 9.7 and "1.0" s: the last is the end of the cycle (t mod T = 0); its
  last three columns are the systolic, diastolic and mean pressure. The site
  `RightPOstInteross` of the supplement is `RightPostInteross` here.
- **Rest state.** With properties varying along a vessel, the rest state is a
  uniform pressure with no flow; it is `A = Ad(x)` only if `Pd + Pext` is
  uniform. `fluidAtRest` tests both a rest at `Pd` and one below it.
- **Order test amplitude.** At the paper's amplitude the pulse is not linear
  enough for a convergence study: it steepens by a few millimetres over 5 m
  (the peak travels 0.06% faster than `c0`), and its error against the linear
  solution levels off at 3% of the peak. `linear/` scales the inflow by 1e-3,
  `full/` keeps it, to show the plateau.

---

## Results (phases 1 to 3)

- `boileau2015/singlePulseInviscid`: peak decay 0.64% at x = 9 m; wave speed 6.178 m/s
  (`c0 = 6.172` m/s: the peak also travels at `U`).
- `boileau2015/singlePulseViscous`: peak within 0.43% of Eq. 15 to x = 9 m.
- `singlePulseOrder`: least-squares order over the four finest meshes 1.98
  (L1), 1.98 (L2), 1.99 (Linf) for `linear`; 0.5 for `full` (plateau).
- `boileau2015/commonCarotid`: all 16 errors of Table IV in the range of the six schemes
  (one below it).
- `boileau2015/upperThoracicAorta`: 15 of 16 errors of Table VI in or below the range;
  `E^MAX_P` is 3.57% against 3.18-3.55%, the same on a finer mesh or without
  limiter.
- `junctionReflection`: reflected and transmitted pressure within 0.0005 of
  linear theory for R = +0.41, 0 and -1/3.
- `boileau2015/aorticBifurcation`: all 36 errors of Table VIII in the range of the six
  schemes (two below it).
- `boileau2015/thirtySevenArtery`: 60 of the 64 errors of Table X in the range; the
  others are `E^DIAS_Q` at aortic arch II (18.45 against 16.93-17.99), thoracic
  aorta II (40.01 against 35.47-39.65) and the right anterior tibial (15.24
  against 13.25-15.09), and `E^RMS_P` there (3.66 against 3.21-3.43). The
  published anterior tibial row does not reproduce from the published data
  (its window runs past the end of the 1D runs).
- `fluidAtRest`: at `Pd` the network stays exactly at rest (Q and P - Pd
  zero); below `Pd`, |Q| stays below 2e-16 of the flow scale and |P - P_rest|
  below 5e-15 of the pressure, over 6000 steps.
- `taperedOrder`: orders 1.95, 1.99, 2.00 and 1.96 from successive meshes of
  750 to 24000 cells.
- `boileau2015/ADAN56`: at the twelve sites of Figs. 11-14, solids4foam differs from each
  of the five consistent schemes (DCG, LCG, FEM, FVM, McC) by less than the
  largest difference between two of them (relative RMS of P 0.16-0.64%, of Q
  0.11-0.57%). Along the aorta it is within 0.73% of FVM at the instants of
  Fig. 15 and within 0.26% for the systolic, mean and diastolic pressure of
  Fig. 16, where the other four schemes differ from FVM by up to 2.0% and
  0.51%.
- The volume balance closes to round-off (1e-14 of the volume) in every case.

## Results (higher order)

- `singlePulseHighOrder`: the single pulse at full amplitude, refined at CFL
  0.5 from 250 to 8000 cells, each run against the 32000-cell run of its set.
  Order between the two finest meshes (Q, L1 / L2 / max): MUSCL + SSPRK2 2.05 /
  2.04 / 2.04; kExact 1 + SSPRK2 the same; kExact 2 + SSPRK3 3.00 / 2.99 / 2.98;
  kExact 3 + RK4 4.01 / 4.01 / 4.00, and with friction 4.01 / 4.01 / 4.01. The
  finest runs of the different schemes differ by the coarser one's own error
  (1.03e-3 and 1.05e-5 of the peak), so they converge to the same solution.
  At 8000 cells (dx = 1.25 mm) the error (L2) is 1.2e-3 at second order and
  9.6e-7 at fourth, for 2.7 times the cost per run.
- The rigid segment with the inertance node (kExact 3 + RK4, 200 to 3200 cells
  per vessel, against 12800): the pipe flow, both pipe pressures and the
  profiles converge at 4.00-4.01, so the nodes at the vessel ends keep the
  order. At the 1e-3 amplitude of the tutorial the errors stop falling at
  about 5e-7 of the peak: round-off in the momentum flux at rest, divided by
  dx, against so small a pulse; at full amplitude they do not.

The 3D-1D coupling cases (the rigid-segment benchmark and the literature
tests) are in [../coupling3D1D](../coupling3D1D).
