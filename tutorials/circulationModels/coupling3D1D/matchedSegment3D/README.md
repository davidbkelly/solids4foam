# Matched segment, 3D-1D: `matchedSegment3D`

Prepared by David Kelly

---

## Tutorial Aims

- Test the 3D solver carrying a pressure wave along its own length, which
  no other coupled case does: in `rigidSegment3D`, `dobroserdova2016` and
  `womersleySegment` the 3D domain is rigid and incompressible, so the wave
  crosses it instantly, and in `formaggia2007` the elastic wall carries it.
- Couple a compressible 3D solver (stock `rhoPimpleFoam`) to 1D vessels: the
  coupler samples the patch flux as `phi/rho` and imposes the pressure in Pa.
- Measure the order of the coupled system along the direction of the wave,
  against an exact solution in which nothing is reflected.

---

## Case Overview

The two vessels of `rigidSegment3D` (1 m each, the single-pulse vessel of
Boileau et al. 2015: `c0` = 6.17 m/s) are joined by a rigid 3D pipe 20 cm
long, radius 1 cm, slip walls, inviscid. The pipe holds a compressible fluid
whose sound speed and impedance are the vessels':

- `adiabaticPerfectFluid` with `gamma` = 1 and `p0` = 0 gives
  `rho = rho0 + p/c^2` with `c^2 = B/rho0`; with `rho0` = 1050 and
  `B = rho0 c0^2` = 40000 Pa the sound speed is exactly `c0`, independent of
  the temperature;
- the patches are 32-sided polygons with the area `pi R^2` exactly, so the
  3D impedance `rho c/A` equals the vessels' `Zc`.

The whole system is then one uniform lossless line. A Gaussian pulse (Boileau
Eq. 14, at amplitude 1e-3) enters the proximal vessel, passes through the
pipe with **no reflection**, and leaves into the distal vessel delayed by
`0.2/c0` = 32.4 ms. The vessels solve the linearised 1D equations (k-exact 3
with RK4, 1600 cells each), so this exact solution holds at any amplitude and
their own error, about 1e-5 of the peak, stays below the 3D one. The 3D mesh
is `rigidSegment3D`'s O-grid with one cell across the outer ring: 96 cells
per section, 40 along the pipe (`dz` = 5 mm), `dt` = 0.4 ms (acoustic
Courant number `c0 dt/dz` = 0.49), BDF2, implicit coupling of order 2.

Two practical points:

- The run has two phases. Nothing reaches the pipe before 0.16 s, even at
  1e-12 of the peak, so the first 0.12 s run with one outer corrector per
  step, and the run continues from the written state (the coupler restarts
  from its checkpoint) with 20.
- The implicit coupling contracts by about the acoustic Courant number per
  outer corrector (0.55 per iteration here), the 3D side answering a change
  of its port pressure through one cell's inertia within a step. Twenty
  correctors leave 1e-5 of a step's change of the port flow. A residual
  control is not used: while the pipe is at rest the normalised pressure
  residual is round-off and never meets a tolerance.

---

## Running the Case

```bash
./Allrun                    # one case, a few minutes, checked against exact
./Allrun-study [NPROC]      # dz = 5, 2.5, 1.25 mm; about 20 min on 3 cores
./Allrun-study 1 --levels 4 # adds dz = 0.625 mm (about 1 h serial)
./Allclean                  # keeps results.md and plots/
```

The case needs the circulation libraries of `src/circulationModels` with the
compressible support of `circulation1DCoupling`.

---

## Expected Results

`./Allrun` prints the largest errors against the exact solution, relative to
the incident peak, and plots the pressures at both ends of the pipe against
it (`matchedSegment3D.png`):

```text
max error / incident peak: P inlet 4.16e-03, P outlet 1.10e-02,
reflection 3.89e-03
```

At 5 mm cells the transmitted pulse is within 1.1 % of the exact one, its
peak arrives 0.5 % of the 32.4 ms delay late (the second-order dispersion of
the 3D scheme), and 0.39 % of the pulse is reflected back into the proximal
vessel, where the exact answer has none.

---

## Verification and Convergence Study

`./Allrun-study` halves the 3D cells along the pipe and the 3D step
together (acoustic Courant number 0.49), the vessels staying at fourth
order. Errors relative to the incident peak, the delay relative to 32.4 ms,
the observed order in brackets (`results.md`, `plots/matchedSegment3D.png`):

| dz [mm] | P, inlet | reflection | P, outlet | delay |
| --- | --- | --- | --- | --- |
| 5 | 4.16e-3 | 3.89e-3 | 1.10e-2 | 5.3e-3 |
| 2.5 | 1.03e-3 (2.02) | 9.59e-4 (2.02) | 2.66e-3 (2.05) | 1.3e-3 (2.07) |
| 1.25 | 2.56e-4 (2.00) | 2.38e-4 (2.01) | 6.51e-4 (2.03) | 3.1e-4 (2.04) |

- Every measure converges at second order, the order of the 3D solver
  (BDF2, `Gauss linear`): this is the 3D solver carrying the wave along the
  pipe, which the other cases never test. The coupling stays transparent
  up to that error: at 1.25 mm 0.024 % of the pulse is reflected, in the
  shape of the pulse's second time derivative (the discrete 3D line's
  impedance differs from `Zc` at second order in `k dz`).
- The 1D side is not involved: doubling its cells changes no figure to
  three digits.
- Explicit coupling of order 1 is stable at Courant 0.49 and reflects
  1.7 % of the pulse at 5 mm. Order 2 diverges at 0.49 and is stable at
  0.25, with the errors of the implicit coupling. A simple model in which
  the 3D side answers within a step with the gain `c0 dt/dz` (the measured
  implicit contraction) gives exactly that: order 1 stable below 1, order 2
  below 1/3.
- The Aitken-relaxed interface does not work here yet: it measures its
  relaxation from round-off while the pipe is at rest and freezes the
  interface. Use the Dirichlet interface with this case.

---

## References

- Boileau E, Nithiarasu P, Blanco PJ, et al. (2015) A benchmark study of
  numerical schemes for one-dimensional arterial blood flow modelling.
  Int J Numer Meth Biomed Engng 31:e02732.
