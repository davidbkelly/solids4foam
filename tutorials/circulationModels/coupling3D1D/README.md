# 3D-1D coupling benchmarks: `coupling3D1D`

Prepared by David Kelly

The test cases of the coupling of 3D CFD domains to 1D arterial networks
(`circulation1DCoupling`): the rigid-segment benchmark with its exact
solution, and the tests of the literature. The theory and the closed forms are
in `src/circulationModels/docs/BENCHMARKS3D1D.md`; the 1D solver's own
benchmarks are in [../circulation1D](../circulation1D).

| Case | What it is | Reference | What it demonstrates |
|---|---|---|---|
| `rigidSegmentInertance` | a rigid pipe (the `inertance` node) between two 1D vessels, 1D only | closed form | the exact solution of the benchmark; second order in dx, and fourth order to 1e-9 with the linearised equations at full amplitude |
| `rigidSegment3D` | the same pipe in pimpleFoam, coupled at both ends | closed form | the coupling's error, order in dt for BDF2/3/4 and coupling orders 1-4, stability limits, fourth order in space and time by self-convergence and, with the linearised vessels, against the closed form to 2e-8 |
| `dobroserdova2016` | the hard- and soft-coupling tests of Dobroserdova, Olshanskii & Simakov (2016): a 1D vessel feeding a rigid pipe with an open end, sinusoidal and step inflows, with and without their absorbing 0D model at the interface | exact delay-equation solution of the lossless line with its lumped load | that the "spurious reflection" of the hard coupling is the physics of the problem posed (|Gamma| = 1, resonance of the vessel), and what the soft coupling's buffer does to it |
| `formaggia2007` | the pressure-pulse test of Formaggia, Moura & Nobile (2007), section 5 (Formaggia et al. 2001): the compliant FSI tube of the `3dTube` tutorial continued by a 1D vessel on its outlet | the same FSI tube three times as long | the coupling on a solids4foam FSI fluid (`circulation1DPressure` on the outlet of a Robin-Neumann FSI case), the 1D wall matched to the 3D wall's static compliance, the reflection at the interface |
| `womersleySegment` | a rigid no-slip pipe between two 1D vessels with a single-frequency incident wave at Womersley numbers 3, 6 and 12 (an axisymmetric wedge in pimpleFoam) | the transmission line with the pipe as Womersley's impedance: reflection, transmission, port pressures and flow, and the velocity profile, in closed form | the coupling with a 3D side that develops a real viscous profile: exact to 1e-4 in amplitude and 0.004 degrees in phase, profiles to 1 % |
| `thoracicAorta37` | the thoracic aorta II of the 37-artery in vitro network (Boileau et al. 2015) as a rigid no-slip cone in pimpleFoam, the other 36 vessels in 1D, at physiological amplitude, restarted from the periodic 1D state | the original 1D network, its exact inviscid rigid limit (an `inertance` node) and the in vitro measurements | a 3D segment in a validated physiological network (the configuration of Blanco et al. 2013 and Malossi et al. 2013 in finite volumes): the rigid segment's effect on the whole network (4 % of the pulse) against its exact limit, the coupling's cost at physiological amplitude |

```bash
./Allrun            # every case (the rigidSegment3D sweeps are separate: see its README)
python3 orderStudy.py --case dobroserdova2016|womersleySegment|thoracicAorta37
python3 orderStudy.py --summary        # timeOrders.md, plots/timeOrders.png
```

## Order in time on every case

`orderStudy.py` runs each pimpleFoam case from rest over a short window at
four halving time steps for the scheme families Euler, `backward`, `BDF 3`
and `BDF 4` (the fluid ddt scheme and the implicit 3D-1D coupling order at
the family's order q), the 1D vessels sub-stepping at CFL 0.1 and the PIMPLE
loop to a pressure residual of 1e-10, and measures the observed order by
self-convergence of the port pressures and flows (`<case>/orders.md`,
`<case>/plots/orders.png`; the summary over the cases is `timeOrders.md`
and `plots/timeOrders.png`). The FSI case has its own study
(`formaggia2007/orderStudy.py`). See each case's README for its table; the
summary:

| case | window | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|---|
| `rigidSegment3D` (exact solution, `sweep.py`) | pulse | 1.00 | 2.00 | 2.98 | 3.78 (implicit); explicit orders 1-4: 1.02 / 2.00 / 2.97 / 3.89 |
| `dobroserdova2016` (sine inflow, rigid pipe) | first period | 0.94 | 2.00 | 2.97 | 3.82 |
| `womersleySegment` (alpha = 6, viscous wedge) | ramped wave, 0.6 s | 1.00 | 1.99 | 2.92 | 3.76 |
| `thoracicAorta37` (37-artery network, cone) | first systole, 0.3 s | 0.98 | 1.95 | 2.34, then a floor | floor at 3e-6 of the signal (see its README) |
| `formaggia2007` (FSI tube, IQN-ILS) | pulse, 12 ms | 1 | 2 | 2 (fluid BDF 3, wall BDF 2) | 2 (fluid BDF 4, wall BDF 2) |

(finest triple, port flow; the rigid-segment entries are its BDF 2/3/4
columns with the coupling at the same order.) The coupling takes the order
of the 3D time scheme up to 4 on a rigid 3D domain, viscous or inviscid,
alone or in a 37-vessel network; on the FSI case it takes the FSI solver's
order, 2.

Needs OpenFOAM.com, the `circulationModels` module built (`libcirculation1D`,
`libcirculation1DCoupling`, and `libsolids4FoamModels` for the `BDF` time
schemes), and a Python with numpy, scipy and matplotlib. `formaggia2007` is
an FSI case and needs solids4foam built with PETSc (it is skipped otherwise).
