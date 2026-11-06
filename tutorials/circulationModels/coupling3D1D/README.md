# 3D-1D coupling benchmarks: `coupling3D1D`

Prepared by David Kelly

The test cases of the coupling of 3D CFD domains to 1D arterial networks
(`circulation1DCoupling`): the rigid-segment benchmark with its exact
solution, and the tests of the literature. The theory and the closed forms are
in `src/circulationModels/docs/BENCHMARKS3D1D.md`; the 1D solver's own
benchmarks are in [../circulation1D](../circulation1D).

| Case | What it is | Reference | What it demonstrates |
|---|---|---|---|
| `rigidSegmentInertance` | a rigid pipe (the `inertance` node) between two 1D vessels, 1D only | closed form | the exact solution of the benchmark; second and fourth order in dx |
| `rigidSegment3D` | the same pipe in pimpleFoam, coupled at both ends | closed form | the coupling's error, order in dt for BDF2/3/4 and coupling orders 1-4, stability limits, fourth order in space and time |
| `dobroserdova2016` | the hard- and soft-coupling tests of Dobroserdova, Olshanskii & Simakov (2016): a 1D vessel feeding a rigid pipe with an open end, sinusoidal and step inflows, with and without their absorbing 0D model at the interface | exact delay-equation solution of the lossless line with its lumped load | that the "spurious reflection" of the hard coupling is the physics of the problem posed (|Gamma| = 1, resonance of the vessel), and what the soft coupling's buffer does to it |
| `formaggia2007` | the pressure-pulse test of Formaggia, Moura & Nobile (2007), section 5 (Formaggia et al. 2001): the compliant FSI tube of the `3dTube` tutorial continued by a 1D vessel on its outlet | the same FSI tube three times as long | the coupling on a solids4foam FSI fluid (`circulation1DPressure` on the outlet of a Robin-Neumann FSI case), the 1D wall matched to the 3D wall's static compliance, the reflection at the interface |

```bash
./Allrun            # every case (the rigidSegment3D sweeps are separate: see its README)
```

Needs OpenFOAM.com, the `circulationModels` module built (`libcirculation1D`,
`libcirculation1DCoupling`, and `libsolids4FoamModels` for the `BDF` time
schemes), and a Python with numpy, scipy and matplotlib. `formaggia2007` is
an FSI case and needs solids4foam built with PETSc (it is skipped otherwise).
