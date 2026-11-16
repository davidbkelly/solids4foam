# Test-circulation1DPorts

A check of the 3D-1D coupling ports of `circulation1D` (`coupledFlow` nodes)
before any CFD: a rigid pipe in plug flow, `L dQ/dt = Pin - Pout`, stands in
for the 3D domain and is integrated with the 3D time step by BDF1 to BDF4
(`Euler`, `backward`, `BDF3`, `BDF4`), as the CFD would be. Each step exchanges with the 1D network through the ports
exactly as the 3D-0D coupler does (`circulationCoupler`): the same forcing
polynomials, explicit or implicit (plain fixed point), orders 1 to 4.

The case (`case/`) is `tutorials/circulationModels/coupling3D1D/
rigidSegmentInertance` with the `inertance` node replaced by two ports, so the
monolithic run of that tutorial is the coupled one's limit.

```bash
cd case && Test-circulation1DPorts          # one run, settings in controlDict
python3 sweep.py WORKDIR                     # accuracy and stability, backward
python3 sweep.py WORKDIR --ddt BDF3,BDF4 --cfl 0.1
```

With `tau = l/(2 c0)` (6.48 ms here) the coupling number is `dt/tau`. The
sweep measures:

- **explicit stability**: the amplification per step at 0.9, 1.1 and 1.25
  times each predicted limit, against the model problem's largest root (the
  limits are tabulated in `sweep.py` and `BENCHMARKS3D1D.md`);
- **the implicit iteration's contraction rate**, against `dt/(a0 tau)`, `a0`
  the leading BDF coefficient;
- **accuracy**: the error of the pipe flow against the same coupling at a step
  eight times finer, so its order in `dt`, and the difference from the
  monolithic run.

`coupling { interface robin; }` in `system/controlDict` makes the pipe see
the ports through the coupler's Robin (impedance) interface, `P = P* -/+
Zc (Q - Q*)` solved implicitly with the pipe (`Q*` the flow the network was
integrated with, `Zc = rho c/A`): for this surrogate the ports' response
within a step is exactly that slope, so one exchange per step reproduces
the converged Dirichlet iteration (to 3e-9 of the peak at `dt/tau` =
0.077) and the explicit scheme is stable at every order to `dt/tau = 16`
(tested), where the Dirichlet limits are 4, 1.33, 0.57 and 0.27.

See `src/circulationModels/docs/BENCHMARKS3D1D.md`.
