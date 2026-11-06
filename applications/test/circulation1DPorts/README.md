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

See `src/circulationModels/docs/BENCHMARKS3D1D.md`.
