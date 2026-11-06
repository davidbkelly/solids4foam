# Hard and soft 3D-1D coupling: `dobroserdova2016`

Prepared by David Kelly

---

## Tutorial Aims

- Reproduce the tests of Dobroserdova, Olshanskii and Simakov, *Multiscale
  coupling of compliant and rigid walls blood flow models*, Int. J. Numer.
  Meth. Fluids 82 (2016) 799-817, sections 4.2 (sinusoidal inflow) and 4.3
  (step inflow): a 1D vessel coupled to a rigid 3D pipe, with the "hard"
  coupling (continuity of pressure and flow) and their "soft" coupling (an
  absorbing 0D model at the interface).
- Compare both with the exact linear solution of the problem they pose, which
  the paper does not have (its reference is a 1D-only vessel with an absorbing
  end).

---

## Case Overview

Their set-up, in cgs: a 1D vessel of 5 cm and radius 1 cm joined to a rigid
3D cylinder of 5 cm and radius 1 cm whose far end is open (zero traction);
rho = 1 g/cm3; an inflow Q = 0.1 sin(2 pi t/0.3) cm3/s, or 0.1 cm3/s for
t <= 0.1 s and then 0; c0 = 350, 700 and 1050 cm/s. The soft coupling puts an
elastic sphere (a compliance C) between two resistances at the interface,
with R1D0D + R0D3D = rho c0/S (their Eq. 13) and C = l S/(rho c0^2), l = 0.1 cm
(Eq. 14). Their reference is a 10 cm vessel with a non-reflecting outlet.

Here, in SI and inviscid (their vessel has Poiseuille friction, with a damping
time of about 3 s; the slip-wall 3D pipe is inviscid anyway), `makeCases.py`
writes, for each inflow and c0:

| run | what |
| --- | --- |
| `reference` | their reference: a 10 cm vessel, absorbing outlet |
| `hard` | the 5 cm vessel ending in the rigid pipe as an `inertance` node open to P = 0 (exact: plug flow in a rigid pipe is an inductor L = rho l/S) |
| `soft` | the same with their 0D model, a `lumped0D` circuit `(port R1 C R2 L ground)`, R1 = R2 = Zc/2 |
| `hard3D` (c0 = 3.5 m/s) | the hard coupling with the pipe in pimpleFoam: slip walls, `circulation1DPressure` at the inlet, p = 0 at the outlet |

**The exact solution.** The vessel is a lossless transmission line of
impedance `Zc = rho c0/S` and transit time `td = L/c0` = 14 ms at 3.5 m/s,
far shorter than the pulses, so the waves the interface reflects return from
the flow-source inlet (which reflects them whole) many times per period: the
vessel is a resonator between the inlet and the pipe. With `f` and `g` the
forward and backward waves at the interface (`P = f + g`, `Q = (f - g)/Zc`),
the inlet gives `f(t) = Zc Qin(t - td) + g(t - 2 td)`, and the pipe
`L dQ/dt = P` (hard) or the RCRL circuit (soft). `analyse.py` integrates these
delay equations to 1e-9. The open-ended rigid pipe reflects everything,
`Gamma = (i w L - Zc)/(i w L + Zc)`, |Gamma| = 1 with a phase lag; the
reference has no reflection at all.

```bash
./Allrun            # 18 1D runs (seconds) + 2 pimpleFoam runs (minutes); results.md, plots/
```

---

## Results

See `results.md` and `plots/dobroserdova_sine.png`, `plots/dobroserdova_step.png`.

- **Sinusoidal inflow.** The 1D runs with the lumped pipe agree with the exact
  solution to 1e-3 of the inflow amplitude (hard and soft, all three c0). The
  hard coupling gives interface flows of 1.59, 1.23 and 1.14 times the inflow
  at c0 = 3.5, 7 and 10.5 m/s, the soft 1.13, 1.02 and 1.01; the reference
  1.00. In the periodic state the hard-coupled amplitude would be only 1.04,
  1.01 and 1.005 (`1/|cos(w td) + i (w L/Zc) sin(w td)|`): the rest is the
  undamped resonance of the vessel (about `c0/(4 L)` = 17.5 Hz at 3.5 m/s),
  started from rest, which their inviscid-in-effect set-up never damps. The
  soft coupling's resistance R = Zc damps that resonance, which is why their
  Figure 3 shows it "much closer to the reference". Neither is wrong: the
  hard coupling is the exact solution of the problem posed, an open-ended
  rigid pipe fed by a flow source through a short compliant vessel; the
  reference is a different problem.
- **Step inflow.** Peak flow 2.27 times the inflow with hard coupling (the
  exact value: 2.27), 1.35 with soft (exact 1.47: the 1D scheme smears the
  jump the soft coupling passes instantly through R1). The RMS errors are
  listed in `results.md`; the max norms are set by the smeared jump.
- **3D runs.** The pimpleFoam pipe with the hard coupling reproduces the
  exact solution: peak flow 1.596 against 1.594 (sine) and 2.263 against
  2.267 (step). Its RMS error in the flow, 3.8e-2 of the inflow for the sine
  at dt = 1 ms, is the phase error of the BDF2 time scheme on the undamped
  17.5 Hz resonance over 21 of its periods (57 steps per period), not a
  coupling error: the 1D run with the lumped pipe, whose CFL step is 50 times
  smaller, is at 3.3e-4. The 3D rigid pipe is the inductor, and the coupling
  adds nothing of its own.

- **Against their figures.** Their Figure 3 (sine, hard) peaks at about 1.55,
  1.10 and 1.05 times the inflow at c0 = 350, 700 and 1050 cm/s and their
  soft coupling at about 1.1; their Figure 5 (step) at about 2.1 (hard) and
  1.3 (soft). Ours are a few percent higher: rerunning the hard cases with
  their Poiseuille friction (`fluid { mu 4e-3; }`) gives 1.593, 1.177 and
  1.109 for the sine (the friction damps the resonance over the first
  period) and leaves the step's first peak at 2.27; the rest is their
  scheme's dissipation and the reading of their figures. The pattern -
  hard well above soft and reference, the excess falling with c0, the
  sustained oscillation after the step - is the same.

What this adds to the paper's test: a reference that is exact instead of a
different problem, and the statement of what the soft coupling's buffer does
(it damps the vessel's resonance by terminating it with its characteristic
impedance at high frequency) rather than that it "reduces spurious
reflections".
