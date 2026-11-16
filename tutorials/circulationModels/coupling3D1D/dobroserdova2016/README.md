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
python3 analyse.py --csv     # the exact solutions as CSV in reference/ (for other codes)
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
  exact solution: peak flow 1.596 against 1.596 (sine) and 2.263 against
  2.265 (step). The exact solution of the 3D run takes the pipe's inductance
  with the area of its 32-sided inlet patch (`A3D` = 0.9936 pi R^2, as the
  rigid-segment benchmark does), which the 3D solution has exactly; with
  pi R^2 instead, the 0.64 % error in L shifts the resonance and the RMS
  flow error over 21 of its periods reads 3.8e-2 instead of 7.7e-3. That
  7.7e-3 of the inflow (sine, dt = 1 ms) is the phase error of the BDF2 time
  scheme on the undamped 17.5 Hz resonance (57 steps per period), not a
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

## The step input: limiters and the unlimited schemes

`limiterStudy.py` (`limiters.md`, `plots/limiters.png`) runs the step case
1D-only (the pipe as an inertance node, c0 = 3.5 m/s, linearised) with the
MUSCL reconstruction and each limiter (minmod, vanLeer, MC, none) under
SSPRK2, and with the unlimited k-exact reconstructions of degree 2 (SSPRK3)
and 3 (RK4), at CFL 0.5 on 100 and 400 cells, against the exact
delay-equation solution, whose flow has a jump at every round trip
(amplified to 2.27 times the inflow) and whose total variation over 0.5 s
is 19.1 times the inflow:

| cells | scheme | overshoot / inflow | TV / TV exact | RMS error / inflow | max error / inflow |
|---|---|---|---|---|---|
| 100 | minmod | -0.003 | 0.859 | 1.55e-02 | 6.73e-02 |
| 100 | vanLeer | -0.001 | 0.929 | 1.08e-02 | 4.69e-02 |
| 100 | MC | -0.001 | 0.951 | 1.10e-02 | 4.98e-02 |
| 100 | none | -0.000 | 0.968 | 1.17e-02 | 5.08e-02 |
| 100 | kExact 2 | +0.000 | 0.969 | 4.00e-03 | 2.29e-02 |
| 100 | kExact 3 | +0.000 | 0.999 | 2.37e-03 | 2.00e-02 |
| 400 | minmod | -0.000 | 0.953 | 4.62e-03 | 2.53e-02 |
| 400 | vanLeer | -0.000 | 0.984 | 3.00e-03 | 2.31e-02 |
| 400 | MC | -0.000 | 0.992 | 3.13e-03 | 2.51e-02 |
| 400 | none | -0.000 | 1.003 | 3.48e-03 | 2.74e-02 |
| 400 | kExact 2 | +0.000 | 0.996 | 8.86e-04 | 1.05e-02 |
| 400 | kExact 3 | +0.000 | 1.002 | 4.69e-04 | 6.54e-03 |

No scheme overshoots the exact peak (the flow at the interface is set by
the inductor's integration of the pressure, which smooths what the 1D
scheme sends it). The limiters show as a deficit of total variation, the
smeared jumps: 14 % with minmod and 5 % with MC at 100 cells, 5 % and 1 %
at 400. The unlimited schemes keep the total variation within 0.3 % of
the exact one (a ringing of that size behind each jump: `none` 1.003 and
kExact 3 1.002 at 400 cells) and have a quarter to a sixth of the RMS
error; the jumps themselves are resolved at first order in dx whatever
the scheme (max error 2 to 7 % of the inflow at 100 cells, 0.7 to 2.7 % at
400). On this problem the unlimited high-order reconstruction is the
better choice; a limiter costs accuracy for no overshoot it would prevent.

## Order in time

`../orderStudy.py --case dobroserdova2016` repeats the 3D run of the hard
coupling (sine inflow, c0 = 3.5 m/s) from rest over the first period at dt
= 2, 1, 0.5, 0.25 and 0.125 ms, with the fluid ddt scheme and the implicit
coupling order both at q = 1 (`Euler`), 2 (`backward`), 3 (`BDF 3`) and 4
(`BDF 4`), the 1D vessel at CFL 0.1 and the PIMPLE loop to a pressure
residual of 1e-10. Observed orders by self-convergence (`orders.md`,
`plots/orders.png`), the RMS over the period of the port flow and pressure:

| dt triple [ms] | Euler | backward | BDF 3 | BDF 4 |
|---|---|---|---|---|
| 2 / 1 / 0.5 | 0.41 / 0.41 | 1.78 / 1.71 | 2.87 / 2.45 | 3.33 / 2.82 |
| 1 / 0.5 / 0.25 | 0.87 / 0.85 | 1.99 / 1.92 | 2.91 / 2.73 | 3.47 / 3.43 |
| 0.5 / 0.25 / 0.125 | 0.94 / 0.92 | 2.00 / 1.97 | 2.97 / 2.94 | 3.82 / 3.92 |

(Q / P.) Every family reaches its nominal order; the coarse triples are
pre-asymptotic because the first period from rest is the vessel's resonant
transient (tau = 2.5 ms against dt = 2 ms). The finest Euler, backward and
BDF 3 runs differ from the finest BDF 4 run by 1.2e-2, 3.4e-5 and 4.4e-7 of
the flow.

**Against the exact solution** (`../orderStudy.py --case dobroserdova2016
--exact`; `orders_linear.md`, `plots/orders_linear.png`). The same runs
with the vessel on the linearised 1D equations (`linearised true`, of which
the delay equations are the exact solution at any amplitude) at fourth order
in space and time (kExact 3, RK4, 200 cells), every run compared with the
exact solution (RMS error / RMS of the signal over the period; the delay
equations integrated with the pipe's patch area and 64000 steps per round
trip, accurate to 1e-8):

| dt [ms] | Euler Q / P | backward Q / P | BDF 3 Q / P | BDF 4 Q / P |
|---|---|---|---|---|
| 2 | 1.36e-1 / 2.72e-1 | 1.06e-2 / 3.30e-2 | 1.48e-3 / 1.37e-2 | 1.28e-3 / 1.12e-2 |
| 1 | 8.55e-2 / 1.58e-1 | 2.42e-3 / 9.17e-3 | 2.19e-4 / 2.61e-3 | 9.63e-5 / 1.89e-3 |
| 0.5 | 4.64e-2 / 8.67e-2 | 5.51e-4 / 2.36e-3 | 2.91e-5 / 6.71e-4 | 1.22e-5 / 5.90e-4 |
| 0.25 | 2.38e-2 / 4.49e-2 | 1.38e-4 / 7.79e-4 | 4.48e-6 / 5.20e-4 | 2.15e-6 / 5.30e-4 |
| 0.125 | 1.21e-2 / 2.29e-2 | 3.44e-5 / 5.37e-4 | 2.21e-6 / 5.25e-4 | 2.26e-6 / 5.29e-4 |
| orders | 0.67, 0.88, 0.96, 0.98 | 2.13, 2.14, 2.00, 2.00 | 2.76, 2.91, 2.70, 1.02 | 3.74, 2.98, 2.50, -0.07 |

(orders of Q.) Euler and backward converge at their orders to the finest
step; BDF 3 and BDF 4 at 2.9 and 3.7 until a floor of 2.2e-6 in Q, and
every family meets a floor of 5.3e-4 in P. Both floors are the 1D vessel's,
not the coupling's: the 1D-only run with the pipe as an `inertance` node,
same linearised equations, same scheme and cells, has 2.33e-6 and 5.29e-4
against the exact solution. The exact solution is not smooth: the sine
inflow starts from rest with a finite slope, so `P = L dQ/dt` at the pipe
jumps at every arrival of the wave front (every round trip, 28 ms), and Q
has a kink there; the 1D scheme resolves the jumps at first order in dx
(P: 1.19e-3, 5.3e-4, 2.4e-4, 1.1e-4 over 100 to 800 cells with kExact 3)
and the kinks at second (Q: 9.2e-6, 2.3e-6, 5.9e-7, 1.5e-7). The time
orders are therefore as the self-convergence says, and the comparison with
the exact solution is limited by the discontinuities of the problem posed,
which no time scheme sees.
