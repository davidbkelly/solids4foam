# Higher-order 1D scheme: `singlePulseHighOrder`

Prepared by David Kelly

---

## Tutorial Aims

- Verify the order of the k-exact reconstruction (`reconstruction kExact`,
  `polynomialOrder` 1-3) with a time scheme of the same order (`SSPRK2`,
  `SSPRK3`, `RK4`): orders 2, 3 and 4 in space and time.

---

## Case Overview

The inviscid single pulse of `singlePulseInviscid` (Boileau et al. 2015,
section 4.1) at full amplitude, refined at CFL 0.5 from 250 to 8000 cells, each
run against the 32000-cell run of its own set (self-convergence): the finest
profile is averaged onto the coarser cells, which is exact for cell averages.
A comparison with the linear solution would not do at third and fourth order:
the pulse steepens (2e-5 of the peak at an amplitude of 1e-3), and at smaller
amplitudes round-off takes over.

| set | reconstruction | time | expected |
| --- | --- | --- | --- |
| `MUSCL_SSPRK2` | MUSCL, limiter none | SSPRK2 | 2 |
| `kExact1_SSPRK2` | kExact, degree 1 | SSPRK2 | 2 |
| `kExact2_SSPRK3` | kExact, degree 2 | SSPRK3 | 3 |
| `kExact3_RK4` | kExact, degree 3 | RK4 | 4 |
| `kExact3_RK4_visc` | as above, with friction (mu = 4e-3) | RK4 | 4 |

```bash
./Allrun            # NPROC=4 runs at a time; orders.md, plots/highOrder.png
```

The runs are started as background jobs of the script's shell, not through
`xargs`: on macOS `xargs` and `sh` drop `DYLD_LIBRARY_PATH` and the runs would
not find the OpenFOAM libraries.

## Results

Order between the two finest meshes (Q; L1 / L2 / max): 2.05 / 2.04 / 2.04
(MUSCL and kExact 1), 3.00 / 2.99 / 2.98 (kExact 2), 4.01 / 4.01 / 4.00
(kExact 3), 4.01 / 4.01 / 4.01 (kExact 3 with friction). See `orders.md`.
