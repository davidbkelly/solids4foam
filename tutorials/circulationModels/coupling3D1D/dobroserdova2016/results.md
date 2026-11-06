# Dobroserdova et al. (2016), sections 4.2-4.3: hard and soft coupling against the exact linear solution

## sine inflow

| c0 [m/s] | config | peak Q / peak Qin | peak Q, exact | Q error, RMS / max | P error, RMS / max |
|---|---|---|---|---|---|
| 3.5 | reference (10 cm, absorbing) | 1.000 | - | - | - |
| 3.5 | hard (1D + lumped pipe) | 1.594 | 1.594 | 3.3e-04 / 1.1e-03 | 3.7e-03 / 1.4e-02 |
| 3.5 | soft (1D + lumped pipe) | 1.128 | 1.128 | 2.2e-04 / 1.9e-03 | 1.3e-03 / 3.4e-03 |
| 3.5 | hard (1D + pimpleFoam pipe) | 1.596 | 1.594 | 3.8e-02 / 9.0e-02 | 3.5e-02 / 1.1e-01 |
| 7.0 | reference (10 cm, absorbing) | 1.000 | - | - | - |
| 7.0 | hard (1D + lumped pipe) | 1.228 | 1.229 | 2.9e-04 / 1.2e-03 | 2.7e-03 / 9.8e-03 |
| 7.0 | soft (1D + lumped pipe) | 1.022 | 1.022 | 8.9e-05 / 7.8e-04 | 6.5e-04 / 1.7e-03 |
| 10.5 | reference (10 cm, absorbing) | 1.000 | - | - | - |
| 10.5 | hard (1D + lumped pipe) | 1.143 | 1.143 | 2.6e-04 / 1.1e-03 | 2.2e-03 / 7.6e-03 |
| 10.5 | soft (1D + lumped pipe) | 1.010 | 1.010 | 5.1e-05 / 5.2e-04 | 4.0e-04 / 1.1e-03 |

Hard coupling, sine, periodic state: |Gamma| = 1 at every c0; with the inlet's own reflection the interface flow is Qin/|cos(w td) + i (w L/Zc) sin(w td)| in amplitude: 1.042 at 3.5 m/s, 1.011 at 7.0 m/s, 1.005 at 10.5 m/s; the runs start from rest and the undamped resonance of the vessel (about c0/(4 L1D) Hz) persists on top of it.

## step inflow

| c0 [m/s] | config | peak Q / peak Qin | peak Q, exact | Q error, RMS / max | P error, RMS / max |
|---|---|---|---|---|---|
| 3.5 | reference (10 cm, absorbing) | 1.000 | - | - | - |
| 3.5 | hard (1D + lumped pipe) | 2.266 | 2.267 | 9.5e-03 / 4.4e-02 | 2.6e-01 / 2.3e+00 |
| 3.5 | soft (1D + lumped pipe) | 1.348 | 1.466 | 3.1e-02 / 3.5e-01 | 6.5e-02 / 4.1e-01 |
| 3.5 | hard (1D + pimpleFoam pipe) | 2.263 | 2.267 | 1.7e-02 / 6.1e-02 | 2.5e-01 / 2.1e+00 |
| 7.0 | reference (10 cm, absorbing) | 1.000 | - | - | - |
| 7.0 | hard (1D + lumped pipe) | 2.317 | 2.326 | 1.1e-02 / 3.9e-02 | 2.2e-01 / 2.7e+00 |
| 7.0 | soft (1D + lumped pipe) | 1.348 | 1.348 | 2.2e-02 / 3.0e-01 | 7.4e-02 / 4.0e-01 |
| 10.5 | reference (10 cm, absorbing) | 1.000 | - | - | - |
| 10.5 | hard (1D + lumped pipe) | 2.317 | 2.326 | 1.8e-02 / 6.0e-02 | 3.2e-01 / 2.2e+00 |
| 10.5 | soft (1D + lumped pipe) | 1.346 | 1.466 | 1.8e-02 / 1.8e-01 | 7.0e-02 / 3.2e-01 |

Hard coupling, sine, periodic state: |Gamma| = 1 at every c0; with the inlet's own reflection the interface flow is Qin/|cos(w td) + i (w L/Zc) sin(w td)| in amplitude: 1.042 at 3.5 m/s, 1.011 at 7.0 m/s, 1.005 at 10.5 m/s; the runs start from rest and the undamped resonance of the vessel (about c0/(4 L1D) Hz) persists on top of it.
