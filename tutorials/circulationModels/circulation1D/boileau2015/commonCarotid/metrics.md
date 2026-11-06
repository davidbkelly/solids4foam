# Common carotid artery: errors against the 3D solution [%]

Boileau et al. (2015) Table IV with a solids4foam column. solids4foam cycle: t = 8.8560 to 9.9560 s; cycle-to-cycle relative RMS change of P at the probe 3.1e-08 (periodic at 1e-4).

| Error | DCG | LCG | FEM | FVM | McC | STM | solids4foam | vs range |
|---|---|---|---|---|---|---|---|---|
| E^RMS_P | 0.20 | 0.24 | 0.23 | 0.28 | 0.24 | 0.24 | **0.24** | in range |
| E^RMS_Q | 0.39 | 0.26 | 0.28 | 0.29 | 0.30 | 0.29 | **0.28** | in range |
| E^RMS_dr | 0.96 | 1.00 | 1.00 | 1.07 | 1.01 | 0.98 | **1.00** | in range |
| E^RMS_dP | 4.42 | 4.35 | 4.27 | 4.21 | 4.47 | 4.47 | **4.27** | in range |
| E^MAX_P | 0.34 | 0.52 | 0.50 | 0.66 | 0.52 | 0.53 | **0.51** | in range |
| E^MAX_Q | 1.23 | 1.02 | 1.12 | 1.07 | 1.20 | 1.13 | **1.12** | in range |
| E^MAX_dr | 1.75 | 1.93 | 1.88 | 2.28 | 1.93 | 1.96 | **1.91** | in range |
| E^MAX_dP | 16.69 | 16.17 | 15.88 | 15.60 | 16.58 | 16.90 | **15.89** | in range |
| E^SYS_P | 0.27 | 0.27 | 0.26 | 0.29 | 0.27 | 0.26 | **0.26** | in range |
| E^SYS_Q | 0.55 | 0.46 | 0.54 | 0.58 | 0.62 | 0.55 | **0.54** | in range |
| E^SYS_dr | 1.65 | 1.63 | 1.61 | 1.68 | 1.63 | 1.63 | **1.61** | in range |
| E^SYS_dP | 15.05 | 15.37 | 15.16 | 14.83 | 15.73 | 15.76 | **15.13** | in range |
| E^DIAS_P | 0.28 | 0.28 | 0.27 | 0.29 | 0.28 | 0.27 | **0.27** | in range |
| E^DIAS_Q | 0.25 | 0.22 | 0.26 | 0.24 | 0.23 | 0.27 | **0.26** | in range |
| E^DIAS_dr | 0.08 | 0.11 | 0.08 | 0.12 | 0.09 | 0.08 | **0.07** | below range |
| E^DIAS_dP | 4.99 | 4.76 | 5.03 | 4.91 | 4.45 | 4.85 | **5.00** | in range |
