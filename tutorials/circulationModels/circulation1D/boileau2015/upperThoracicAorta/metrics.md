# Upper thoracic aorta: errors against the 3D solution [%]

Boileau et al. (2015) Table VI with a solids4foam column. solids4foam cycle: t = 18.1470 to 19.1020 s; cycle-to-cycle relative RMS change of P at the probe 6.5e-06 (periodic at 1e-4).

| Error | DCG | LCG | FEM | FVM | McC | STM | solids4foam | vs range |
|---|---|---|---|---|---|---|---|---|
| E^RMS_P | 1.09 | 1.21 | 1.09 | 1.14 | 1.08 | 1.12 | **1.20** | in range |
| E^RMS_Q | 2.55 | 2.19 | 2.23 | 2.17 | 2.22 | 2.28 | **2.20** | in range |
| E^RMS_dr | 2.32 | 2.41 | 2.33 | 2.44 | 2.33 | 2.90 | **2.40** | in range |
| E^RMS_dP | 7.70 | 7.23 | 7.23 | 7.14 | 7.23 | 7.52 | **7.24** | in range |
| E^MAX_P | 3.54 | 3.55 | 3.26 | 3.18 | 3.27 | 3.30 | **3.57** | above (+0.02) |
| E^MAX_Q | 9.13 | 7.20 | 7.04 | 7.07 | 7.04 | 7.40 | **7.22** | in range |
| E^MAX_dr | 7.80 | 7.67 | 7.32 | 7.20 | 7.09 | 8.04 | **7.73** | in range |
| E^MAX_dP | 31.83 | 29.31 | 29.15 | 29.13 | 29.31 | 31.19 | **29.34** | in range |
| E^SYS_P | 0.57 | 0.32 | 0.56 | 0.71 | 0.55 | 0.59 | **0.30** | below range |
| E^SYS_Q | 5.56 | 5.36 | 5.63 | 5.29 | 5.58 | 5.56 | **5.37** | in range |
| E^SYS_dr | 2.54 | 1.93 | 2.50 | 2.86 | 2.49 | 4.48 | **1.91** | below range |
| E^SYS_dP | 8.92 | 8.89 | 9.29 | 8.53 | 9.04 | 9.70 | **8.84** | in range |
| E^DIAS_P | 0.85 | 1.12 | 0.88 | 0.99 | 0.87 | 0.85 | **1.11** | in range |
| E^DIAS_Q | 2.75 | 2.67 | 3.09 | 3.48 | 2.67 | 12.95 | **2.69** | in range |
| E^DIAS_dr | 1.85 | 2.24 | 1.89 | 2.05 | 1.89 | 0.00 | **2.20** | in range |
| E^DIAS_dP | 7.32 | 7.24 | 7.60 | 6.37 | 7.32 | 7.47 | **7.21** | in range |
