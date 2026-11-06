# Aortic bifurcation: errors against the 3D solution [%]

Boileau et al. (2015) Table VIII with a solids4foam column. solids4foam cycle: t = 28.6010 to 29.7010 s; cycle-to-cycle relative RMS change of P at the aorta midpoint 3.6e-06.


## Midpoint of the aorta

| Error | DCG | LCG | FEM | FVM | McC | STM | solids4foam | vs range |
|---|---|---|---|---|---|---|---|---|
| E^RMS_P | 0.37 | 0.39 | 0.37 | 0.40 | 0.38 | 0.38 | **0.37** | in range |
| E^RMS_Q | 0.89 | 0.71 | 0.82 | 0.81 | 0.83 | 0.82 | **0.82** | in range |
| E^RMS_dr | 2.48 | 2.43 | 2.55 | 2.49 | 2.48 | 2.45 | **2.47** | in range |
| E^MAX_P | 0.61 | 0.67 | 0.65 | 0.68 | 0.66 | 0.66 | **0.66** | in range |
| E^MAX_Q | 2.51 | 2.24 | 2.53 | 2.44 | 2.56 | 2.52 | **2.53** | in range |
| E^MAX_dr | 4.05 | 3.90 | 4.03 | 4.01 | 3.98 | 3.94 | **3.95** | in range |
| E^SYS_P | 0.53 | 0.54 | 0.51 | 0.55 | 0.53 | 0.51 | **0.51** | in range |
| E^SYS_Q | 2.47 | 2.23 | 2.52 | 2.44 | 2.55 | 2.52 | **2.52** | in range |
| E^SYS_dr | 3.96 | 3.86 | 3.99 | 4.00 | 3.94 | 3.93 | **3.91** | in range |
| E^DIAS_P | 0.50 | 0.51 | 0.47 | 0.53 | 0.49 | 0.48 | **0.48** | in range |
| E^DIAS_Q | 1.16 | 1.03 | 1.16 | 1.16 | 1.17 | 1.17 | **1.16** | in range |
| E^DIAS_dr | 1.39 | 1.37 | 1.49 | 1.35 | 1.40 | 1.41 | **1.41** | in range |

## End of the aorta

| Error | DCG | LCG | FEM | FVM | McC | STM | solids4foam | vs range |
|---|---|---|---|---|---|---|---|---|
| E^RMS_P | 0.42 | 0.43 | 0.41 | 0.44 | 0.42 | 0.42 | **0.42** | in range |
| E^RMS_Q | 1.20 | 1.13 | 1.12 | 1.08 | 1.12 | 0.46 | **1.12** | in range |
| E^RMS_dr | 4.07 | 4.02 | 4.14 | 4.08 | 4.08 | 4.05 | **4.06** | in range |
| E^MAX_P | 0.73 | 0.75 | 0.72 | 0.80 | 0.74 | 0.73 | **0.73** | in range |
| E^MAX_Q | 3.55 | 3.54 | 3.53 | 3.28 | 3.53 | 1.48 | **3.53** | in range |
| E^MAX_dr | 6.87 | 6.76 | 6.89 | 6.95 | 6.85 | 6.98 | **6.82** | in range |
| E^SYS_P | 0.71 | 0.70 | 0.69 | 0.72 | 0.70 | 0.69 | **0.69** | in range |
| E^SYS_Q | 3.47 | 3.53 | 3.52 | 3.26 | 3.52 | 1.08 | **3.52** | in range |
| E^SYS_dr | 6.84 | 6.73 | 6.87 | 6.88 | 6.83 | 6.90 | **6.80** | in range |
| E^DIAS_P | 0.54 | 0.55 | 0.53 | 0.57 | 0.54 | 0.53 | **0.53** | in range |
| E^DIAS_Q | 1.76 | 1.78 | 1.76 | 1.74 | 1.77 | 0.69 | **1.76** | in range |
| E^DIAS_dr | 1.93 | 1.92 | 2.02 | 1.89 | 1.94 | 1.88 | **1.95** | in range |

## Midpoint of an iliac artery

| Error | DCG | LCG | FEM | FVM | McC | STM | solids4foam | vs range |
|---|---|---|---|---|---|---|---|---|
| E^RMS_P | 0.45 | 0.45 | 0.44 | 0.47 | 0.45 | 0.44 | **0.44** | in range |
| E^RMS_Q | 0.92 | 0.66 | 0.65 | 0.68 | 0.67 | 0.66 | **0.64** | below range |
| E^RMS_dr | 4.29 | 4.30 | 4.35 | 4.31 | 4.31 | 4.26 | **4.30** | in range |
| E^MAX_P | 0.84 | 0.86 | 0.84 | 0.92 | 0.86 | 0.85 | **0.85** | in range |
| E^MAX_Q | 2.33 | 2.00 | 2.01 | 1.80 | 2.01 | 2.02 | **2.01** | in range |
| E^MAX_dr | 7.29 | 7.36 | 7.37 | 7.48 | 7.35 | 7.34 | **7.33** | in range |
| E^SYS_P | 0.83 | 0.83 | 0.81 | 0.85 | 0.82 | 0.81 | **0.81** | in range |
| E^SYS_Q | 1.68 | 1.84 | 1.81 | 1.58 | 1.89 | 1.79 | **1.81** | in range |
| E^SYS_dr | 7.25 | 7.26 | 7.27 | 7.31 | 7.25 | 7.22 | **7.23** | in range |
| E^DIAS_P | 0.55 | 0.55 | 0.54 | 0.57 | 0.54 | 0.54 | **0.53** | below range |
| E^DIAS_Q | 1.17 | 1.19 | 1.17 | 1.18 | 1.19 | 1.18 | **1.17** | in range |
| E^DIAS_dr | 2.36 | 2.37 | 2.45 | 2.35 | 2.38 | 2.41 | **2.39** | in range |
