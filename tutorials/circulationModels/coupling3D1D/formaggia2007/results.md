# Formaggia et al. (2007), section 5: the pulse tube continued by a 1D vessel

3D wave speed of full3D (half-peak arrival at z = 1.25, 2.5, 3.75, 5, 6.25, 7.5 cm): 4.917 m/s

3D wave speed of full3D_lib2 (half-peak arrival at z = 1.25, 2.5, 3.75, 5, 6.25, 7.5 cm): 4.912 m/s

1D wave speed, same measure (1D alone): 5.415 m/s

## coupled: 1D c0 = 5.25 m/s (against full3D)

| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] | arrival full3D / coupled / 1D alone [ms] | RMS(P coupled - full3D) / pulse | RMS(Q coupled - full3D) / peak Q |
|---|---|---|---|---|---|
| 1.25 | 3D | 1388 / 1388 / 1338 | 2.74 / 2.74 / 2.30 | 3.02e-02 | 5.35e-02 |
| 2.5 | 3D | 1288 / 1322 / 1332 | 5.26 / 5.19 / 4.60 | 2.55e-02 | 9.30e-02 |
| 3.75 | 3D | 1133 / 1174 / 1330 | 7.90 / 7.98 / 6.91 | 4.08e-02 | 8.55e-02 |
| 5 | port | 1093 / 1165 / 1328 | 10.42 / 10.20 / 9.22 | 4.16e-02 | 9.91e-02 |
| 6.25 | 1D | 1053 / 1156 / 1327 | 12.96 / 12.52 / 11.53 | 8.96e-02 | 1.51e-01 |
| 7.5 | 1D | 993 / 1152 / 1325 | 15.41 / 14.83 / 13.84 | 9.66e-02 | 2.35e-01 |
| 10 | 1D | 486 / 879 / 1322 | nan / 19.46 / 18.46 | 5.56e-02 | 2.43e-01 |
| 12.5 | 1D | 38 / 45 / 0 | nan / nan / nan | 2.35e-02 | 4.58e-01 |

Radial wall displacement at z = 2.5 cm: peak 150.32 um (full3D), 150.18 um (coupled); RMS difference 3.500 um

## coupled_c4p92: 1D c0 = 4.92 m/s (against full3D)

| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] | arrival full3D / coupled / 1D alone [ms] | RMS(P coupled - full3D) / pulse | RMS(Q coupled - full3D) / peak Q |
|---|---|---|---|---|---|
| 1.25 | 3D | 1388 / 1388 / 1338 | 2.74 / 2.74 / 2.30 | 2.75e-02 | 5.58e-02 |
| 2.5 | 3D | 1288 / 1321 / 1332 | 5.26 / 5.19 / 4.60 | 2.75e-02 | 9.45e-02 |
| 3.75 | 3D | 1133 / 1171 / 1330 | 7.90 / 7.98 / 6.91 | 4.27e-02 | 8.74e-02 |
| 5 | port | 1093 / 1124 / 1328 | 10.42 / 10.28 / 9.22 | 3.81e-02 | 1.07e-01 |
| 6.25 | 1D | 1053 / 1115 / 1327 | 12.96 / 12.74 / 11.53 | 6.41e-02 | 1.56e-01 |
| 7.5 | 1D | 993 / 1110 / 1325 | 15.41 / 15.21 / 13.84 | 5.44e-02 | 2.03e-01 |
| 10 | 1D | 486 / 620 / 1322 | nan / nan / 18.46 | 2.09e-02 | 1.97e-01 |
| 12.5 | 1D | 38 / 42 / 0 | nan / nan / nan | 2.13e-02 | 4.58e-01 |

Radial wall displacement at z = 2.5 cm: peak 150.32 um (full3D), 150.16 um (coupled_c4p92); RMS difference 3.315 um

## coupled_c4p92_closed_20ms: 1D c0 = 4.92 m/s (against full3D)

| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] | arrival full3D / coupled / 1D alone [ms] | RMS(P coupled - full3D) / pulse | RMS(Q coupled - full3D) / peak Q |
|---|---|---|---|---|---|
| 1.25 | 3D | 1388 / 1352 / 1338 | 2.74 / 2.71 / 2.30 | 3.13e-02 | 4.92e-02 |
| 2.5 | 3D | 1288 / 1303 / 1332 | 5.26 / 5.19 / 4.60 | 3.15e-02 | 9.24e-02 |
| 3.75 | 3D | 1133 / 1169 / 1330 | 7.90 / 7.93 / 6.91 | 4.81e-02 | 8.25e-02 |
| 5 | port | 1093 / 1104 / 1328 | 10.42 / 10.27 / 9.22 | 5.10e-02 | 1.19e-01 |
| 6.25 | 1D | 1053 / 1096 / 1327 | 12.96 / 12.72 / 11.53 | 7.55e-02 | 1.54e-01 |
| 7.5 | 1D | 993 / 1091 / 1325 | 15.41 / 15.18 / 13.84 | 6.48e-02 | 2.10e-01 |
| 10 | 1D | 486 / 632 / 1322 | nan / nan / 18.46 | 2.09e-02 | 1.96e-01 |
| 12.5 | 1D | 38 / 38 / 0 | nan / nan / nan | 2.12e-02 | 4.58e-01 |

Radial wall displacement at z = 2.5 cm: peak 150.32 um (full3D), 146.76 um (coupled_c4p92_closed_20ms); RMS difference 4.242 um

## coupled_c4p92_free: 1D c0 = 4.92 m/s (against full3D)

| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] | arrival full3D / coupled / 1D alone [ms] | RMS(P coupled - full3D) / pulse | RMS(Q coupled - full3D) / peak Q |
|---|---|---|---|---|---|
| 1.25 | 3D | 1388 / 1388 / 1338 | 2.74 / 2.74 / 2.30 | 4.69e-02 | 6.59e-02 |
| 2.5 | 3D | 1288 / 1249 / 1332 | 5.26 / 5.34 / 4.60 | 4.79e-02 | 1.36e-01 |
| 3.75 | 3D | 1133 / 1107 / 1330 | 7.90 / 7.84 / 6.91 | 2.66e-02 | 2.15e-01 |
| 5 | port | 1093 / 1110 / 1328 | 10.42 / 10.47 / 9.22 | 2.03e-02 | 1.15e-01 |
| 6.25 | 1D | 1053 / 1102 / 1327 | 12.96 / 12.93 / 11.53 | 3.70e-02 | 1.54e-01 |
| 7.5 | 1D | 993 / 1098 / 1325 | 15.41 / 15.39 / 13.84 | 3.43e-02 | 1.77e-01 |
| 10 | 1D | 486 / 549 / 1322 | nan / nan / 18.46 | 2.11e-02 | 1.83e-01 |
| 12.5 | 1D | 38 / 11 / 0 | nan / nan / nan | 2.12e-02 | 4.56e-01 |

Radial wall displacement at z = 2.5 cm: peak 150.32 um (full3D), 150.14 um (coupled_c4p92_free); RMS difference 3.525 um

## coupled_c4p92_free_lib2: 1D c0 = 4.92 m/s (against full3D_lib2)

| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] | arrival full3D / coupled / 1D alone [ms] | RMS(P coupled - full3D) / pulse | RMS(Q coupled - full3D) / peak Q |
|---|---|---|---|---|---|
| 1.25 | 3D | 1352 / 1352 / 1338 | 2.71 / 2.71 / 2.30 | 4.29e-02 | 6.01e-02 |
| 2.5 | 3D | 1271 / 1233 / 1332 | 5.26 / 5.33 / 4.60 | 4.39e-02 | 1.26e-01 |
| 3.75 | 3D | 1132 / 1108 / 1330 | 7.89 / 7.82 / 6.91 | 2.50e-02 | 1.98e-01 |
| 5 | port | 1077 / 1093 / 1328 | 10.42 / 10.47 / 9.22 | 1.84e-02 | 1.10e-01 |
| 6.25 | 1D | 1045 / 1085 / 1327 | 12.94 / 12.93 / 11.53 | 3.63e-02 | 1.44e-01 |
| 7.5 | 1D | 984 / 1081 / 1325 | 15.41 / 15.38 / 13.84 | 3.30e-02 | 1.73e-01 |
| 10 | 1D | 492 / 555 / 1322 | nan / nan / 18.46 | 2.09e-02 | 1.71e-01 |
| 12.5 | 1D | 36 / 11 / 0 | nan / nan / nan | 1.98e-02 | 4.55e-01 |

Radial wall displacement at z = 2.5 cm: peak 146.89 um (full3D_lib2), 146.74 um (coupled_c4p92_free_lib2); RMS difference 3.220 um

## coupled_c4p92_lib2: 1D c0 = 4.92 m/s (against full3D_lib2)

| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] | arrival full3D / coupled / 1D alone [ms] | RMS(P coupled - full3D) / pulse | RMS(Q coupled - full3D) / peak Q |
|---|---|---|---|---|---|
| 1.25 | 3D | 1352 / 1352 / 1338 | 2.71 / 2.71 / 2.30 | 2.59e-02 | 5.04e-02 |
| 2.5 | 3D | 1271 / 1303 / 1332 | 5.26 / 5.19 / 4.60 | 3.01e-02 | 9.04e-02 |
| 3.75 | 3D | 1132 / 1169 / 1330 | 7.89 / 7.93 / 6.91 | 4.59e-02 | 8.24e-02 |
| 5 | port | 1077 / 1104 / 1328 | 10.42 / 10.27 / 9.22 | 4.40e-02 | 1.10e-01 |
| 6.25 | 1D | 1045 / 1096 / 1327 | 12.94 / 12.72 / 11.53 | 6.87e-02 | 1.48e-01 |
| 7.5 | 1D | 984 / 1091 / 1325 | 15.41 / 15.18 / 13.84 | 5.99e-02 | 2.05e-01 |
| 10 | 1D | 492 / 632 / 1322 | nan / nan / 18.46 | 2.03e-02 | 1.84e-01 |
| 12.5 | 1D | 36 / 38 / 0 | nan / nan / nan | 1.99e-02 | 4.56e-01 |

Radial wall displacement at z = 2.5 cm: peak 146.89 um (full3D_lib2), 146.76 um (coupled_c4p92_lib2); RMS difference 3.599 um

## coupled_c4p92_spring: 1D c0 = 4.92 m/s (against full3D_lib2)

| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] | arrival full3D / coupled / 1D alone [ms] | RMS(P coupled - full3D) / pulse | RMS(Q coupled - full3D) / peak Q |
|---|---|---|---|---|---|
| 1.25 | 3D | 1352 / 1352 / 1338 | 2.71 / 2.71 / 2.30 | 4.23e-02 | 5.89e-02 |
| 2.5 | 3D | 1271 / 1233 / 1332 | 5.26 / 5.33 / 4.60 | 5.11e-02 | 1.31e-01 |
| 3.75 | 3D | 1132 / 1111 / 1330 | 7.89 / 7.82 / 6.91 | 2.52e-02 | 2.13e-01 |
| 5 | port | 1077 / 1098 / 1328 | 10.42 / 10.43 / 9.22 | 2.21e-02 | 1.13e-01 |
| 6.25 | 1D | 1045 / 1089 / 1327 | 12.94 / 12.89 / 11.53 | 4.62e-02 | 1.43e-01 |
| 7.5 | 1D | 984 / 1084 / 1325 | 15.41 / 15.35 / 13.84 | 3.76e-02 | 1.78e-01 |
| 10 | 1D | 492 / 572 / 1322 | nan / nan / 18.46 | 2.04e-02 | 1.71e-01 |
| 12.5 | 1D | 36 / 11 / 0 | nan / nan / nan | 1.98e-02 | 4.55e-01 |

Radial wall displacement at z = 2.5 cm: peak 146.89 um (full3D_lib2), 146.74 um (coupled_c4p92_spring); RMS difference 4.068 um
