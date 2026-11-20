# Mesh refinement of the coupled tube, per end condition of the wall

Level 2: twice the cells in every direction of both regions, same step (25 us). RMS(P coupled - P full) over the window / pulse at the sections (z < 5 cm in the 3D part, 5 cm the port, beyond it the 1D vessel), each level against the full tube of the same level; D the radial wall displacement at z = 2.5 cm.

| end | level | run / reference | z = 1.25 | z = 2.5 | z = 3.75 | z = 5 | z = 6.25 | z = 7.5 | D peak full / coupled [um] | RMS dD / peak |
|---|---|---|---|---|---|---|---|---|---|---|
| clamped | 1 | coupled_c5p25 / full3D | 2.88e-02 | 2.80e-02 | 4.55e-02 | 4.75e-02 | 9.41e-02 | 1.03e-01 | 146.9 / 146.8 | 2.50e-02 |
| clamped | 2 | coupled_c5p25_x2 / full3D_x2 | 2.79e-02 | 2.74e-02 | 4.47e-02 | 4.57e-02 | 9.20e-02 | 1.01e-01 | 145.2 / 145.0 | 2.50e-02 |
| free | 1 | coupled_c5p25_free / full3D | 3.82e-02 | 4.06e-02 | 3.24e-02 | 2.20e-02 | 5.77e-02 | 6.55e-02 | 146.9 / 146.7 | 1.83e-02 |
| free | 2 | coupled_c5p25_free_x2 / full3D_x2 | 3.48e-02 | 3.70e-02 | 3.05e-02 | 2.08e-02 | 5.72e-02 | 6.66e-02 | 145.2 / 145.1 | 1.67e-02 |
| spring | 1 | coupled_c5p25_spring / full3D | 3.78e-02 | 4.85e-02 | 3.28e-02 | 2.52e-02 | 6.85e-02 | 7.46e-02 | 146.9 / 146.7 | 2.48e-02 |
| spring | 2 | coupled_c5p25_spring_x2 / full3D_x2 | 3.45e-02 | 4.44e-02 | 3.11e-02 | 2.39e-02 | 6.70e-02 | 7.38e-02 | 145.2 / 145.1 | 2.27e-02 |

Mesh sensitivity, level 2 against level 1 of the same run (RMS / pulse), and the differences between end conditions at level 2 (RMS / pulse):

| what | z = 1.25 | z = 2.5 | z = 3.75 | z = 5 | z = 6.25 | z = 7.5 |
|---|---|---|---|---|---|---|
| full tube, level 2 - level 1 | 3.13e-03 | 4.99e-03 | 6.31e-03 | 6.93e-03 | 8.14e-03 | 6.17e-03 |
| clamped, level 2 - level 1 | 3.55e-03 | 5.39e-03 | 6.93e-03 | 1.08e-02 | 1.01e-02 | 9.46e-03 |
| free, level 2 - level 1 | 4.24e-03 | 5.71e-03 | 6.50e-03 | 7.82e-03 | 7.39e-03 | 6.10e-03 |
| spring, level 2 - level 1 | 4.08e-03 | 5.89e-03 | 5.85e-03 | 1.05e-02 | 9.74e-03 | 8.77e-03 |
| clamped - free, level 2 | 6.16e-02 | 4.11e-02 | 5.87e-02 | 5.06e-02 | 4.26e-02 | 3.97e-02 |
| clamped - spring, level 2 | 6.12e-02 | 4.39e-02 | 5.20e-02 | 4.27e-02 | 3.32e-02 | 3.12e-02 |
| free - spring, level 2 | 2.41e-03 | 9.18e-03 | 1.02e-02 | 1.32e-02 | 1.26e-02 | 9.81e-03 |
