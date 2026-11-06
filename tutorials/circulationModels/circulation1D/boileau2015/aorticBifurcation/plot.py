#!/usr/bin/env python3
"""Errors against the 3D solution (Boileau et al. 2015, Table VIII) and Fig. 5.

    python3 plot.py

Writes metrics.md and the figure in plots/. Needs numpy and matplotlib.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "common"))

import boileauMetrics  # noqa: E402

boileauMetrics.analyse_bifurcation(HERE)
