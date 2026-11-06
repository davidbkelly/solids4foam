#!/usr/bin/env python3
"""Errors against the 3D solution (Boileau et al. 2015) and the midpoint waveforms.

    python3 plot.py

Writes metrics.md (the paper's table with a solids4foam column) and the figure in
plots/. Needs numpy and matplotlib.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "common"))

import boileauMetrics  # noqa: E402

boileauMetrics.analyse_case(HERE.name, HERE)
