#!/usr/bin/env python3
"""Pressure along the vessel and peak decay of the single pulse (Boileau et al. 2015).

    python3 plot.py

Writes peaks.md and the figure in plots/. Needs numpy and matplotlib.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "common"))

import singlePulse  # noqa: E402

singlePulse.analyse(HERE)
