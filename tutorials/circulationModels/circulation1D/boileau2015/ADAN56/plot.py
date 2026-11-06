#!/usr/bin/env python3
"""Waveforms of the ADAN56 network against the six schemes (Boileau et al. 2015,
Figs. 11-16) and their differences.

    python3 plot.py

Writes metrics.md and the figures in plots/. Needs numpy and matplotlib.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "common"))

import boileauMetrics  # noqa: E402

boileauMetrics.analyse_adan56(HERE)
