#!/usr/bin/env python3
"""Errors against the in vitro data (Boileau et al. 2015, Table X) and Figs. 6-9.

    python3 plot.py

Writes metrics.md and the figures in plots/. Needs numpy and matplotlib.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "common"))

import boileauMetrics  # noqa: E402

boileauMetrics.analyse_network(HERE)
