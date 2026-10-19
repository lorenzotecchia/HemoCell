#!/usr/bin/env python3
"""Margination vs hematocrit from the HemoCell CSV output. Usage: ./analyze.py <sweep_dir>  (needs matplotlib)"""
import csv
import glob
import math
import os
import sys

import matplotlib.pyplot as plt

import matplotlib; matplotlib.use("Agg")

R, C, UM = 12.5, 12.5, 1e6          # tube radius / axis (y,z) in um; CSV positions are in metres
WALL = 2.5                           # "near wall" = PLT centre within this distance of the wall [um]
sweep = sys.argv[1]

def series(run, cell):
    for f in sorted(glob.glob(f"{run}/output/csv/{cell}.*.csv")):
        step = int(f.split(".")[-2])
        rows = list(csv.DictReader(open(f)))
        yield step, [R - math.hypot(float(r["Y"]) * UM - C, float(r["Z"]) * UM - C) for r in rows]  # dist to wall

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
out = open(f"{sweep}/margination.csv", "w"); out.write("run,step,mean_dist_wall_um,frac_near_wall\n")
for run in sorted(glob.glob(f"{sweep}/ht*")):
    name = os.path.basename(run); t, d, fr = [], [], []
    for step, dist in series(run, "PLT"):
        t.append(step * 1e-7 * 1e3)  # ms; dt = 1e-7 s
        d.append(sum(dist) / len(dist)); fr.append(sum(x < WALL for x in dist) / len(dist))
        out.write(f"{name},{step},{d[-1]:.3f},{fr[-1]:.3f}\n")
    ax[0].plot(t, d, label=name); ax[1].plot(t, fr, label=name)
ax[0].set(xlabel="time [ms]", ylabel="mean PLT distance to wall [µm]"); ax[0].axhline(R / 2, ls=":", c="k")  # dotted = uniform
ax[1].set(xlabel="time [ms]", ylabel=f"fraction of PLT within {WALL} µm of wall")
ax[1].axhline(1 - ((R - WALL) / R) ** 2, ls=":", c="k")   # uniform-distribution expectation
ax[0].legend(); plt.tight_layout(); plt.savefig(f"{sweep}/margination.png", dpi=150)
print("wrote margination.png / margination.csv")
