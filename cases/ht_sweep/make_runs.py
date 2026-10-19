#!/usr/bin/env python3
"""Build one run dir per hematocrit (pos files + edited config.xml) and an sbatch array script.

Usage: ./make_runs.py [--base examples/pipeflow] [--out $SCRATCH/ht_sweep] [--ht 0.1 0.2 0.3 0.4 0.45]
Then:  sbatch <out>/submit.sh
"""
import argparse, math, os, shutil, subprocess, xml.etree.ElementTree as ET

ap = argparse.ArgumentParser()
ap.add_argument("--base", default="../../examples/pipeflow")      # dir with pipeflow, tube.stl, *.xml
ap.add_argument("--packcells", default="../../tools/packCells/packCells")
ap.add_argument("--out", default=f"/projects/prjs2303/{os.environ.get('USER','me')}/ht_sweep")
ap.add_argument("--ht", type=float, nargs="+", default=[0.05, 0.10, 0.20, 0.30, 0.40])
ap.add_argument("--nplt", type=int, default=40)          # same PLT count in every run -> only Ht changes
ap.add_argument("--tmax", type=int, default=1_000_000)   # dt=1e-7 s -> 0.1 s; ~44 passes of the 50 um box
ap.add_argument("--tcsv", type=int, default=5000)        # cell positions (margination data)
ap.add_argument("--tmeas", type=int, default=5000)       # hdf5 for ParaView: 200 frames over the run
ap.add_argument("--ntasks", type=int, default=16)
ap.add_argument("--time", default="24:00:00")
a = ap.parse_args()

# Geometry from tube.stl / config: 103x53x53 lattice @ 0.5 um -> 50 x 25 x 25 um box, axis along x
L, R = 50.0, 12.5
V = math.pi * R**2 * L
RBC_VOL = 97.0
BOX_TO_PIPE = math.pi / 4   # pipe/box volume: HemoCell drops cells touching the wall, so pack the full box

def pack(ht, d):
    n = round(ht * V / RBC_VOL / BOX_TO_PIPE)           # RBCs packed in the box
    subprocess.run([a.packcells, "--rbc", str(n), "--plt", str(round(a.nplt / BOX_TO_PIPE)),
                    str(int(L)), str(int(2 * R)), str(int(2 * R))], cwd=d, check=True, stdout=subprocess.DEVNULL)
    return n

os.makedirs(a.out, exist_ok=True)
names = []
for ht in a.ht:
    name = f"ht{round(ht*100):02d}"; names.append(name)
    d = f"{a.out}/{name}"; os.makedirs(d, exist_ok=True)
    for f in ("pipeflow", "tube.stl", "RBC.xml", "PLT.xml"):
        shutil.copy(f"{a.base}/{f}", d)
    n = pack(ht, d)
    t = ET.parse(f"{a.base}/config.xml"); r = t.getroot()
    for path, v in (("sim/tmax", a.tmax), ("sim/tcsv", a.tcsv), ("sim/tmeas", a.tmeas),
                    ("sim/tcheckpoint", 100000), ("parameters/outputDirectory", "output")):
        r.find(path).text = str(v)
    t.write(f"{d}/config.xml")
    print(f"{name}: nominal Ht={ht} -> {n} RBC packed (real Ht = RBCs kept * 97 / {V:.0f}, see log)")

open(f"{a.out}/submit.sh", "w").write(f"""#!/bin/bash
#SBATCH -J ht_sweep
#SBATCH --account prjs2303
#SBATCH --partition rome
#SBATCH --ntasks {a.ntasks}
#SBATCH --time {a.time}
#SBATCH --array 0-{len(names)-1}
#SBATCH --output {a.out}/slurm-%a.out
source {os.path.abspath('../../scripts/snellius_env.sh')}
RUNS=({' '.join(names)})
cd {a.out}/${{RUNS[$SLURM_ARRAY_TASK_ID]}}
# resubmitting the same array resumes from the last checkpoint (checkpoint.xml replaces config.xml as argv[1])
if [ -f output/checkpoint/checkpoint.xml ]; then srun ./pipeflow output/checkpoint/checkpoint.xml; else srun ./pipeflow config.xml; fi
""")
print(f"sbatch {a.out}/submit.sh")
