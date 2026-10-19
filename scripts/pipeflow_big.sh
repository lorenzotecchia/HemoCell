#!/bin/bash
# Submit a bigger/longer pipeflow run on Snellius. Run from the repo root:
#   S=2 HT=0.3 NODES=2 TMAX=2000000 ./scripts/pipeflow_big.sh
# S scales the 50x25x25 um tube (S=2 -> 100x50x50 um, ~8x the cells).
set -e
S=${S:-2}; HT=${HT:-0.3}; NODES=${NODES:-2}; TMAX=${TMAX:-2000000}
ACCOUNT=${ACCOUNT:-ltecchia}; TIME=${TIME:-48:00:00}
ROOT=$(pwd)
D=${D:-/projects/prjs2303/$USER/pipeflow_big_S$S}
mkdir -p "$D"
cp examples/pipeflow/{pipeflow,RBC.xml,PLT.xml,config.xml} "$D"

# scale the ASCII STL
awk -v s=$S '/vertex/{printf "vertex %g %g %g\n",$2*s,$3*s,$4*s; next}1' examples/pipeflow/tube.stl > "$D/tube.stl"

# RBCs packed in the full box (pi/4 corrects for cells dropped at the wall, as in ht_sweep); PLT is a guess
RBC=$(python3 -c "import math;S=$S;print(round($HT*math.pi*(12.5*S)**2*50*S/97/(math.pi/4)))")
(cd "$D" && "$ROOT/tools/packCells/packCells" --rbc $RBC --plt $((RBC/20)) $((50*S)) $((25*S)) $((25*S)))

sed -i -e "s|<refDirN>.*</refDirN>|<refDirN> $((50*S)) </refDirN>|" \
       -e "s|<tmax>.*</tmax>|<tmax> $TMAX </tmax>|" \
       -e "s|<tcheckpoint>.*</tcheckpoint>|<tcheckpoint> 100000 </tcheckpoint>|" "$D/config.xml"

cd "$D"
sbatch -J pipe_big -A $ACCOUNT -p rome -N $NODES --ntasks-per-node 128 -t $TIME \
  --wrap "source $ROOT/scripts/snellius_env.sh; srun ./pipeflow config.xml"
