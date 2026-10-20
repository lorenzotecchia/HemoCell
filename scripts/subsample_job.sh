#!/bin/bash
# Slurm job: pick every Nth hdf5 snapshot of a run, convert it to .xmf, so a small subset can be rsynced to a laptop.
# Usage (repo root): EVERY=20 NOFLUID=1 ./scripts/subsample_job.sh <run dir containing output/>
# Result: <run dir>/subset/{hdf5,*.xmf}; copy with  rsync -avPL snellius:<run dir>/subset/ ./subset/
# The last snapshot is skipped (may be truncated if the run was cancelled).
set -e
RUN=$(cd "${1:?run dir}" && pwd)
ROOT=$(pwd)
EVERY=${EVERY:-20}; NOFLUID=${NOFLUID:-}
sbatch -J pipe_subset -A ${ACCOUNT:-ltecchia} -p rome -n 1 -t ${TIME:-04:00:00} -o "$RUN/subset-%j.out" <<EOF
#!/bin/bash
set -e
module load 2025 SciPy-bundle h5py 2>/dev/null
python3 -c 'import h5py, numpy' 2>/dev/null || pip install --user h5py numpy
OUT=\$(ls -d $RUN/output* | head -1)   # ponytail: first output dir; set it by hand if there are several
S=$RUN/subset; rm -rf \$S; mkdir -p \$S/hdf5
for snap in \$(ls \$OUT/hdf5 | head -n -1 | awk 'NR%$EVERY==1'); do
  mkdir \$S/hdf5/\$snap
  for f in \$OUT/hdf5/\$snap/*; do
    [ -n "$NOFLUID" ] && [[ \$(basename \$f) == Fluid* ]] && continue
    ln -s \$f \$S/hdf5/\$snap/
  done
done
echo "linked \$(ls \$S/hdf5 | wc -l) snapshots"
cd \$S && $ROOT/scripts/batchPostProcess.sh
ls \$S/*.xmf
EOF
