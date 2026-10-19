#!/bin/bash
# Run hdf5 -> xmf conversion as a Slurm job instead of on the login node.
# Usage (repo root): ./scripts/postprocess_job.sh <run dir containing output/>
#   ACCOUNT=ltecchia TIME=04:00:00 ./scripts/postprocess_job.sh $D
set -e
RUN=$(cd "${1:?run dir}" && pwd)
ROOT=$(pwd)
sbatch -J pipe_post -A ${ACCOUNT:-ltecchia} -p rome -n 1 -t ${TIME:-04:00:00} -o "$RUN/post-%j.out" --wrap "
module load 2025 SciPy-bundle h5py 2>/dev/null
python3 -c 'import h5py, numpy' 2>/dev/null || pip install --user h5py numpy
cd $RUN/output* && $ROOT/scripts/batchPostProcess.sh"
