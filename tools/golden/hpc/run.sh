#!/bin/bash
#SBATCH -p serial
#SBATCH -J golden-run
# One arm of the golden battery on Wayland, or one element of it. The RiskMap installed in GOLDEN_HOME/library/COMMIT
# runs the battery of GOLDEN_HOME/src/BATTERY/tools/golden (BATTERY defaults to COMMIT, and the arms of a comparison
# should use the same battery), and the outputs go to GOLDEN_HOME/arms/ARM. TASKS is all (the default), smoke or a
# comma-separated list. With GROUP_FILE (groups.txt of the battery), element i of an array runs the tasks on line i
# instead. Tasks whose output exists are skipped, so an element can be resubmitted. Submitted by submit.sh,
# which adds the node type, memory, time and log files.
source /etc/profile >/dev/null 2>&1
set -euo pipefail
module add apptainer/1.5.3
GOLDEN_HOME=$(readlink -f "${GOLDEN_HOME:-$HOME/storage/riskmap-golden}")
CONTAINER=$(readlink -f "${CONTAINER:-$HOME/storage/containers/geospatial_4.6.1.sif}")
SHARED=$(readlink -f "${SHARED_LIBRARY:-$HOME/storage/riskmap-sim-env/library}")
BATTERY_DIR="${GOLDEN_HOME}/src/${BATTERY:-$COMMIT}/tools/golden"
LIBRARY="${GOLDEN_HOME}/library/${COMMIT}"
LIBRARIES="${LIBRARY}:${GOLDEN_HOME}/library/deps:${SHARED}"
TASKS=${TASKS:-all}
if [ -n "${GROUP_FILE:-}" ]; then
  TASKS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BATTERY_DIR}/${GROUP_FILE}")
fi
echo "node $(hostname), element ${SLURM_ARRAY_TASK_ID:-none}, RiskMap ${COMMIT}, battery ${BATTERY:-$COMMIT}, tasks ${TASKS}"
apptainer exec --bind /mmfs1 --pwd "${GOLDEN_HOME}" \
  --env R_LIBS_USER="${LIBRARIES}",OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,R_GC_MEM_GROW=3 \
  "${CONTAINER}" Rscript --vanilla "${BATTERY_DIR}/run.R" --lib="${LIBRARY}" --commit="${COMMIT}" \
  --out="${GOLDEN_HOME}/arms/${ARM}" --tasks="${TASKS}"
