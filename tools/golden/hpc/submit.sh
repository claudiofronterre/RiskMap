#!/bin/bash
# Submits one arm of the golden battery on Wayland, from the login node: the install of RiskMap at COMMIT, unless its
# library is there already (REINSTALL=1 forces it), and the run after it, both on one node type. The battery is that of
# BATTERY, by default COMMIT, and both commits must have been copied there by push.sh.
#   ~/storage/riskmap-golden/src/<battery>/tools/golden/hpc/submit.sh <commit> <arm> [<battery> [<tasks>]]
# <tasks> is all (the default), smoke, a comma-separated list, or groups for an array over the lines of groups.txt.
# NODE_TYPE (64core256G), EXCLUDE (comp06-[00-24]), MEMORY (4G) and TIME (02:00:00) can be set in the environment.
source /etc/profile >/dev/null 2>&1
set -euo pipefail
if [ $# -lt 2 ]; then echo "usage: submit.sh <commit> <arm> [<battery> [<tasks>]]" >&2; exit 2; fi
export COMMIT=$1 ARM=$2 BATTERY=${3:-$1} TASKS=${4:-all}
export GOLDEN_HOME=$(readlink -f "${GOLDEN_HOME:-$HOME/storage/riskmap-golden}")
HPC="${GOLDEN_HOME}/src/${BATTERY}/tools/golden/hpc"
for sha in "${COMMIT}" "${BATTERY}"; do
  if [ ! -d "${GOLDEN_HOME}/src/${sha}" ]; then echo "push ${sha} first with tools/golden/hpc/push.sh" >&2; exit 1; fi
done
mkdir -p "${GOLDEN_HOME}/logs" "${GOLDEN_HOME}/arms"
placement=(-C "node_type=${NODE_TYPE:-64core256G}" --exclude="${EXCLUDE:-comp06-[00-24]}")
after=()
if [ ! -d "${GOLDEN_HOME}/library/${COMMIT}/RiskMap" ] || [ -n "${REINSTALL:-}" ]; then
  install=$(sbatch --parsable "${placement[@]}" --export=ALL -o "${GOLDEN_HOME}/logs/install-${COMMIT}.%j.out" \
    -e "${GOLDEN_HOME}/logs/install-${COMMIT}.%j.err" "${HPC}/install.sh")
  after=(--dependency="afterok:${install}")
  echo "install of ${COMMIT}: job ${install}"
fi
array=()
log="${GOLDEN_HOME}/logs/%x.%j"
if [ "${TASKS}" = groups ]; then
  export GROUP_FILE=groups.txt
  array=(--array="1-$(grep -c . "${GOLDEN_HOME}/src/${BATTERY}/tools/golden/groups.txt")")
  log="${GOLDEN_HOME}/logs/%x.%A_%a"
fi
run=$(sbatch --parsable "${placement[@]}" ${after[@]+"${after[@]}"} ${array[@]+"${array[@]}"} -J "golden-${ARM}" \
  --mem="${MEMORY:-4G}" --time="${TIME:-02:00:00}" --export=ALL -o "${log}.out" -e "${log}.err" "${HPC}/run.sh")
echo "arm ${ARM}: job ${run}, outputs in ${GOLDEN_HOME}/arms/${ARM}"
