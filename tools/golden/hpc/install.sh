#!/bin/bash
#SBATCH -p serial
#SBATCH -J golden-install
#SBATCH --mem=4G
#SBATCH --time=01:00:00
# Installs RiskMap at COMMIT on Wayland, from the clean archive that push.sh unpacked in GOLDEN_HOME/src/COMMIT, into
# GOLDEN_HOME/library/COMMIT inside the container. Imports of RiskMap missing from the shared library and from
# GOLDEN_HOME/library/deps are installed into deps first, so that every arm loads the same dependencies. Submitted by
# submit.sh, which adds the node type and the log files, or by hand as
#   sbatch -C node_type=64core256G --exclude='comp06-[00-24]' -o <log> --export=ALL,COMMIT=<sha> install.sh
source /etc/profile >/dev/null 2>&1
set -euo pipefail
module add apptainer/1.5.3
GOLDEN_HOME=$(readlink -f "${GOLDEN_HOME:-$HOME/storage/riskmap-golden}")
CONTAINER=$(readlink -f "${CONTAINER:-$HOME/storage/containers/geospatial_4.6.1.sif}")
SHARED=$(readlink -f "${SHARED_LIBRARY:-$HOME/storage/riskmap-sim-env/library}")
SOURCE="${GOLDEN_HOME}/src/${COMMIT}"
LIBRARY="${GOLDEN_HOME}/library/${COMMIT}"
DEPS="${GOLDEN_HOME}/library/deps"
LIBRARIES="${LIBRARY}:${DEPS}:${SHARED}"
mkdir -p "${LIBRARY}" "${DEPS}"
echo "node $(hostname), installing ${SOURCE} into ${LIBRARY}"
in_container() {
  apptainer exec --bind /mmfs1 --pwd "${GOLDEN_HOME}" --env R_LIBS_USER="${LIBRARIES}" "${CONTAINER}" "$@"
}
in_container Rscript --vanilla -e "
  imports <- read.dcf('${SOURCE}/DESCRIPTION', 'Imports')[1, 1]
  imports <- trimws(sub('[(].*', '', strsplit(gsub('\\n', ' ', imports), ',')[[1]]))
  missing <- imports[!vapply(imports, requireNamespace, logical(1), quietly = TRUE)]
  cat('missing imports:', if (length(missing)) missing else 'none', '\n')
  # Downloads from the compute nodes can be slow, and R gives up after 60 s.
  options(timeout = 600)
  if (length(missing)) install.packages(missing, lib = '${DEPS}', repos = 'https://cloud.r-project.org')"
# --preclean so that no object file of an earlier build is reused, --clean so that the source stays as unpacked.
in_container R CMD INSTALL --preclean --clean -l "${LIBRARY}" "${SOURCE}" 2>&1 | tee "${LIBRARY}.install.log"
echo "compiled $(grep -c -E '^(g\+\+|gcc) .* -c ' "${LIBRARY}.install.log") source files"
in_container Rscript --vanilla -e "
  library(RiskMap, lib.loc = '${LIBRARY}')
  cat('RiskMap', as.character(packageVersion('RiskMap')), 'from', find.package('RiskMap'), '\n')"
