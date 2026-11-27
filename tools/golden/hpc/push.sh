#!/bin/bash
# Copies a clean archive of each commit given, from the RiskMap checkout it is run in, to Wayland's
# GOLDEN_HOME/src/<sha> (GOLDEN_HOME is ~/storage/riskmap-golden there), and checks that the md5 sums of its files
# agree on both sides. A commit already there is only checked, since a commit's archive never changes.
#   tools/golden/hpc/push.sh <commit> [<commit> ...]
set -euo pipefail
HOST=${GOLDEN_HOST:-fronter2@wayland-2022.hec.lancaster.ac.uk}
REMOTE=${GOLDEN_REMOTE:-storage/riskmap-golden}
remote() { timeout 300 ssh -o BatchMode=yes -o ConnectTimeout=20 "${HOST}" "$@"; }
checksum='find . -type f -print0 | LC_ALL=C sort -z | xargs -0 md5sum | md5sum | cut -d" " -f1'
scratch=$(mktemp -d)
trap 'rm -rf "${scratch}"' EXIT
for ref in "$@"; do
  sha=$(git rev-parse --short=7 "${ref}^{commit}")
  mkdir -p "${scratch}/${sha}"
  git archive --format=tar "${sha}" | tar -x -C "${scratch}/${sha}"
  local_sum=$(cd "${scratch}/${sha}" && eval "${checksum}")
  target="${REMOTE}/src/${sha}"
  state=$(remote "if [ -d '${target}' ]; then echo present; else echo absent; fi") || state=unreachable
  if [ "${state}" = present ]; then
    echo "${sha} is already on Wayland"
  elif [ "${state}" != absent ]; then
    echo "cannot reach Wayland to check ${target}" >&2
    exit 1
  else
    git archive --format=tar "${sha}" |
      remote "mkdir -p '${target}.partial' && tar -x -C '${target}.partial' && mv '${target}.partial' '${target}'"
    echo "${sha} copied to ${target}"
  fi
  remote_sum=$(remote "cd '${target}' && ${checksum}")
  if [ "${local_sum}" != "${remote_sum}" ]; then
    echo "md5 of ${sha} differs: ${local_sum} here, ${remote_sum} on Wayland" >&2
    exit 1
  fi
  echo "${sha} md5 ${local_sum} on both sides"
done
