#!/bin/bash
# Public entry point. Downloads/setup/preflight on login; all processing in SLURM.
set -euo pipefail
cd "$(dirname "$0")"
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
  cat <<'HELP'
Usage: bash run.sh [options]
  --config FILE         Configuration (relative to this project's root)
  --download-sources    Download missing HapMap3+/Catalog files, then exit
  --check-only          Check existing environment/resources; do not submit
  --chromosomes LIST    Autosomes, e.g. 22 or 1-3,22 (default: 1-22)
  --max-concurrent N    Maximum simultaneous chromosome jobs (default: 4)
Default: validate staged local references, then submit candidates -> chromosomes -> merge.
No large references are downloaded automatically. See README.md for staging details.
Environment: PYTHON_MODULE=python/3.13; HVS_TOOL_MODULES="bcftools picard".
HELP
  exit 0
fi
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "Invoke run.sh on the login node; processing is submitted to SLURM." >&2
  exit 1
fi
check_only=0
download_only=0
for arg in "$@"; do
  [[ "$arg" != --check-only ]] || check_only=1
  [[ "$arg" != --download-sources ]] || download_only=1
done
if type module >/dev/null 2>&1; then
  module load "${PYTHON_MODULE:-python/3.13}"
  if [[ "$download_only" == 0 ]]; then
    read -r -a tool_modules <<< "${HVS_TOOL_MODULES-bcftools picard}"
    for tool_module in "${tool_modules[@]}"; do module load "$tool_module"; done
  fi
fi
python3 -c 'import sys; assert sys.version_info >= (3, 11), "Requires Python 3.11 or newer"'
python_bin=env/harmonized/bin/python
if ! [[ -x "$python_bin" ]]; then
  if [[ "$check_only" == 1 ]]; then
    echo "Missing env/harmonized; run without --check-only to prepare it." >&2
    exit 1
  fi
  python3 -m venv env/harmonized
fi
[[ "$(python3 -c 'import sys; print(sys.version_info[:2])')" == \
   "$("$python_bin" -c 'import sys; print(sys.version_info[:2])')" ]] || {
  echo "Recreate env/harmonized with the selected PYTHON_MODULE." >&2
  exit 1
}
if [[ "$check_only" == 0 ]]; then
  "$python_bin" -m pip install --only-binary=:all: -r requirements.txt
fi
"$python_bin" -c 'import pandas, numpy, pysam, pyreadr'
"$python_bin" -m pip check
exec "$python_bin" -u workflow.py "$@"
