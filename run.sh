#!/bin/bash
# Public entry point. Downloads/setup/preflight on login; all processing in SLURM.
set -euo pipefail
cd "$(dirname "$0")"
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
  cat <<'HELP'
Usage: bash run.sh [options]
  --time HH:MM:SS      Wall time per job (default: 01:00:00)
  --account ACCOUNT     SLURM billing account (default: config/environment/job.sbatch)
  --config FILE         Configuration (relative to this project's root)
  --download-sources    Download missing HapMap3+/Catalog files, then exit
  --download-references Download reference files/indexes, then exit (no jobs)
  --setup-only          Create environment, install/check dependencies, then exit
  --check-only          Check existing environment/resources; do not submit
  --status              Show controller progress and last error; do not submit
  --chromosomes LIST    Autosomes, e.g. 22 or 1-3,22 (default: saved run, otherwise 1-22)
  --max-concurrent N    Maximum simultaneous chromosome jobs (default: 2)
  --download-workers N  Parallel slices per chromosome download (default: 4)
  --scratch-dir DIR     Raw gnomAD download storage (default: /scratch/szabad)
Default: start/resume the login-node download controller and offline SLURM analysis jobs.
The login controller downloads to scratch; offline chromosome jobs cache counts and clean up VCFs.
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
setup_only=0
status_only=0
for arg in "$@"; do
  [[ "$arg" != --check-only && "$arg" != --status ]] || check_only=1
  [[ "$arg" != --download-sources && "$arg" != --download-references ]] || download_only=1
  [[ "$arg" != --setup-only ]] || setup_only=1
  [[ "$arg" != --status ]] || status_only=1
done
if (( check_only + download_only + setup_only > 1 )); then
  echo "Use only one of --setup-only, --check-only, --download-sources or --download-references." >&2
  exit 1
fi
if type module >/dev/null 2>&1; then
  module load "${PYTHON_MODULE:-python/3.13}"
  if [[ "$status_only" == 0 && ( "$download_only" == 0 || " $* " == *" --download-references "* ) ]]; then
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
  # Older pip releases may not recognize the cluster's compatible wheel tags.
  "$python_bin" -m pip install --upgrade pip
  # Prefer wheels, but permit pyreadr source builds on platforms without one.
  if ! "$python_bin" -m pip install --prefer-binary --only-binary=numpy,pandas,pysam -r requirements.txt; then
    echo "Dependency installation failed. pyreadr source builds require a C compiler, Python headers and zlib/bzip2/lzma development headers." >&2
    echo "If pip still reports no matching distribution, check access to your configured package index and its pyreadr packages." >&2
    exit 1
  fi
fi
"$python_bin" -c 'import pandas, numpy, pysam, pyreadr'
"$python_bin" -m pip check
if [[ "$setup_only" == 1 ]]; then
  # Validate tools without requiring staged reference data or a completed config.
  "$python_bin" - "$@" <<'PYTHON'
import json
from pathlib import Path
import shutil
import sys

args = sys.argv[1:]
config = Path('config.json')
for index, arg in enumerate(args):
    if arg == '--config':
        if index + 1 == len(args):
            raise SystemExit('--config requires a file')
        config = Path(args[index + 1]).expanduser()
    elif arg.startswith('--config='):
        config = Path(arg.split('=', 1)[1]).expanduser()
from workflow import picard_command
command = picard_command(json.loads(config.read_text()))
for tool in ('bcftools', command[0], 'sbatch', 'curl'):
    if not shutil.which(tool):
        raise SystemExit(f'Missing executable: {tool}. Load its cluster module or add it to PATH, then rerun setup.')
PYTHON
  echo "Environment ready: env/harmonized"
  exit 0
fi
exec "$python_bin" -u workflow.py "$@"
