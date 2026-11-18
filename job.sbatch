#!/bin/bash
#SBATCH --job-name=hvs
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=05:00:00
# Internal worker, submitted by run.sh. All inputs must already be local.
set -euo pipefail
export HVS_ROOT="$1"
config_file="$2"
workflow_file="$3"
stage="$4"
chromosomes="$5"
cd "$HVS_ROOT"
module load "${PYTHON_MODULE:-python/3.13}"
read -r -a tool_modules <<< "${HVS_TOOL_MODULES-bcftools picard}"
for tool_module in "${tool_modules[@]}"; do module load "$tool_module"; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx4g"
: "${SLURM_TMPDIR:?SLURM must provide local scratch storage}"
export TMPDIR="$SLURM_TMPDIR"
args=(--config "$config_file" --stage "$stage" --chromosomes "$chromosomes")
if [[ "$stage" == chromosome ]]; then args+=(--chrom "$SLURM_ARRAY_TASK_ID"); fi
exec env/harmonized/bin/python -u "$workflow_file" "${args[@]}"
