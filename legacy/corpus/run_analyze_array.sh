#!/usr/bin/env bash
#SBATCH --job-name=hsc-pair-analysis
#SBATCH --array=0-15
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
set -euo pipefail

: "${HSC_OUT:?Set HSC_OUT to the RNA-first corpus plan directory}"
: "${HSC_DATA_DIR:?Set HSC_DATA_DIR to sample-id-named H5AD files}"
: "${HSC_PYTHON:?Set HSC_PYTHON to the installed Python executable}"
: "${SLURM_ARRAY_TASK_ID:?Run as a Slurm array}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../../mining" && pwd)"
cd "${PROJECT_DIR}"
"${HSC_PYTHON}" -m human_scatac_corpus_integration.pipeline analyze \
  --out "${HSC_OUT}" --data-dir "${HSC_DATA_DIR}" \
  --num-shards 16 --shard-index "${SLURM_ARRAY_TASK_ID}"
