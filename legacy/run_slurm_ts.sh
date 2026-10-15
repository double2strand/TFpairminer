#!/usr/bin/env bash
#SBATCH --job-name=tf-ts-validate
#SBATCH --array=0-99
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/tf-ts-%A_%a.out

set -euo pipefail
: "${TF_CENSUS_ROOT:?Set absolute project root}"
: "${TS_VALIDATION_OUT:?Set absolute validation output path}"
: "${TF_CENSUS_PYTHON:?Set absolute Python executable path}"

"${TF_CENSUS_PYTHON}" "${TF_CENSUS_ROOT}/cellxgene_tf_mining/tabula_sapiens_validate.py" run \
  --out "${TS_VALIDATION_OUT}" \
  --num-shards 100 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --batch-size 1000 \
  --tmpdir "${SLURM_TMPDIR:-/tmp}"
