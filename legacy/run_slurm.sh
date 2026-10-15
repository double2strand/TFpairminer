#!/usr/bin/env bash
#SBATCH --job-name=tf-census
#SBATCH --array=0-99
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/tf-census-%A_%a.out

set -euo pipefail

: "${TF_CENSUS_ROOT:?Set TF_CENSUS_ROOT to the copied project directory}"
: "${TF_CENSUS_OUT:?Set TF_CENSUS_OUT to the shared plan/output directory}"
: "${TF_CENSUS_PYTHON:?Set TF_CENSUS_PYTHON to the virtualenv Python executable}"

"${TF_CENSUS_PYTHON}" "${TF_CENSUS_ROOT}/cellxgene_tf_mining/tf_census.py" run \
  --out "${TF_CENSUS_OUT}" \
  --num-shards 100 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --batch-size 1000 \
  --tmpdir "${SLURM_TMPDIR:-/tmp}"
