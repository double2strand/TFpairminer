#!/usr/bin/env bash
#SBATCH --job-name=catlas-scan
#SBATCH --array=0-31
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
set -euo pipefail

: "${CATLAS_OUT:?Set CATLAS_OUT to the plan directory}"
: "${CATLAS_H5AD:?Set CATLAS_H5AD to catlas_binary_atac.h5ad}"
: "${CATLAS_GENOME_FASTA:?Set CATLAS_GENOME_FASTA to an uncompressed hg38 FASTA}"
: "${CATLAS_PYTHON:?Set CATLAS_PYTHON to the installed Python executable}"
: "${SLURM_ARRAY_TASK_ID:?Run as a Slurm array}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${CATLAS_PYTHON}" "${SCRIPT_DIR}/../mining/catlas_pipeline.py" scan \
  --out "${CATLAS_OUT}" \
  --catlas-h5ad "${CATLAS_H5AD}" \
  --genome-fasta "${CATLAS_GENOME_FASTA}" \
  --num-shards 32 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --tmpdir "${SLURM_TMPDIR:-/tmp}"
