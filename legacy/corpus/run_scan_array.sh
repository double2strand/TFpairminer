#!/usr/bin/env bash
#SBATCH --job-name=hsc-motif-scan
#SBATCH --array=0-31
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
set -euo pipefail

: "${HSC_OUT:?Set HSC_OUT to the RNA-first corpus plan directory}"
: "${HSC_DATA_DIR:?Set HSC_DATA_DIR to sample-id-named H5AD files}"
: "${HSC_GENOME_FASTA:?Set HSC_GENOME_FASTA to an uncompressed hg38 FASTA}"
: "${HSC_PYTHON:?Set HSC_PYTHON to the installed Python executable}"
: "${SLURM_ARRAY_TASK_ID:?Run as a Slurm array}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../../mining" && pwd)"
cd "${PROJECT_DIR}"
"${HSC_PYTHON}" -m human_scatac_corpus_integration.pipeline scan \
  --out "${HSC_OUT}" --data-dir "${HSC_DATA_DIR}" \
  --genome-fasta "${HSC_GENOME_FASTA}" \
  --num-shards 32 --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --tmpdir "${SLURM_TMPDIR:-/tmp}"
