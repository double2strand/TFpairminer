#!/usr/bin/env bash
# Run against a pinned TFpairminer checkout on a remote node or scheduler.
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 plan|worker|aggregate" >&2
  exit 2
fi
: "${TFPM_REPO:?Set TFPM_REPO to a pinned TFpairminer checkout}"
: "${WORK_ROOT:?Set WORK_ROOT to the output root}"
: "${REGISTER:?Set REGISTER to gene_register.tsv}"
: "${CENSUS_VERSION:?Set CENSUS_VERSION to an explicit release date}"
: "${TABULA_DATASET_ID:?Set TABULA_DATASET_ID for All Cells in the pinned release}"

miner="$TFPM_REPO/mining/tf_census.py"
if [[ ! -f "$miner" ]]; then
  echo "TFpairminer miner not found: $miner" >&2
  exit 2
fi
if [[ ! -f "$REGISTER" ]]; then
  echo "Gene register not found: $REGISTER" >&2
  exit 2
fi

case "$1" in
  plan)
    python "$miner" plan --register "$REGISTER" --out "$WORK_ROOT/cellxgene" \
      --census-version "$CENSUS_VERSION" --tf-set exploratory \
      --exclude-collection-name "Tabula Sapiens"
    python "$miner" plan --register "$REGISTER" --out "$WORK_ROOT/tabula" \
      --census-version "$CENSUS_VERSION" --tf-set exploratory \
      --dataset-id "$TABULA_DATASET_ID"
    ;;
  worker)
    : "${SCREEN:?Set SCREEN=cellxgene or tabula}"
    : "${SHARD_INDEX:?Set SHARD_INDEX to a zero-based index}"
    : "${NUM_SHARDS:?Set NUM_SHARDS to the worker count}"
    if [[ "$SCREEN" != "cellxgene" && "$SCREEN" != "tabula" ]]; then
      echo "SCREEN must be cellxgene or tabula" >&2
      exit 2
    fi
    python "$miner" run --out "$WORK_ROOT/$SCREEN" \
      --shard-index "$SHARD_INDEX" --num-shards "$NUM_SHARDS" \
      --primary-stat metacell --max-q 0.05 --min-rho 0.25 \
      --min-lift 1.1 --min-codetected 25 --min-detection 0.05 \
      --min-cells 100 --top-k 0 --max-strata 0
    ;;
  aggregate)
    : "${SCREEN:?Set SCREEN=cellxgene or tabula}"
    if [[ "$SCREEN" != "cellxgene" && "$SCREEN" != "tabula" ]]; then
      echo "SCREEN must be cellxgene or tabula" >&2
      exit 2
    fi
    python "$miner" aggregate --out "$WORK_ROOT/$SCREEN"
    ;;
  *) echo "Unknown stage: $1" >&2; exit 2 ;;
esac
