# TFpairminer

Mine candidate **functional TF pairs per cell type** from single-cell atlases, and browse them on a searchable site.

**Site:** https://double2strand.github.io/TFpairminer/ — search a TF (`FOXA2`), a pair (`FOXA2:HNF4A`) or a cell type (`cholangiocyte`).
- **Visual:** pair network, TF × cell-context heatmap, per-TF pages.
- **Scientific:** a sortable, filterable table with cell r, donor r, lift, detection, co-detected cells, dataset and collection support, prior flag, Pair_ID and PMIDs. It also shows the thresholds and the Census release, and has TSV download.
- **Run info:** provenance, which covers the manifest, the compendium hash, completion status and thresholds.

Layer 1 (this repo today) is CELLxGENE Census discovery. SCEA and Tabula Sapiens replication, plus scATAC evidence levels A–D, come next.

## Quickstart: 4 commands

```bash
git clone https://github.com/double2strand/TFpairminer.git && cd TFpairminer
./tfpm setup                 # venv + deps + unit tests + Census plan (compendium TSVs are in compendium/)
./tfpm mine --runner wr      # cluster: 100 shard jobs via wr + 1 dependent aggregate job
./tfpm status                # progress
./tfpm publish               # aggregate -> docs/data JSON -> git commit + push (site updates)
```

**Workstation** (for example 256 GB RAM): same, but `./tfpm mine --runner local --parallel 6`. That runs 6 shards at a time at about 32 GB each, in the background via `nohup`, so you can log out.

Notes
- **wr:** start a manager first (`wr manager start`, with `-s lsf` if needed). Jobs: `wr status -i tfpm-mine -o s` · failures: `wr status -i tfpm-mine -b` · retry: `wr retry -i tfpm-mine -b`. The aggregate job (`tfpm-aggregate`) waits on dep-group `tfpm_shards`. Optional settings: `TFPM_WR_QUEUE`, `TFPM_WR_MEM`, `TFPM_WR_TIME`, `TFPM_WR_CPUS`.
- **Settings:** env vars, or a `tfpm.conf` with `KEY=value` lines. `TFPM_OUT` (default `./run`, which should be on shared storage for wr), `TFPM_CENSUS` (`2025-11-08`), `TFPM_SHARDS` (100), `TFPM_TF_SET` (`exploratory`), `TFPM_PYTHON` (a Python 3.11 or 3.12 binary).
- **New compendium version:** `./tfpm setup /path/to/compendium.xlsx` re-exports `compendium/gene_register.tsv` and `interaction_priors.tsv`. Use a new `TFPM_OUT` so the plan provenance stays immutable.
- **Small test first:** `./tfpm mine --runner wr --pilot` runs 1 stratum per dataset with ≥20 cells.
- **Publish:** `publish` needs push rights in this clone (an SSH remote or a token). The site shows a **pilot** banner until `summary.json` reports `atlas_wide: true`.
- **Resuming:** shards are resumable. Re-running `mine` skips datasets that already have `_SUCCESS`.

## Statistics (v2)
The default primary statistic is a **metacell partial correlation**: within-donor metacells of 20 cells, adjusted for donor and log library size, with a t-test p-value and **BH q within each stratum**. Pseudobulk (donor-mean), residualised cell-level and the original raw cell-level r (`cell_raw`, v1) are also reported. Aggregation adds a **Fisher-z fixed and random-effects meta-analysis** across datasets with q-values, plus a **positive-control check** (MYC–MAX, JUN–FOS and others). Pass options through `tfpm`, for example `./tfpm mine --runner wr -- --primary-stat pseudobulk --max-q 0.1`, or set `TFPM_RUN_ARGS`. Details are in [mining/METHOD.md](mining/METHOD.md#v2-statistics-schema_version-2).

## Publishing and versioning
`./tfpm publish [--run DIR] [--label TEXT] [--no-latest]` runs these steps:
1. Validates the run.
2. Aggregates it if that hasn't been done.
3. Exports a versioned release to `docs/data/runs/<run_id>/`, sharded as per-TF (`tf/<SYMBOL>.json`) and per-cell-context (`ct/<i>.json`) files plus the top pairs, so it scales to the full atlas.
4. Updates `docs/data/runs/index.json`.
5. Commits and pushes.

The latest run is the site default. Older runs can be picked from the selector, or opened with `?run=<run_id>`. Data schema: `schema_version: 2`.

## Layout
| Path | What |
|---|---|
| `tfpm` | CLI wrapper |
| `mining/tf_census.py` | Census miner (`plan` / `run` / `aggregate`). The method is in [mining/METHOD.md](mining/METHOD.md) |
| `mining/test_tf_census.py` | Unit tests (`cd mining && python -m unittest -v`) |
| `scripts/export_compendium.py` | xlsx → `compendium/gene_register.tsv`, `interaction_priors.tsv` |
| `scripts/export_site.py` | Aggregate outputs → compact JSON in `docs/data/` |
| `docs/` | Static GitHub Pages site (vanilla JS) |
| `legacy/run_slurm.sh` | Old Slurm array script (not used) |

Candidates are co-expression associations, not verified interactions. Absence of a prior is not evidence that two TFs don't interact.
