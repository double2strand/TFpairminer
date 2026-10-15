# Tabula Sapiens validation of mined TF pairs

`tabula_sapiens_validate.py` tests **every row** of the CELLxGENE discovery table `aggregate/cell_type_pair_candidates.tsv` against a distinct Tabula Sapiens study. It does not search all pair combinations again. Candidate pairs are matched by Ensembl gene IDs and exact cell-type, general-tissue, and disease ontology IDs; Tabula Sapiens assays are evaluated separately. This is biological replication in a different study, but both studies are accessed through the same CELLxGENE Census portal—not independent-portal replication. The pinned release is inherited from discovery.

Use a discovery run that excludes the complete `Tabula Sapiens` collection. The validator also rejects a discovery run with a completed Tabula Sapiens dataset that contributed profiles. Within the validation collection, it uses only `Tabula Sapiens - All Cells` and `is_primary_data == True`. The organ-specific Tabula Sapiens views in the tested Census release reuse those cells and are marked non-primary; adding them would double-count observations. See the official [dataset-presence example](https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_dataset_presence.html) and [Census schema](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_schema.html).

## Run remotely

Install the [shared requirements](requirements.txt) as in the main README. Copy the discovery directory, including its `manifest.json`, `datasets.tsv`, completed-dataset markers, `aggregate/summary.json`, and pair-candidate TSV. Source and validation paths should be on shared storage. From the project directory:

```bash
python cellxgene_tf_mining/tabula_sapiens_validate.py plan \
  --discovery-run tf_census_run \
  --out ts_validation_run

python cellxgene_tf_mining/tabula_sapiens_validate.py run \
  --out ts_validation_run \
  --num-shards 100 --shard-index 0 \
  --batch-size 1000 --tmpdir /path/to/node-local-scratch

python cellxgene_tf_mining/tabula_sapiens_validate.py aggregate \
  --out ts_validation_run
```

Submit shard indices `0`–`99`, not only shard `0`, for a full run; [run_slurm_ts.sh](run_slurm_ts.sh) is an example. Before `sbatch`, create `logs/` and export absolute `TF_CENSUS_ROOT`, `TS_VALIDATION_OUT`, and `TF_CENSUS_PYTHON` paths, as shown for discovery in the main README. Shards are deterministic by cell-type/tissue/disease context, so all pairs for one context share an expression read and are not repeatedly downloaded. Per-context `.done` files make interrupted shards resumable. The aggregate reports its completion fraction; never interpret a partial aggregate as all-candidate coverage. `--assay-id` can be repeated during `plan` to restrict validation to specified assay ontology IDs. For a one-context pilot, use one shard (`--num-shards 1 --shard-index 0`). A fresh output directory is required for a new plan or a changed threshold set.

## Decision rule and interpretation

For each matched Tabula Sapiens stratum, raw counts are normalized as `log1p(10,000 × count / obs.raw_sum)`. Every candidate gets its cell-level Pearson correlation, individual detection fractions, co-detected cell count, and co-detection lift. A pair passes the **cell-level** screen when each gene is detected in at least 5% of cells, at least 25 cells co-detect it, `r ≥ 0.25`, and lift `≥ 1.1` (or stricter thresholds inherited from the discovery aggregate). The default validation stratum needs at least 100 primary cells. These cutoffs are operational screens, not p-values or controlled FDR.

For a stronger `replicated_donor_supported` call, the cell-level screen must pass and at least two donors with at least 20 cells each must have calculable within-donor correlations, with at least 75% of those donors showing positive correlation. `cell_level_support_only` means the pooled-cell screen passes without that donor consistency. `not_replicated` means the pair was measurable but fell below cell-level `r` or lift thresholds. `inconclusive_low_expression` means detection/co-detection was insufficient; it is **not** evidence against the pair. Separate statuses identify unavailable genes, insufficient cells, and missing matched cell-type/tissue/disease/assay strata. A donor-mean correlation is also reported when at least four sufficiently sampled donors exist, but it is descriptive, not part of the call.

The source and validation assays need not match by default; cross-platform effects and cell-state composition can alter correlations. Use an assay filter or review assay-specific rows when comparing a focused biological subset. Statistical uncertainty, donor dependence, compositional confounding, and multiple testing still require follow-up modeling. Co-expression alone does not establish protein interaction or cooperative binding.

## Outputs and scaling

| File | Meaning |
|---|---|
| `manifest.json` | Pinned release, Tabula Sapiens dataset ID, candidate-table hash, thresholds, and discovery completeness. |
| `features.tsv` | Which candidate Ensembl genes are present in the Tabula Sapiens dataset; absence is distinct from measured zero. |
| `contexts.tsv`, `candidates.sqlite` | Exact ontology contexts and indexed input candidates. |
| `groups/<context_id>.tsv.gz` | One row per pair × Tabula Sapiens assay, including untestable pairs and reasons. |
| `aggregate/replication_results.tsv.gz` | All assay-specific rows together. |
| `aggregate/candidate_summary.tsv.gz` | One best-supported assay row per discovery candidate, with numbers of assays tested/reported. |
| `aggregate/summary.json` | Planned/completed contexts, candidate counts, and status tallies. |

Input candidate count is not hard-coded. The plan indexes all rows in SQLite; each shard streams only the primary-cell metadata and the genes required by its assigned contexts. A full-atlas discovery may produce many contexts, so choose shard count and node-local scratch size according to the number of contexts and genes; run a small shard first to estimate memory and Census I/O. More shards do not divide a *single* very large context; if one dominates, partition that context's candidate-gene computation in a future extension. Preserve the manifest and all source files for provenance when moving runs between machines.
