# CELLxGENE Census TF-expression discovery

This is the first, transcriptomic-discovery layer of the TF-pair toolkit. It takes the [v0.2 human TF/regulator gene register](../outputs/tf_compendium_v0_2/gene_register.tsv), measures expression of all mapped candidates across annotated human cell populations in the [CZ CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/), and screens co-expressed TF pairs. It does **not** infer binding, direct cooperation, or regulatory causality. Independent RNA-atlas replication and scATAC motif/co-accessibility checks remain separate downstream layers.

## Remote setup

Use Python 3.11 or 3.12 on a Linux VM or cluster node with outbound HTTPS/S3 access and enough shared storage for outputs. Copy this directory and `gene_register.tsv` to the remote system. A node-local scratch directory is recommended; 32 GB RAM per worker is a conservative starting allocation. The tested client version is pinned in `requirements.txt`.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r cellxgene_tf_mining/requirements.txt
.venv/bin/python -m unittest discover -s cellxgene_tf_mining -p 'test_*.py' -v
```

Run from the directory containing `cellxgene_tf_mining` and `outputs`. The `plan` stage resolves `stable` to an immutable Census build date, maps Ensembl IDs, and records source provenance. It does not download the complete expression atlas.

```bash
.venv/bin/python cellxgene_tf_mining/tf_census.py plan \
  --register outputs/tf_compendium_v0_2/gene_register.tsv \
  --out tf_census_run \
  --census-version stable \
  --tf-set exploratory
```

For a small, known-human pilot, add `--dataset-id 0c9a8cfb-6649-4d52-b418-6d8e56bd7afe` when making a **new** plan directory. Then run it with `--min-cells 20 --max-strata 1`. The ID is only an example; discoverability and metadata can change across releases. A plan directory cannot be overwritten, so a changed release or TF set requires a new directory.

The full run can be split into independent deterministic shards. All jobs share one output directory. Each dataset gets an `_SUCCESS` marker only after its outputs are complete; rerunning the same shard skips completed datasets. Different run thresholds for a completed dataset are rejected.

```bash
.venv/bin/python cellxgene_tf_mining/tf_census.py run \
  --out tf_census_run \
  --num-shards 100 \
  --shard-index 0 \
  --batch-size 1000 \
  --tmpdir /path/to/node-local-scratch
```

Submit shard indices `0` through `99` in parallel, or use [run_slurm.sh](run_slurm.sh). The dataset catalog is Census-wide (including datasets with no human primary cells); such datasets produce empty outputs and are counted as completed, not as evidence of absent expression. Large datasets with many cell-type strata can take hours. Lower `--batch-size` if memory is limited, increase it if network round trips dominate. Start with one pilot shard and inspect runtime before scaling.

For the included Slurm example, create the log directory **before** submission and export absolute paths:

```bash
mkdir -p logs
export TF_CENSUS_ROOT="$PWD"
export TF_CENSUS_OUT="$PWD/tf_census_run"
export TF_CENSUS_PYTHON="$PWD/.venv/bin/python"
sbatch cellxgene_tf_mining/run_slurm.sh
```

Adjust the `#SBATCH --array=0-99` line and matching `--num-shards 100` together if using a different shard count.

Once the intended shards finish:

```bash
.venv/bin/python cellxgene_tf_mining/tf_census.py aggregate \
  --out tf_census_run \
  --interaction-priors outputs/tf_compendium_v0_2/interaction_priors.tsv
```

The `--interaction-priors` input is optional. It labels Census pair hits with existing physical-association evidence from the compendium; absence of a recorded prior is **not** evidence that the proteins do not interact. Check `aggregate/summary.json`: `complete: true` means every **planned** dataset has a completion marker; `atlas_wide: true` additionally requires an unrestricted dataset plan and no per-dataset stratum cap. Aggregating an incomplete plan is allowed for inspection, but the summary marks it incomplete. The aggregate database is intentionally not overwritten; move the old `aggregate/` aside or use a new plan directory before re-aggregating.

## Exact analysis method

1. Read the compendium's Ensembl IDs, resolve them against the pinned human RNA feature table, and report unmapped IDs. All mapped compendium genes get expression profiles; pair screening uses the selected TF set (`strict`, `exploratory`, or explicitly `expanded` regulator set).
2. For each Census dataset, query only `Homo sapiens` non-spatial RNA cells marked `is_primary_data == True` and with positive `obs.raw_sum`. This prevents duplicate-cell reuse and undefined library-size normalization. Record the Census feature-presence flag separately; genes absent from this flag are not treated as measured zeros. The Census schema's exact interpretation of this flag is described [here](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_schema.html).
3. Keep cell-type ontology ID, general-tissue ontology ID, disease ontology ID, and assay ontology ID separate. No cells from unlike disease states, tissues, or assay types are pooled inside a stratum. Donor IDs are retained for donor-level diagnostics. A stratum must contain at least `--min-cells` primary cells (default 100).
4. Read raw counts in bounded cell batches. For each gene, calculate detection fraction, mean raw count, and mean `log1p(10,000 × raw count / raw_sum)`. These are **within-dataset** summaries; the default assay-stratified output should not be interpreted as a batch-corrected cross-study expression matrix.
5. For each TF pair within a stratum, calculate Pearson correlation of log-normalized cell-level expression and co-detection lift `n_both × n_cells / (n_detect_TF1 × n_detect_TF2)`. Default screening requires each TF in ≥5% of cells, Pearson `r ≥ 0.25`, lift `≥1.1`, and ≥25 co-detected cells. Every passing pair is retained by default; `--top-k 25` is an optional per-TF output cap. If ≥4 donors each contribute ≥10 cells, also report Pearson correlation across donor mean-expression profiles. This donor statistic is descriptive and is not used as an automatic significance test.
6. Aggregate gene profiles over datasets sharing the same ontology-defined stratum. Rank screened pairs by the number of distinct source **collections** in which they pass, followed by cell-count-weighted mean correlation. A pair supported by ≥2 collections is a stronger follow-up than a single-dataset hit, but the table does **not** estimate an FDR or formal cross-study meta-analysis.

If comparing platform technologies, use repeated `--assay-id` to restrict to comparable UMI assays (for example, `EFO:0009899` is 10x 3′ v2 in the tested release). The default processes all assays **separately**; raw-count library-size scaling does not length-correct full-gene protocols such as Smart-seq. Even within an assay, cell-level correlations can reflect cell state, library depth, or hidden covariates. The resulting pairs are *candidate co-expression associations*, not verified heterodimers or functional cooperativity.

This first pass does not stratify developmental stage, sex, or treatment, and it treats any multi-valued disease annotation as one exact source string. Add those covariates or use donor-level models for a focused biological comparison before claiming cell-type-specific cooperation.

## Files

| File | Meaning |
|---|---|
| `manifest.json` | Resolved Census release, compendium hash, mapping counts and analysis defaults. |
| `features.tsv`, `unmapped.tsv`, `datasets.tsv` | Gene and dataset mapping/provenance. |
| `datasets/<dataset_id>/feature_availability.tsv.gz` | Census feature-presence flag for each mapped compendium gene. |
| `datasets/<dataset_id>/profiles.tsv.gz` | One gene × study × cell-type/tissue/disease/assay stratum. |
| `datasets/<dataset_id>/pairs.tsv.gz` | Within-stratum screened TF-pair hits, with cell and optional donor correlations. |
| `aggregate/cell_type_profiles.tsv` | Cross-dataset expression summary for each gene and stratum. |
| `aggregate/cell_type_pair_candidates.tsv` | TF-pair hits with study/collection support and optional physical-prior annotation for the next validation layer. |
| `aggregate/summary.json` | Completion status and run thresholds. |

Source APIs and rationale: [Census Python API](https://chanzuckerberg.github.io/cellxgene-census/python-api.html), [query/extraction tutorial](https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_query_extract.html), [dataset feature-presence guide](https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_dataset_presence.html), [duplicate-cell guidance](https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_duplicated_cells.html). The pipeline uses public Census data; no user credentials are required.

## Local validation performed

The unit tests exercise normalization, pair statistics, gzip output, and aggregation. On Census release `2025-11-08`, `plan` mapped 3,868 of 3,874 compendium rows, including 1,542 exploratory TFs. Three source rows had no Ensembl ID and three Ensembl IDs were absent from the pinned Census; all six are reported in `unmapped.tsv`. An end-to-end live pilot on one healthy-liver cholangiocyte stratum (375 cells) generated 3,864 gene-profile rows and 137 pair hits with pilot thresholds and no output cap; 10 hits had a recorded physical-interaction prior. This is an **API smoke test**, not an atlas-wide candidate result; the full run has not been launched here.

## v2 statistics (schema_version 2)

The cell-level screen above is kept as `--primary-stat cell_raw`, which reproduces v1 selection exactly. The default is now `--primary-stat metacell`. Within each stratum, `run` streams sufficient statistics and computes four correlation estimates for every TF pair:

| Statistic | Units | Adjustment | df |
|---|---|---|---|
| `cell_pearson_r` (v1) | cells | none | n − 2 |
| `resid_cell_r` | cells | donor fixed effects (within-donor centring; unknown donors form one group) + log library size | n − G − 2 |
| `metacell_r` | within-donor metacells of `--metacell-size` (20) consecutive cells; metacells with fewer than half that many cells are dropped | donor fixed effects + mean log library size | m − G − 2 |
| `pseudobulk_r` (`donor_pearson_r`) | donor means (donors with ≥ `--min-donor-cells`, at least 4 donors) | log library size | D − 3 |

Dataset and assay are constant within a stratum, so they are controlled by design. The primary statistic gets a two-sided t-test p-value. If it can't be computed (for example, too few metacells or donors), it falls back down the chain pseudobulk → metacell → resid_cell → cell_raw, and `primary_stat` records which one was used. Benjamini–Hochberg q-values are computed over **all** pairs tested in that stratum, meaning both TFs have detection ≥ `--min-detection` (`n_tested_pairs`). A pair passes when primary r ≥ `--min-rho`, q ≤ `--max-q`, lift ≥ `--min-lift`, and co-detected cells ≥ `--min-codetected`. Non-passing tested pairs with p ≤ `--store-p` and r ≥ `--store-min-r` are also stored (`passed=False`) for the meta-analysis.

`aggregate` performs a Fisher-z meta-analysis for each pair × stratum across datasets. It weights each dataset by df − 1 and reports fixed-effect and DerSimonian–Laird random-effects r and p, τ² and I², then BH q-values over random-effects p within the stratum. Limitation: datasets where a pair was not stored, because it wasn't significant, are missing. The meta-statistics are therefore optimistic, and collection support remains the primary replication signal. `aggregate/positive_controls.tsv` checks known partners (MYC–MAX, JUN–FOS, JUNB–FOS, JUN–ATF3, HIF1A–ARNT, CLOCK–BMAL1, STAT1–STAT2, SMAD3–SMAD4, RELA–NFKB1, TCF3–TCF12, HNF1A–HNF4A, FOXA2–HNF4A, SPI1–CEBPA, GATA1–TAL1, PBX1–MEIS1, SOX2–POU5F1, MAX–MXD4). `summary.json` reports how many were recovered. Co-expression is not expected for every physical partner, so treat this as a sanity check, not a benchmark.

Unit tests cover the BH reference values, agreement of p-values with scipy, equivalence between the streamed residualised correlation and an explicit regression, recovery of a simulated donor-confounded true pair, rejection of independent TFs, FDR calibration under the null, and meta-analysis behaviour.
