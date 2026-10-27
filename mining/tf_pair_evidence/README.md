# ChromBPNet × two-scRNA TF-pair evidence layer

This stage takes **tested** `pairs_cbp.tsv` outputs from [TFpairminer's ChromBPNet module](https://github.com/double2strand/TFpairminer/tree/main/mining/motif_cooccurrence) and runs two independent, release-pinned RNA screens using its donor-aware Census miner. It then writes the four requested evidence views. These are not four mutually exclusive bins: a pair with both RNA hits appears in both ATAC+RNA views and in `06_atac_plus_both_rna.tsv`.

## 1. Remote RNA screens

Use a pinned TFpairminer checkout and one explicit CELLxGENE Census release for both screens. In a fresh Python 3.11/3.12 Linux environment, install that checkout's `mining/requirements.txt`. The script defaults to the exploratory TF set from the compendium, includes all complete strata, and uses within-donor metacell partial correlation with BH q≤0.05, r≥0.25, co-detection lift≥1.1 and ≥25 co-detected cells. It excludes the Tabula Sapiens collection from the broad discovery screen and uses only Tabula Sapiens **All Cells** for RNA2, so the two datasets do not overlap. Verify the Tabula dataset ID in the chosen Census release; `53d208b0-2cfd-4366-9866-c3c6114081bc` is the ID recorded in the local 2025-11-08 pilot.

```bash
export TFPM_REPO=/path/to/TFpairminer
export WORK_ROOT=/scratch/tf_pair_evidence
export REGISTER=/path/to/gene_register.tsv
export CENSUS_VERSION=2025-11-08
export TABULA_DATASET_ID=53d208b0-2cfd-4366-9866-c3c6114081bc
bash tf_pair_evidence/run_remote.sh plan

# Submit 100 separate jobs with SHARD_INDEX=0,...,99:
SCREEN=cellxgene NUM_SHARDS=100 SHARD_INDEX=0 bash tf_pair_evidence/run_remote.sh worker
# Tabula can use one job (or be sharded too):
SCREEN=tabula NUM_SHARDS=1 SHARD_INDEX=0 bash tf_pair_evidence/run_remote.sh worker

SCREEN=cellxgene bash tf_pair_evidence/run_remote.sh aggregate
SCREEN=tabula bash tf_pair_evidence/run_remote.sh aggregate
```

On Slurm, submit `worker` as an array for each screen and run `aggregate` only after all shards finish. The source checkout should be pinned to a commit (audited against `812cf1f338942739ec681deeebd93a8b34162a4e`). The two planned `datasets.tsv` files should be inspected to confirm that RNA1 excludes all Tabula views and RNA2 contains the intended All Cells dataset. Keep the whole planned run directory; the integrator checks manifest provenance, completion, no stratum or top-k cap, same Census release, same register and non-overlapping dataset IDs.

## 2. Register ChromBPNet outputs and context matches

Make an ATAC manifest with one row per **actual ChromBPNet model/context run**:

```tsv
atac_dataset	atac_context	pairs_cbp_path
scATAC_Corpus	cell_type_A	/absolute/path/to/pairs_cbp.tsv
ATAC_atlas_2	cell_type_B	/absolute/path/to/pairs_cbp.tsv
```

Make a separate context map with ontology-based RNA matches:

```tsv
atac_dataset	atac_context	cell_type_id	tissue_id	disease_id	mapping_level
scATAC_Corpus	cell_type_A	CL:0000000	UBERON:0000000	PATO:0000461	exact
ATAC_atlas_2	cell_type_B	CL:0000000	UBERON:0000000	PATO:0000461	proxy
```

Only `exact` matches enter the four primary views. `proxy` and `unmapped` contexts remain in `05_unassessed.tsv`. A cancer cell line should **not** be declared an exact match to a normal-tissue RNA cell type just because its name or lineage is similar. For those cases, obtain matched cell-line scRNA data or report the cross-modal result as exploratory.

Run:

```bash
python -m tf_pair_evidence.integrate \
  --atac-manifest /path/to/atac_manifest.tsv \
  --context-map /path/to/context_map.tsv \
  --register /path/to/gene_register.tsv \
  --rna1 /scratch/tf_pair_evidence/cellxgene \
  --rna2 /scratch/tf_pair_evidence/tabula \
  --out /scratch/tf_pair_evidence/integrated
```

## 3. Output meanings

| Output | Operational criterion |
| --- | --- |
| `01_atac_only.tsv` | ChromBPNet `status=tested`, `k=2`, `contribution_supported=yes`; both complete RNA screens cover the pair but have no qualifying hit. This means **no detected RNA hit under these settings**, not no transcriptional relationship. |
| `02_atac_plus_rna1.tsv` | ChromBPNet-supported plus significant CELLxGENE non-Tabula co-expression in the exact context. |
| `03_atac_plus_rna2.tsv` | ChromBPNet-supported plus significant Tabula Sapiens co-expression in the exact context. |
| `04_rna_replicated_cbp_no_support.tsv` | Significant in both RNA screens; ChromBPNet explicitly tested the pair but `contribution_supported=no`. This means **no support by the selected ChromBPNet criterion**, not proof of no chromatin interaction. |

`master.tsv` contains every evaluated context-pair, including missingness and synergy metadata. `05_unassessed.tsv` contains untested ATAC pairs, inadequate RNA coverage, proxy context matches, and inconclusive tests; do not relabel these as ATAC- or RNA-negative. `06_atac_plus_both_rna.tsv` is the overlap of views 2 and 3. `07_rna_replicated_atac_unassessed.tsv` is a priority queue for RNA-replicated pairs that still need a valid ChromBPNet test; it is **not** the requested RNA-only validated set. `modules.tsv` reports connected components with at least three TFs within each evidence view, ATAC dataset and exact cell context. A component is a network summary, not proof that all members form one physical complex. `rejected_cbp_labels.tsv` records ambiguous family/composite TF labels that cannot be mapped one-to-one to the gene register. `summary.json` records counts and input hashes.

RNA significance is based on `min_q_within_stratum`, not the aggregate random-effects meta q (which can be biased when non-significant datasets are not stored). Only metacell or pseudobulk hits qualify; a `cell_raw` fallback does not. A missing RNA hit is considered screened only if both genes pass detection in the same assay/context and the uncapped run is complete. Since TFpairminer preselects ChromBPNet candidate pairs, absence from `pairs_cbp.tsv` means **ATAC unassessed**, not ATAC-negative. The RNA-only output is accordingly limited to pairs explicitly tested by ChromBPNet.

## Current input audit

The existing local HepG2 chr22 output is a **motif co-occurrence pilot**, not a ChromBPNet `pairs_cbp.tsv`. The local CELLxGENE/Tabula pilot uses `cell_raw` or donor-direction thresholds without donor-aware q-values. The public TFpairminer exploratory ATAC view is explicitly marked `atac_catlas_proxy_fetal` and the repository does not publish a `pairs_cbp.tsv` result. Those files are deliberately rejected as validated inputs. The integration code is runnable and unit-tested, but the four biological result tables cannot yet be populated honestly from currently available files.

Run tests: `python -m unittest tf_pair_evidence.test_integrate -v`.

## Targeted ATAC + RNA layer for TFpairminer v2 (`targeted_rna.py`, `targeted_integrate.py`)

This is a box-scale alternative to the atlas-wide `integrate.py`. It tests only the ChromBPNet candidates (`pairs_cbp.tsv`, k=2 and k=3) for co-expression in a matched normal RNA context. It runs two screens separately:
- **CELLxGENE Census** for discovery. Datasets from Tabula Sapiens collections are excluded.
- **Tabula Sapiens** All Cells (`53d208b0…`) for replication.

**Statistic.** The test reuses the `tf_census` within-donor metacell partial r, with pseudobulk as the fallback.
- Strata are dataset × cell type × tissue × assay, with at least 100 cells.
- Both genes must be detected in ≥5% of cells.
- A hit needs r≥0.25, BH q≤0.05, lift≥1.1 and ≥25 co-detected cells.
- BH is applied only over the candidate pairs in each stratum. These q values are therefore less conservative than tf_census q values.

**Contexts.** Every context match is labelled **proxy**, because a cancer cell line is never an exact match to normal tissue.
- HepG2 → hepatocytes.
- K562 → erythroid, MEP and CMP progenitors. Mature erythrocytes are excluded.

**Rules.**
- A triplet is an RNA hit only when all three sub-pairs are hits.
- A pair or triplet is **unassessed** when any of these is true:
  - a sub-pair is untested (not detected or no stratum);
  - a label is composite;
  - a screen is missing;
  - the profile has no ChromBPNet output, which is true of every Corpus profile.
- Unassessed never counts as negative.
- The ATAC atlases (CCLA and Corpus) are never pooled.

Run: `./tfpm v2 evidence --runner wr|local [--only CCLA_HepG2]`. Output goes to `$TFPM_V2_OUT/rna/<screen>/<id>` (checkpointed per stratum) and `$TFPM_V2_OUT/atac_rna/`, which contains:
- `master.tsv`;
- `01_chrombpnet_only.tsv` … `04_rna_replicated_cbp_unsupported.tsv`;
- `05_unassessed.tsv`;
- `modules.tsv`;
- `summary.json`.

Tests: `python -m unittest discover -s tf_pair_evidence/tests -t .` (run from `mining/`).
