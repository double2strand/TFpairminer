# CATlas integration for replicated TF pairs

This is the third, independent evidence layer of the TF-pair screen. It takes the
`candidate_summary.tsv.gz` output of Tabula Sapiens replication, maps each RNA
cell-type context to CATlas clusters, and asks whether **both TF sequence motifs
are enriched in accessible distal cCREs in the same mapped CATlas cluster**.
An inaccessible distal-cCRE control set is matched by GC decile and width
quartile within each cluster. The pipeline is designed to run on a remote Linux
machine or Slurm cluster. Its stages are `fetch → plan → scan → compile →
analyze → aggregate`.

## Inputs and provenance

`fetch` obtains three pinned bulk inputs and checks SHA-256 before exposing them:

| Input | Source | SHA-256 |
| --- | --- | --- |
| 222-cluster × 1,154,611-cCRE binary hg38 H5AD | [CATlas-derived Genentech mirror](https://huggingface.co/datasets/Genentech/binary-atac-tutorial-data/tree/6b6d749d445605f2cdb98e5ca83018b63271b3c4) | `bbe4343c308c83efee0c4ec085229c82a6a7077df1bd66cb4096801a4578aa60` |
| Cell Ontology table | [CATlas official table](https://decoder-genetics.wustl.edu/catlasv1/humanenhancer/data/Cell_ontology.tsv) | `7ef6945a19298bd24c12f42768df64d5d1374029b28755f551f7cd0f14ad1062` |
| CORE vertebrate non-redundant MEME profiles | [JASPAR 2026](https://jaspar2026.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt) | `ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15` |

The H5AD mirror is not the official CATlas hosting service; the source paper
and its cell-type ontology are the reference for biological interpretation.
The legacy CATlas bulk matrix links currently return HTTP 404, so this mirror
is pinned and checksummed. Retain `input_manifest.json` and the run's
`manifest.json` for provenance. A valid **uncompressed hg38/GRCh38 FASTA with
`chr`-prefixed contigs** is additionally required for motif scanning. Obtain it
from the [GENCODE primary-assembly FASTA release](https://www.gencodegenes.org/human/)
or an equivalent documented hg38 source; keep its checksum in your project log.
Do not use GRCh37 or a `1`/`2`/`X`-named FASTA without reconciliation.

## Install on a remote machine

The Python versions in `requirements.txt` are those used for the local synthetic
test. `bedtools` and MEME Suite's `fimo` are external executables.

```bash
conda create -n catlas_pairs -c conda-forge -c bioconda python=3.11 bedtools meme
conda activate catlas_pairs
python -m pip install -r catlas_tf_integration/requirements.txt
python -m unittest discover -s catlas_tf_integration -p 'test_*.py' -v
```

## Full run

Run from the parent directory of `catlas_tf_integration`, substituting your
own absolute input/output paths. The `candidate_summary.tsv.gz` may contain
many cell types and candidate pairs; only `replicated_donor_supported` rows are
selected by default. The work is sharded across the CATlas cCREs for scanning
and across clusters for analysis; each expensive motif scan is shared across
all candidate pairs.

```bash
python catlas_tf_integration/catlas_pipeline.py fetch --out /data/catlas_inputs

python catlas_tf_integration/catlas_pipeline.py plan \
  --candidate-summary /data/ts/aggregate/candidate_summary.tsv.gz \
  --catlas-h5ad /data/catlas_inputs/catlas_binary_atac.h5ad \
  --ontology /data/catlas_inputs/catlas_cell_ontology.tsv \
  --meme /data/catlas_inputs/jaspar2026_core_vertebrates.meme \
  --out /data/catlas_pairs_run
```

First inspect `mapping_review.tsv`, `clusters.tsv`, `motif_map.tsv` and
`manifest.json`. A cell-type ID must match exactly before automatic mapping;
such ontology-only matches are labelled **provisional**. For a claim of strict
CATlas confirmation, provide a reviewed, tissue- and stage-consistent
crosswalk and make a **new plan directory**:

```tsv
catlas_cluster	cell_type_id	tissue_id	stage	mapping_tier	reviewed
Hepatocyte	CL:0000182	UBERON:0002107	adult	exact_reviewed	Yes
```

That row is only an illustration: verify the ontology ID, tissue and CATlas
cluster before use. Supported tiers are `exact_reviewed`, `proxy_fetal` and
`proxy_related_celltype`; proxies are never called exact validation. The
pipeline does not silently infer a proxy from a similar name. Disease-context
RNA candidates cannot be called exact matches to this healthy CATlas atlas;
automatic ontology matches for them are marked `proxy_disease_context`. If
the manifest reports `n_mapped_rows: 0`, skip scanning and run `aggregate`
immediately; the output will record `no_exact_celltype` or
`ontology_present_cluster_unresolved` for each screened pair.

For mapped candidates, run the remaining stages. On Slurm, set
`CATLAS_OUT`, `CATLAS_H5AD`, `CATLAS_GENOME_FASTA`, and `CATLAS_PYTHON` to
absolute paths and submit the supplied array scripts. Only proceed to the next
stage after the previous array has completed successfully.

```bash
sbatch --export=ALL,CATLAS_OUT=/data/catlas_pairs_run,CATLAS_H5AD=/data/catlas_inputs/catlas_binary_atac.h5ad,CATLAS_GENOME_FASTA=/data/hg38.fa,CATLAS_PYTHON=/path/to/env/bin/python catlas_tf_integration/run_scan_array.sh

python catlas_tf_integration/catlas_pipeline.py compile --out /data/catlas_pairs_run

sbatch --export=ALL,CATLAS_OUT=/data/catlas_pairs_run,CATLAS_H5AD=/data/catlas_inputs/catlas_binary_atac.h5ad,CATLAS_PYTHON=/path/to/env/bin/python catlas_tf_integration/run_analyze_array.sh

python catlas_tf_integration/catlas_pipeline.py aggregate --out /data/catlas_pairs_run
```

On a single node, replace each Slurm array with a loop over shard indices
(`scan --num-shards 32 --shard-index i`, then `analyze --num-shards 16
--shard-index i`). The job scripts default to 32 scan and 16 analysis shards;
change both the `#SBATCH --array` range and corresponding `--num-shards`
argument together. Shards have completion markers and can be rerun after a
failure. Do not reuse a plan directory if the candidate set, crosswalk,
reference data, or thresholds change. Allow several GB of shared scratch;
the complete genome FASTA and FIMO sequence scans are the costly steps.

## Test and interpretation

`aggregate/catlas_pair_results.tsv.gz` contains one row per candidate ×
matched CATlas cluster. `aggregate/summary.json` reports coverage and result
tiers. We test each TF's motif presence and joint same-cCRE motif presence in
matched foreground versus control cCREs using one-sided Fisher tests. BH FDR
is calculated across all tested cluster–TF and cluster–pair hypotheses.
`catlas_coenriched` requires a reviewed exact crosswalk, both TFs and the
joint motif hit enriched at `q ≤ 0.05` and odds ratio > 1, at least 20 joint
accessible cCREs, and greater joint motif lift in accessible than control
cCREs. `provisional_*` and `proxy_*` use the same motif criteria but are
separate, weaker mapping tiers. `not_coenriched` means the test ran and the
predefined screen was not met; `no_exact_celltype`,
`ontology_present_cluster_unresolved`, `motif_unavailable`,
`shared_motif_ambiguous`, `inconclusive_sparse_motif`, and
`insufficient_matched_cCREs` are **not negative biological results**.

The test is a prioritization screen, not biochemical proof that two proteins
bind together or simultaneously occupy a cCRE. JASPAR family members can
have nearly indistinguishable motifs, and two motif calls may overlap at one
sequence site. The joint test is **same-cCRE sequence-model co-enrichment**,
not a TF–TF interaction test; positive AP-1-family pairs in particular should
be checked with motif-similarity filtering, footprint/ChIP data, perturbation,
or cis-regulatory assays. This pipeline uses distal cCREs, not promoters, and
does not assign cCREs to target genes or establish expression causality.

## Current pilot coverage

The supplied strict Tabula Sapiens pilot has 20 donor-supported pairs from
adult intrahepatic cholangiocytes (`CL:0002538`). CATlas has **fetal**
cholangiocytes (`CL:1000488`) but no exact adult intrahepatic-cholangiocyte
cluster in this matrix. The current pilot therefore has **0 exact mapped rows**.
Its valid result is a coverage gap, not CATlas rejection or confirmation of
the 20 pairs. A fetal cholangiocyte analysis could be run separately with an
explicit `proxy_fetal` crosswalk, but it would answer a different developmental
question and remain labelled exploratory.
