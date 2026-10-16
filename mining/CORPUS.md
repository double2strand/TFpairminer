# RNA-first Human-scATAC-Corpus integration

This is a **parallel ATAC evidence layer**, not a CATlas-filtered stage:

```text
CELLxGENE Census TF-pair screen
              ↓
Tabula Sapiens RNA replication → fixed set of replicated candidates
              ├─ CATlas motif/accessibility screen
              └─ Human-scATAC-Corpus motif/accessibility screen
```

The corpus pipeline reads the Tabula Sapiens `candidate_summary.tsv.gz`
**directly** and selects `replicated_donor_supported` by default. It never
reads CATlas results or requires a CATlas positive hit. CATlas is used only as
a provenance label: a corpus sample derived from Zhang2021 or Domcke2020
overlaps the primary studies behind the CATlas atlas and cannot be counted as
an independent experiment. A pair can therefore get corpus evidence despite
having no CATlas cell-type match, but only if the corpus itself contains a
reviewed matching cell type and sufficient cells.

[Human-scATAC-Corpus](https://health.tsinghua.edu.cn/human-scatac-corpus/)
currently describes 5,407,621 cells across 35 datasets, with per-sample
cell-by-cCRE H5AD matrices and cell-level metadata. Its public [download
page](https://health.tsinghua.edu.cn/human-scatac-corpus/download.php) provides
individual and grouped downloads. The website does **not** advertise a stable
programmatic API; `catalog` reads its present browse-page data endpoint once
and saves a frozen TSV. If that endpoint changes, use a previously saved
catalog or a manually reviewed TSV with the same columns. Do not download all
5.4 million cells without selecting relevant tissues/studies first.

## Setup

Run from the directory containing both `catlas_tf_integration` and
`human_scatac_corpus_integration`. The corpus module reuses tested motif,
matching and FDR helpers from the former; it does **not** reuse CATlas data or
CATlas-selected pairs.

```bash
conda create -n tf_atac -c conda-forge -c bioconda python=3.11 bedtools meme
conda activate tf_atac
python -m pip install -r catlas_tf_integration/requirements.txt
python -m unittest human_scatac_corpus_integration.test_pipeline -v
```

Download the [JASPAR 2026 CORE vertebrate non-redundant MEME file](https://jaspar2026.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt)
for `--meme`; the locally verified SHA-256 was
`ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15`.
Supply an **uncompressed, chr-prefixed hg38/GRCh38 genome FASTA** for FIMO
scanning. Record its checksum. The HSC cCRE matrix stores coordinates as
`chr:start-end`, interpreted as BED-style half-open intervals; check the
reference FASTA/contig names before submitting large jobs. The pipeline
checks that all mapped samples have exactly the same ordered cCRE universe.

## 1. Catalog, select samples, and index cell types

```bash
python -m human_scatac_corpus_integration.pipeline catalog \
  --out /data/hsc/catalog.tsv

python -m human_scatac_corpus_integration.pipeline recommend \
  --candidate-summary /data/ts/aggregate/candidate_summary.tsv.gz \
  --catalog /data/hsc/catalog.tsv \
  --out /data/hsc/sample_recommendations.tsv

python -m human_scatac_corpus_integration.pipeline fetch \
  --catalog /data/hsc/catalog.tsv --sample-id 11 \
  --out /data/hsc/samples

python -m human_scatac_corpus_integration.pipeline index \
  --catalog /data/hsc/catalog.tsv \
  --data-dir /data/hsc/samples \
  --out /data/hsc/celltypes.tsv
```

`catalog.tsv` records the sample, dataset, tissue, organ, PMID and direct cCRE
matrix URL. `recommend` lists portal samples in tissues named by the
RNA-replicated candidates; it is a sample-download suggestion, **not** a
cell-type match or additional pair-selection step. Review its source-study
overlap and approximate file size before fetching. Each downloaded matrix is
named `<sample_id>.h5ad`;
`download_manifest.json` records its computed SHA-256 and source URL. The
`index` command reads H5AD metadata in backed mode and reports each label's
cell count, number of known donors and batches, developmental-stage hint,
cCRE-universe checksum, and source-study overlap. `0` known donors means
donor metadata were absent/unknown, **not** that a sample has zero donors.
Select further sample IDs based on the RNA contexts you actually need.

For the current cholangiocyte pilot, sample 11 is adult Zhang2021 hepatic
tissue and sample 32 is fetal Domcke2020 hepatic tissue. Sample 32's cCRE
matrix is about 5.9 GB; download it only if an explicit fetal-proxy analysis
is desired. If this host's certificate chain is not trusted by your remote
machine, investigate the cause; `catalog` and `fetch` offer an explicit
`--insecure` fallback, which disables TLS verification and should not be the
normal mode.

## 2. Review cell-type mappings and freeze the RNA candidate set

HSC labels lack stable Cell Ontology IDs in the inspected H5AD, so matching
RNA contexts requires an explicit reviewed crosswalk. Its tab-separated
columns are:

```tsv
cell_type_id	tissue_id	sample_id	hsc_cell_type	mapping_tier	reviewed
CL:0000182	UBERON:0002107	11	Hepatocyte	exact_reviewed	Yes
```

This Hepatocyte row is an **example only**, not a match for the current
cholangiocyte candidates. Verify tissue, cell identity, developmental stage,
and disease state before setting `reviewed=Yes`. Permitted tiers are
`exact_reviewed`, `proxy_fetal`, `proxy_related_celltype`, and
`proxy_disease_context`. Fetal or related-cell proxies are not called exact
confirmation. The pipeline writes `mapping_review.tsv` with candidate-context
suggestions but never auto-approves a fuzzy name match.

```bash
python -m human_scatac_corpus_integration.pipeline plan \
  --candidate-summary /data/ts/aggregate/candidate_summary.tsv.gz \
  --index /data/hsc/celltypes.tsv \
  --crosswalk /data/hsc/reviewed_crosswalk.tsv \
  --meme /data/jaspar2026_core_vertebrates.meme \
  --out /data/hsc/rna_first_run
```

Omit `--crosswalk` to obtain a coverage audit and mapping suggestions first.
If `manifest.json` reports `n_mapped_rows: 0`, skip the costly scan and run
`aggregate` immediately. A new crosswalk or candidate set requires a **new**
plan directory, preserving the previous run's provenance.

## 3. Scan the common cCRE universe once, then analyze each matched cell type

```bash
sbatch --export=ALL,HSC_OUT=/data/hsc/rna_first_run,HSC_DATA_DIR=/data/hsc/samples,HSC_GENOME_FASTA=/data/hg38.fa,HSC_PYTHON=/path/to/env/bin/python human_scatac_corpus_integration/run_scan_array.sh

python -m human_scatac_corpus_integration.pipeline compile \
  --out /data/hsc/rna_first_run

sbatch --export=ALL,HSC_OUT=/data/hsc/rna_first_run,HSC_DATA_DIR=/data/hsc/samples,HSC_PYTHON=/path/to/env/bin/python human_scatac_corpus_integration/run_analyze_array.sh

python -m human_scatac_corpus_integration.pipeline aggregate \
  --out /data/hsc/rna_first_run
```

Wait for each array to complete before the next stage. The scan has 32 cCRE
shards; analysis has 16 sample/cell-type shards. Adjust both the Slurm array
range and corresponding `--num-shards` argument together if changing them.
On a single node, run shard indices sequentially or under your scheduler.
The FIMO site threshold is `p < 10⁻⁴`. Multiple single-TF JASPAR profiles
for one TF are unioned; composite TF-pair profiles are excluded. A cCRE is
called accessible in the cell type if at least 5 cells and at least 5% of its
cells have a nonzero matrix count. Zero-support cCREs form the control pool;
intermediate-support cCREs are excluded. Foreground and controls are matched
by GC decile and cCRE-width quartile. Groups with fewer than 50 cells or fewer
than 1,000 matched cCREs are inconclusive, not negative.

For each TF and the joint same-cCRE motif signature, the analysis performs a
one-sided Fisher test and BH FDR across all tested sample–TF and sample–pair
hypotheses. The default `q ≤ 0.05` screen also requires odds ratios > 1,
at least 20 joint accessible cCREs, and greater joint motif lift in accessible
than control cCREs. Results are explicitly split into
`corpus_coenriched_independent_source`,
`corpus_coenriched_shared_catlas_source`, `proxy_coenriched`, tested negatives,
and coverage/inconclusive states. The main deliverables are:

- `aggregate/corpus_pair_results.tsv.gz`: one RNA candidate × sample × cell-type result;
- `aggregate/candidate_summary.tsv`: pair-level evidence across datasets;
- `aggregate/summary.json`: coverage, counts, thresholds and provenance.

## 4. Assemble the final evidence table without ATAC gating

After the two ATAC branches finish (or report a coverage gap), left-join their
outputs onto **all RNA-replicated candidates**. This step does not rerun or
reselect pairs:

```bash
python -m human_scatac_corpus_integration.pipeline combine \
  --candidate-summary /data/ts/aggregate/candidate_summary.tsv.gz \
  --catlas-results /data/catlas_pairs_run/aggregate/catlas_pair_results.tsv.gz \
  --corpus-summary /data/hsc/rna_first_run/aggregate/candidate_summary.tsv \
  --out /data/tf_pair_rna_atac_evidence.tsv
```

Either ATAC input may be omitted if that branch has not been run; its status
then reads `not_assessed`. The combined TSV keeps every RNA-selected pair and
reports strict CATlas support, proxy CATlas support, corpus support from an
independent source study, shared-source-only support, tested non-support, or
no matched ATAC test in **separate columns**. Its sidecar manifest stores the
input checksums and counts. No CATlas status is consulted when selecting or
analyzing HSC candidates.

Motif co-enrichment is a **prioritization signal**, not simultaneous TF
occupancy, heterodimerization or gene-expression causality. Similar TF-family
motifs and overlapping sequence hits, especially for AP-1 factors, require
additional motif-similarity, footprint/ChIP, or perturbation validation.
The corpus analysis uses its cCRE universe, including any promoter-proximal
elements; it is not numerically interchangeable with the distal-only CATlas
screen.

## Current pilot

The 20 donor-supported Tabula Sapiens pairs are all from adult intrahepatic
cholangiocytes. The indexed adult hepatic corpus sample 11 contains 24 HSC
cell-type labels but no cholangiocyte label; moreover it is the Zhang2021
study used in CATlas. The present corpus pilot therefore has **20 unmapped
RNA candidates, zero ATAC tests**, not 20 negative results. This statement
is for the indexed adult hepatic sample, not a claim that every corpus file
has been exhaustively inspected. The separate fetal hepatic study is a
possible exploratory proxy after its labels are indexed and reviewed, but
would not establish adult-cell replication.
