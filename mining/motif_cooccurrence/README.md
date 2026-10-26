# TFpairminer v2 — scATAC motif co-occurrence mining

This toolkit mines **TF motif pairs, triplet sets and modules ("mixtures")** that co-occur inside the same
accessible chromatin windows. It is an analysis layer that is **independent of, and parallel to,** the RNA (Census)
layers of TFpairminer v1. Results are published to the v2 site at <https://double2strand.github.io/TFpairminer/v2/>.

> Every output is a **motif-sequence co-occurrence hypothesis**. None of it comes from ChromBPNet contribution
> scores, TF occupancy (ChIP/CUT&RUN) or evidence of physical interaction.

## Data sources (`manifest.tsv`)

| kind | dataset | what is used |
|---|---|---|
| `fragments` | CCLA / scatlas HCL cancer cell lines (CNGB SCDS0000040) | Two scATAC fragment files per line (`{line}_1`, `{line}_2`), used as replicates |
| `corpus` | Human scATAC Corpus (Tsinghua) | The cell × cCRE `.h5ad`. Each cell type with ≥ `--min-cells` cells is pseudobulked, and two random half-cell pseudoreplicates are formed |

Manifest columns: `sample_id  dataset  context  kind  inputs(comma-separated paths or URLs)  celltype_col`.
URLs are downloaded on first use into `$OUT/inputs/`. Add rows to extend the run; each row is one shard/job.

## Method, per sample (or per Corpus cell type)

1. **Windows.**
   - *Fragments:* Tn5 insertions (+4/−5) go into 50 bp bins. Smoothed hotspots are called, giving 300 bp windows centred on the summits. Blacklisted windows and windows with more than 5% N are dropped.
   - Three window sets are built: `pooled` (top 20k, requiring ≥ 2 cuts in each replicate), `rep1` and `rep2` (top 10k each, called on one replicate only).
   - *Corpus:* the top cCREs by pseudobulk CPM are re-centred to 300 bp. The pseudoreplicate halves give `rep1`/`rep2`.
2. **Motif scan.** The full JASPAR 2024 CORE vertebrates non-redundant set (879 PWMs, bundled in `resources/`) is scanned on both strands with a numba scanner. Each PWM's threshold is the exact score for p = 1e-4 (dynamic programming over a uniform background).
3. **Pairs.**
   - A stratified hypergeometric z-test is run on window co-occurrence for every motif pair with ≥ 20 windows.
   - The background is matched on GC quintile × accessibility tertile × motif-density tertile × repeat fraction (≥ 50% soft-masked in UCSC hg38). Motif density is the number of distinct motifs hit; it is used only when ≥ 50 motifs are tested.
   - BH correction gives `q`. A second, pilot-style background without motif-density matching gives `lift_basic` and `q_basic`.
   - The 300 strongest pairs also get stratum-preserving permutations (999).
4. **Replicate split.** The same test is repeated in `rep1` and `rep2` windows. `replicate_support = yes` when lift > 1 and p < 0.05 in both.
5. **Family ambiguity.** A pair is flagged for any of:
   - a shared TF family (regex families);
   - PWM similarity ≥ 0.75 (best ungapped column-correlation alignment);
   - ≥ 50% overlapping sites;
   - a heterodimer motif (`::`);
   - joint windows that are ≥ 50% repeat-rich.
6. **Status.**
   - `priority_distinct` requires all of: q ≤ 0.05, matched lift ≥ 1.05, basic lift ≥ 1.25, n_both ≥ 10, ≥ 5 windows with non-overlapping sites ≤ 150 bp apart, and no family / PWM / overlap / repeat flag.
   - `family_level` passes the same thresholds but carries one of those flags.
   - Everything else is `exploratory`.
7. **Higher order.**
   - *Triplets:* built from one representative motif per family/PWM cluster, taken from replicate-supported distinct pairs (up to 40 motifs). Each triplet is tested against the best constituent pair (Poisson with stratified expectation, then BH). `synergistic_set` means q ≤ 0.05 and synergy lift ≥ 1.2.
   - *Modules:* maximal cliques (size ≥ 3) of replicate-supported `priority_distinct` pairs.
8. **Aggregation.**
   - Pairs are combined across all samples and cell types.
   - Ratings:
     - *strong:* distinct and replicate-supported in ≥ 3 profiles from ≥ 2 parent samples;
     - *moderate:* the same in ≥ 2 profiles, or in 1 profile with permutation p ≤ 0.01;
     - *preliminary:* the same in 1 profile;
     - *family-level only*;
     - *weak*.
   - The site shows a cross-sample heatmap and site examples.

The original HepG2 chr22 ChromBPNet-free pilot is in `/workspace/tfpm_v2/src/chrombpnet_pair_discovery`. It is not in the repo, and its numbers are quoted in the v2 site's About tab.

## Workstation (local runner)

```bash
git clone https://github.com/double2strand/TFpairminer && cd TFpairminer
./tfpm v2 setup                       # .venv_v2 + requirements + tests (or: conda env create -f mining/motif_cooccurrence/environment.yml)
./tfpm v2 discover                    # downloads UCSC hg38 (~1 GB), builds code arrays (~3 GB), checks all manifest URLs
./tfpm v2 run --runner local --parallel 1          # one sample at a time; each uses all cores via numba
./tfpm v2 status
./tfpm v2 publish                     # aggregate -> docs/v2/data -> git commit + push
```
Resources: one fragment sample takes about 4 min on 8 cores, peaks at about 6 GB RAM and needs about 0.5 GB of disk for outputs. A Corpus h5ad needs more RAM because it is loaded whole: about 4 GB for 6.5k cells × 1.36M cCREs.

## Cluster (wr; no Slurm)

```bash
wr manager start                      # once, on the head node
export TFPM_V2_OUT=/lustre/.../tfpm_v2_run TFPM_WR_QUEUE=normal   # optional queue
./tfpm v2 setup && ./tfpm v2 discover
./tfpm v2 run --runner wr             # wr add -f jobs.txt -i tfpm-v2 -e tfpm_v2_samples --cpus 8 -m 16G -t 12h -r 2
                                      # + aggregate job: wr add -i tfpm-v2-aggregate -d tfpm_v2_samples
wr status -i tfpm-v2                  # or: ./tfpm v2 status
./tfpm v2 publish
```
Extra flags go to every job after `--`, for example `./tfpm v2 run --runner wr -- --max-windows 40000 --min-cells 200`. Use `--only CCLA_K562,CCLA_A549` to submit a subset. Each manifest row is one shard. For finer sharding, split a sample with `--chroms chr1,chr2,...`, giving one row or job per chromosome group.

## Direct CLI

```bash
cd mining
python -m motif_cooccurrence.cli prep-genome --fasta hg38.fa.gz --out ref/hg38_codes
python -m motif_cooccurrence.cli sample --manifest manifest.tsv --sample CCLA_HepG2 --out run --genome ref/hg38_codes \
   --meme motif_cooccurrence/resources/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
   --blacklist motif_cooccurrence/resources/hg38-blacklist.v2.bed.gz [--force | --redo-windows]
python -m motif_cooccurrence.cli aggregate --out run --site run/site --meme ...
```
Per-sample outputs are `windows.tsv`, `hits.npz`, `pairs.tsv`, `pair_site_examples.tsv`, `motif_sets.tsv`, `modules.tsv` and `summary.json`.

## Tests

`cd mining && python -m pytest -q motif_cooccurrence/tests`

The tests cover an exact-threshold check against brute force, BH, strata, and a synthetic genome with a planted motif pair recovered end-to-end through aggregation.

## Status / untested

- The `wr` runner path was **not executed** on a wr cluster from the build box (no wr manager there). The job lines are the same commands the local runner executes.
- The `discover` URL check and the hg38 download run on the box, but `tfpm v2 setup` was exercised with an existing venv.
- Some CCLA lines (HCT116, MCF7, HeLa, GM12878) return 404 at CNGB and are not in the manifest.
