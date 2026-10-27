# TFpairminer v2 cluster pipeline: RUNBOOK

This pipeline takes every CCLA cell line (39) and every Human scATAC Corpus sample (77, which expand to one
profile per cell type with ≥ `motifs.min_cells` cells) and produces all evidence levels the portal shows:

| level | what | stage |
|---|---|---|
| ATAC motif co-occurrence | pairs / triplets / modules, permutation q, replicate split | 02, 03, 03b |
| ChromBPNet | model (ENCODE pretrained or de-novo), contributions, TF-MoDISco, in-silico pair ablation | 04–06b |
| RNA level 1 | targeted, donor-aware co-expression in CELLxGENE Census (Tabula Sapiens excluded) | 07 |
| RNA level 2 | the same in Tabula Sapiens All Cells | 08 |
| strict (atlas-wide) | uncapped Census/TS screens (`tf_pair_evidence/run_remote.sh`), then `integrate.py` | 07a/07b, 08a/08b, 09 |
| integration | ATAC×RNA classes (01_atac_only … 05_unassessed), `master.tsv`, `site/atac_rna.json` | 09, 10 |
| portal data + QC | `site/` copied atomically to `publish.dir`, `qc/QC.md` | 10, 11 |

Nothing here writes into v1 (`docs/data`, `docs/index.html`). Only stage 10 writes outside `OUT`: it replaces `publish.dir`
(production: `docs/v2/data` in your checkout). It refuses `docs/`, `docs/data` (v1) and any directory that does not end in `data`.
To keep the current portal data untouched until you have reviewed it, point `publish.dir` at a preview dir first.
Committing and pushing is a separate, manual step (see "Publish").

## Files

```
mining/v2_pipeline/
  pipeline.py            driver: run | resume | dry-run | status | stages  (exec/expand are used internally by wr)
  config.yaml            production config (paths through ${ENV:-default})
  config.smoke.yaml      synthetic end-to-end test
  samples.tsv            manifest: 116 rows (39 CCLA + 77 Corpus) with URLs, sizes, md5, pretrained ChromBPNet + sha256
  sources/build_manifest.py   regenerates samples.tsv from CNGB / Tsinghua Corpus / HuggingFace (ENCODE models)
  context_map.tsv        ATAC profile or cell type -> RNA context (CL ids, exact | proxy | unmapped)
  envs/*.yml, containers/tfpm_v2.def   conda envs and an Apptainer image
  smoke/make_smoke.py    tiny synthetic dataset
  tests/test_pipeline.py unit tests
```

Each stage is a set of independent **units** (one per profile, shard or global step). Every unit:
- writes `OUT/state/<stage>/<unit>.{running,done,skipped,failed}` markers and `OUT/logs/<stage>/<unit>.log`;
- only counts as `.done` when its declared outputs exist and are non-empty;
- is skipped on re-run once it is `.done` or `.skipped`, so every command can safely be repeated. Use `--force` to redo.

`.skipped` means "not applicable, with a reason". Examples: no fragments for ChromBPNet training, fewer than
`min_fragments` fragments, TF-MoDISco found no seqlets, or strict screens are disabled. Skipped reasons are listed in `qc/QC.md`.

## Stages

| stage | units | wr resources (config `resources`) | notes |
|---|---|---|---|
| 00 env | 1 | 1 cpu, 2 GB | checks the 3 pythons and the tools (`OUT/env_report.json`) |
| 01 fetch | refs + 1 per sample | 2 cpu, 4 GB, ≤6 h | hg38 (UCSC md5), chrom sizes, inputs (size/md5 checks), ENCODE models (sha256), pinned in `*.lock.tsv` |
| 02 windows | 1 per sample | 8 cpu, 64 GB | CCLA: fragments → replicate-aware hotspots. Corpus: cCRE pseudobulk per cell type, plus per-cell-type fragments split by barcode (for training) |
| 03 motifs | 1 per sample | 8 cpu, 32 GB, ≤12 h | `motif_cooccurrence.cli sample` (scan, 999 permutations, modules) |
| 03b aggregate_atac | 1 | 2 cpu, 16 GB | builds `site/` + `site/index.json`, which 05 needs (candidate list) |
| 04 cbp_model | 1 per profile | GPU, 64 GB, ≤48 h | pretrained ENCODE model if the manifest has one, else de-novo `train.sh` (MACS3 → bias → ChromBPNet, ENCODE fold_0) |
| 05 contrib | `contrib_shards` per profile | GPU, 32 GB | integrated gradients over the candidate windows |
| 05b modisco | 1 per profile | 8 cpu, 32 GB | TF-MoDISco patterns matched to JASPAR |
| 06 ablation | `pairs_shards` per profile | GPU, 32 GB | in-silico pair ablation (synergy) |
| 06b merge | 1 per profile | 1 cpu | `pairs_cbp.tsv`, `cbp_summary.json`, `tracks.json` |
| 07a / 08a | 1 each | 2 cpu | strict screen plans (Census / Tabula Sapiens) |
| 07 / 08 | 1 per mapped profile + `strict.shards` workers | 2 cpu, 16 GB | targeted RNA per profile, plus strict screen workers |
| 07b / 08b | 1 each | 2 cpu | strict aggregate |
| 09 integrate | targeted + strict | 2 cpu, 16 GB | `OUT/atac_rna` (targeted) and `OUT/atac_rna_strict` |
| 10 publish | 1 | 2 cpu | re-aggregates `site/` with the RNA layers, then copies it atomically to `publish.dir` |
| 11 qc | 1 | 1 cpu | `qc/profiles.tsv` and `qc/QC.md` (per profile: windows, pairs, CBP, RNA strata/hits, failures, skips) |

Corpus cell-type profiles (`<sample>__<CellType>`) are only known after stage 02. The `wr` runner therefore
expands each stage just before it runs, so the full unit counts only appear in `status` / `dry-run` after 02.

Runtime estimates for the full set (not yet measured on the cluster):
- 02/03: about 1–4 h per sample.
- 04: pretrained models take minutes. De-novo training takes about 6–24 h per profile on one GPU, so this is the long pole.
- 05/06: about 1–3 h per shard on a GPU.
- 07/08 targeted: minutes per profile. Strict: about 200 workers at about 1–2 h each.

## 1. Set up (once)

```bash
git clone https://github.com/double2strand/TFpairminer.git && cd TFpairminer
conda env create -f mining/v2_pipeline/envs/main.yml
conda env create -f mining/v2_pipeline/envs/rna.yml
conda env create -f mining/v2_pipeline/envs/chrombpnet_gpu.yml          # optional fallback; the container (1b) is recommended
#  or: apptainer build tfpm_v2.sif mining/v2_pipeline/containers/tfpm_v2.def
export CONDA_PREFIX_MAIN=$(conda info --base)/envs/tfpm-v2
export CONDA_PREFIX_RNA=$(conda info --base)/envs/tfpm-v2-rna
export CONDA_PREFIX_GPU=$(conda info --base)/envs/tfpm-v2-chrombpnet
export TFPM_V2P_OUT=/lustre/scratch/$USER/tfpm_v2p                       # needs ~3-5 TB (fragments + genome + models)
export TFPM_V2P_PY=$CONDA_PREFIX_MAIN/bin/python
cp mining/v2_pipeline/config.yaml mining/v2_pipeline/config.site.yaml     # edit wr queues / GPU misc / resources
```

Then check that everything works before running real data.

## 1b. Sanger farm22 (LSF + wr v0.39): the tested path

Base: `/lustre/scratch126/gengen/projects_v2/ncap_selex/tfpm`, with conda prefixes `envs/main` and `envs/rna`.
ChromBPNet runs from the official image rather than `envs/chrombpnet_gpu.yml`, because chrombpnet pins
tensorflow 2.8.0 / numpy 1.23.4 / python 3.8.

```bash
B=/lustre/scratch126/gengen/projects_v2/ncap_selex/tfpm
module load ISG/apptainer/1.4.0
apptainer pull $B/containers/chrombpnet.sif docker://kundajelab/chrombpnet:latest
# GPU check (TF 2.8.0 should list GPU:0):
bsub -q gpu -gpu "num=1" -M 16000 -R "select[mem>16000] rusage[mem=16000]" -Is \
  apptainer exec --nv $B/containers/chrombpnet.sif python -c "import tensorflow as tf; print(tf.__version__, tf.config.list_physical_devices('GPU'))"
cd $B/TFpairminer && git pull
export TFPM_V2P_PY=$B/envs/main/bin/python          # runs the driver / expanders
C=mining/v2_pipeline/config.farm22.yaml               # out defaults to $B/v2p_out (override with TFPM_V2P_OUT)
wr manager start -s lsf
./tfpm v2 pipeline --config $C dry-run --stage all --runner wr | less
./tfpm v2 pipeline --config $C run --stage 00-03b --runner wr        # ATAC first; curate context_map.tsv after 02
./tfpm v2 pipeline --config $C run --stage 04-11 --runner wr
```

What `config.farm22.yaml` sets up:
- **ChromBPNet container:** `chrombpnet.container` makes stages 04/05/06 run
  `module load ISG/apptainer/1.4.0; apptainer exec --nv -B <repo>,<out>,/lustre/scratch126 --env PYTHONPATH=<repo>/mining <sif> python -m motif_cooccurrence.chrombpnet …`.
  Stage 05b (MoDISco, CPU) uses the same image without `--nv`. De-novo training (`train.sh`) runs inside the image too,
  using `macs2` if `macs3` is absent and `gzip` if `bgzip` is absent.
- **Queues:**
  - CPU stages default to `normal`; 02, 03, 07 and 08 go to `long`.
  - 05 and 06 go to `gpu`. 04 goes to `gpu-basement`, with a 72 h walltime, because de-novo training is long.
  - Each stage's queue comes from `resources.<stage>.queue`.
  - GPU jobs pass `--misc '-gpu "num=1"'` to wr. wr itself turns `-m` into `-M`/`-R rusage[mem]`, so do not put memory in `misc`.
  - `override: higher` stops wr's learned estimates from going below the configured memory and time.
- **Expanders:** the expander jobs get no explicit queue, so wr picks one from their walltime.
- **Extra binds:** if inputs live outside the repo, `out` and `/lustre/scratch126`, add them to `chrombpnet.container_binds`.

Not yet verified on farm22: whether the repo's TFpairminer code runs end-to-end under the image's python 3.8 / TF 2.8
(it is 3.8-compatible syntax-wise). The first `04` unit is the test: check `OUT/logs/04/<unit>.log`.

## 2. Smoke test (about 3 minutes, no network, no GPU)

```bash
export TFPM_V2P_SMOKE=/tmp/tfpm_v2p_smoke TFPM_PY_MAIN=$CONDA_PREFIX_MAIN/bin/python \
       TFPM_PY_RNA=$CONDA_PREFIX_RNA/bin/python TFPM_PY_GPU=$CONDA_PREFIX_GPU/bin/python
$TFPM_PY_GPU mining/v2_pipeline/smoke/make_smoke.py --out $TFPM_V2P_SMOKE/data     # needs tensorflow for the toy model
./tfpm v2 pipeline --config mining/v2_pipeline/config.smoke.yaml run --stage all --runner local
./tfpm v2 pipeline --config mining/v2_pipeline/config.smoke.yaml status          # expect 0 failed
$TFPM_PY_MAIN -m pytest mining/v2_pipeline/tests -q
```

The smoke data plants HNF4A+FOXA2 (liver-like) and SPI1+CEBPB (macrophage) motif pairs. Stage 03 should recover
them in the matching profiles. The toy ChromBPNet model is random, so CBP synergy is 0 and MoDISco may report "skipped".

## 3. Full run with wr

```bash
wr manager start -s lsf            # or your site's usual flags (see Troubleshooting)
C=mining/v2_pipeline/config.site.yaml
./tfpm v2 pipeline --config $C dry-run --stage all --runner wr | less     # review the commands; nothing is submitted
./tfpm v2 pipeline --config $C run --stage all --runner wr               # submits the expander chain and returns
./tfpm v2 pipeline --config $C status                                    # sample x stage matrix
wr status -i tfpm-v2-04 -o s                                              # wr's own view of one stage
```

The wr runner submits one small **expander** job per stage (`tfpm-v2-exp-<stage>`). Each expander depends (`-d`) on the
previous one. When an expander starts, it adds its stage's units (`wr add -i tfpm-v2-<stage> -e tfpm-v2-<stage>`) with that
stage's cpus, memory, time and retries (GPU stages also get `--queue gpu_queue --misc gpu_misc`). It then waits for their
markers and exits non-zero if any unit failed, which blocks the following stages. A unit killed with SIGKILL (OOM or
walltime) is turned into `.failed` once its `.running` marker is older than time × (retries+1) + 1 h.

You can also run a subset of stages or samples:
```bash
./tfpm v2 pipeline --config $C run --stage 02-03b --runner wr
./tfpm v2 pipeline --config $C run --stage 04 --runner wr --only 'CCLA_HepG2,CCLA_K562'
./tfpm v2 pipeline --config $C run --stage 05 --runner wr --force --only 'CCLA_.*'   # recompute
```

### Workstation mode (one 256 GB machine, optional GPU)
```bash
# config: runner.local {cpus: 32, mem_gb: 240, gpus: 1}; concurrency = min(cpus/stage cpus, mem/stage mem)
nohup ./tfpm v2 pipeline --config $C run --stage all --runner local > pipeline.out 2>&1 &
./tfpm v2 pipeline --config $C resume --stage all --runner local    # after a crash/reboot: skips finished units
```
Without a GPU, use `chrombpnet.pretrained: manifest` and expect the de-novo profiles to be very slow, or limit them with `--only`.

## 4. Checkpoints (look at these before moving on)

1. **After 01:** `OUT/inputs/*/inputs.lock.tsv` (sha256 of every input) and `OUT/ref/refs.lock.tsv`. A size or md5 mismatch moves the file to `.badsize`/`.badmd5`; re-run to download it again.
2. **After 02:**
   - `OUT/samples/<SID>/celltypes/corpus_windows.json` lists the cell types kept.
   - `fragment_split.json` gives the barcode match rate. Below 30%, training on that sample is skipped.
   - Run `status` to see the newly discovered `<SID>__<CellType>` profiles.
   - Run `./tfpm v2 pipeline --config $C stages`, then check the unmapped labels in `qc/QC.md` after 11, or with `grep unmapped`.
   - **Curate `context_map.tsv` now**: the RNA levels are only computed for mapped profiles.
3. **After 03b:** `OUT/site/index.json` should list every profile. Spot-check `samples/*/pairs.tsv`.
4. **After 04:** each `samples/.../chrombpnet/model/SOURCE.json` records pretrained (ENCODE URL and sha256) vs de-novo. For de-novo models, check the training logs and bigwigs under `chrombpnet/train/`.
5. **After 06b:** `cbp_summary.json` for each profile (tested / synergistic).
6. **After 07–09:**
   - `OUT/rna/<screen>/<profile>/manifest.json` and `strata/*.json`.
   - `OUT/atac_rna/master.tsv` and `summary.json` (labelled TARGETED).
   - `OUT/atac_rna_strict/` once the strict screens finish.
7. **After 10/11:** `qc/QC.md`, then publish.

## 5. Publish to the portal

Stage 10 writes `publish.dir`. The production value is `docs/v2/data`, the v2 portal's data dir; v1 lives in `docs/data` and is never touched.
To publish into the repository:
```bash
./tfpm v2 pipeline --config $C run --stage 10-11 --runner local
git status docs/v2/data            # review: only docs/v2/data/** should change
git add docs/v2/data && git commit -m "v2 data: <date> full CCLA+Corpus run (pipeline <commit>)"
git push origin HEAD:main          # GitHub Pages redeploys
```
For a dry publish, set `publish.dir: /lustre/.../portal_preview/data` and serve that directory with the `docs/v2` HTML.

## 6. Troubleshooting

- **A unit failed:**
  - Read `OUT/logs/<stage>/<unit>.log` (the full command and traceback).
  - Fix the problem, then re-run the same `run` command. Only failed or pending units are redone.
  - wr side: `wr status -i tfpm-v2-<stage> -b` shows buried jobs. `wr kick -i tfpm-v2-<stage>` retries them, or `wr remove -i …` lets the expander re-add them on the next `run`.
- **An expander exited non-zero:** later stages stay blocked by `-d`. Fix the failed units, then
  `wr remove -i tfpm-v2-exp-<stage>` (and the later expanders) and `run --stage <stage>-11 --runner wr` again.
- **wr manager slow, unresponsive or database bloated** (hundreds of thousands of completed jobs in the history):
  1. Check that nothing is running: `wr status -o s`.
  2. Run `wr manager stop`, then `wr manager compact`, then `wr manager start` (the fix that worked on this cluster).
  3. Only if compaction fails: stop the manager, move `~/.wr_production/db` and `db_bk` aside as a backup, and start it again.
  4. Check with `wr status`.
  5. Re-run `./tfpm v2 pipeline … run --stage all --runner wr`.
  
  This is safe because the pipeline's state lives in `OUT/state`, not in wr. Finished units are skipped and only the
  remaining ones are resubmitted. Prefer this compaction over editing the database. Keep `poll_seconds` ≥ 60 to limit
  the load on the manager.
- **Stale `.running` marker** (local run killed): `resume` clears it. With wr, the expander's timeout logic or `--force` handles it.
- **Corpus fragments missing or not matching the barcodes:** `fragments.MISSING`, or `fragment_split.json` below 30%. ChromBPNet de-novo training is
  skipped for that sample. All ATAC motif levels still run.
- **Census access:** stages 07/08 need outbound HTTPS to the CELLxGENE Census S3 bucket (`census_version` is pinned).
  Offline nodes: point `rna.local_h5ad` at a CELLxGENE-schema h5ad export. This disables the strict screens.
- **GPU OOM in 05/06:** lower `--batch` through `contrib_extra`/`pairs_extra`, or raise `gmem` in `gpu_misc`.

## Verified vs not yet verified (as of the commit that added this file)
- Verified on the build box with `config.smoke.yaml`:
  - all stages 00–11 run with the local runner;
  - the same run works through a fake `wr` shim (expander chain, `-i/-e/-d`, per-stage resources);
  - re-runs are idempotent;
  - the production `dry-run --runner wr` is checked;
  - unit tests pass.
- Not yet verified:
  - a real wr manager;
  - GPU de-novo ChromBPNet training (`train.sh`: needs macs3, bedtools and chrombpnet in the GPU env);
  - full-size data;
  - strict atlas screens through the pipeline (`run_remote.sh` itself was used before);
  - the Apptainer build.
