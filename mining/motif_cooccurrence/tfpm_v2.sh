#!/usr/bin/env bash
# tfpm v2 — scATAC motif co-occurrence mining (pairs, triplet sets, modules). Called as: tfpm v2 <cmd>
# Settings (env or <repo>/tfpm.conf):
#   TFPM_V2_OUT (default <repo>/run_v2)        TFPM_V2_MANIFEST (default mining/motif_cooccurrence/manifest.tsv)
#   TFPM_V2_GENOME_FASTA (default: download UCSC hg38.fa.gz into $OUT/ref)
#   TFPM_V2_MEME / TFPM_V2_BLACKLIST (default: bundled JASPAR2024 CORE non-redundant / ENCODE hg38 blacklist v2)
#   TFPM_V2_WR_MEM (16G) TFPM_V2_WR_TIME (12h) TFPM_V2_WR_CPUS (8) TFPM_WR_QUEUE (optional)
#   NUMBA_NUM_THREADS (defaults to the job's CPUs)
set -euo pipefail
HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
[ -f "$ROOT/tfpm.conf" ] && . "$ROOT/tfpm.conf"
OUT="$(readlink -m "${TFPM_V2_OUT:-$ROOT/run_v2}")"
MANIFEST="$(readlink -m "${TFPM_V2_MANIFEST:-$HERE/manifest.tsv}")"
MEME="${TFPM_V2_MEME:-$HERE/resources/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt}"
BL="${TFPM_V2_BLACKLIST:-$HERE/resources/hg38-blacklist.v2.bed.gz}"
GENOME="$OUT/ref/hg38_codes"
VENV="${TFPM_V2_VENV:-$ROOT/.venv_v2}"; PY="$VENV/bin/python"
die(){ echo "tfpm v2: $*" >&2; exit 1; }
usage(){ cat <<U
Usage: tfpm v2 <command>
  setup                      create $VENV, install requirements, run tests
  discover                   fetch hg38 + build 2-bit-like code arrays; check every manifest input is reachable
  run --runner wr|local [--parallel N] [--only ID[,ID]] [-- extra 'sample' flags, e.g. --max-windows 40000]
                             one job per manifest sample (= shard); wr: -i tfpm-v2 rep group, -e dep group,
                             dependent aggregate job via -d
  status                     per-sample state
  aggregate                  cross-sample consistency -> $OUT/site
  publish [--no-push]        aggregate -> docs/v2/data -> git commit (+ push)
Out: $OUT   Manifest: $MANIFEST
U
}
common(){ echo "--manifest $MANIFEST --out $OUT --genome $GENOME --meme $MEME --blacklist $BL"; }
need_py(){ [ -x "$PY" ] || die "run 'tfpm v2 setup' first"; }

cmd_setup(){
  local base="${TFPM_PYTHON:-$(command -v python3.12 || command -v python3.11 || command -v python3)}"
  [ -x "$PY" ] || "$base" -m venv "$VENV"
  "$PY" -m pip install -q --upgrade pip
  "$PY" -m pip install -q -r "$HERE/requirements.txt"
  (cd "$ROOT/mining" && "$PY" -m pytest -q motif_cooccurrence/tests)
}

cmd_discover(){
  need_py; mkdir -p "$OUT/ref"
  local fa="${TFPM_V2_GENOME_FASTA:-$OUT/ref/hg38.fa.gz}"
  if [ ! -s "$fa" ]; then
    curl -fL --retry 5 -o "$fa.part" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && mv "$fa.part" "$fa"
  fi
  [ -s "$GENOME/chr1.npy" ] || (cd "$ROOT/mining" && "$PY" -m motif_cooccurrence.cli prep-genome --fasta "$fa" --out "$GENOME")
  echo "checking manifest inputs..."
  tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r id ds ctx kind inputs rest; do
    for x in ${inputs//,/ }; do
      if [[ "$x" == *://* ]]; then code=$(curl -sIL -o /dev/null -w '%{http_code}' "$x" || true); else [ -e "$x" ] && code=local || code=MISSING; fi
      printf '%s\t%s\t%s\n' "$id" "$code" "$x"
    done
  done
}

cmd_run(){
  need_py; local runner="" par=1 only="" extra=""
  while [ $# -gt 0 ]; do case "$1" in
    --runner) runner="$2"; shift 2;; --parallel|-P) par="$2"; shift 2;; --only) only="$2"; shift 2;;
    --) shift; extra="$*"; break;; *) die "unknown option $1";; esac; done
  [ -s "$GENOME/chr1.npy" ] || die "genome not prepared; run 'tfpm v2 discover'"
  mkdir -p "$OUT/logs"; local jobs="$OUT/jobs.txt"; : > "$jobs"
  tail -n +2 "$MANIFEST" | cut -f1 | while read -r id; do
    [ -z "$only" ] || [[ ",$only," == *",$id,"* ]] || continue
    echo "cd $ROOT/mining && $PY -m motif_cooccurrence.cli sample $(common) --sample $id $extra > $OUT/logs/$id.log 2>&1" >> "$jobs"
  done
  local n; n=$(wc -l < "$jobs")
  case "$runner" in
    wr)
      command -v wr >/dev/null || die "wr not on PATH"
      wr manager status >/dev/null 2>&1 || die "no wr manager running; start it with: wr manager start"
      local q=(); [ -n "${TFPM_WR_QUEUE:-}" ] && q=(--queue "$TFPM_WR_QUEUE")
      local cpus="${TFPM_V2_WR_CPUS:-8}"
      sed -i "s#^cd #export NUMBA_NUM_THREADS=$cpus; cd #" "$jobs"
      wr add -f "$jobs" -i tfpm-v2 -e tfpm_v2_samples -m "${TFPM_V2_WR_MEM:-16G}" -t "${TFPM_V2_WR_TIME:-12h}" \
         --cpus "$cpus" -r 2 "${q[@]}"
      echo "$HERE/tfpm_v2.sh aggregate" | wr add -i tfpm-v2-aggregate -d tfpm_v2_samples -m 8G -t 2h -r 1 "${q[@]}"
      echo "Queued $n sample jobs + 1 dependent aggregate job (wr status -i tfpm-v2)";;
    local)
      setsid nohup bash -c "xargs -P $par -I{} -d '\n' bash -c '{}' < $jobs; '$HERE/tfpm_v2.sh' aggregate > $OUT/logs/aggregate.log 2>&1 && touch $OUT/logs/local.finished" \
        > "$OUT/logs/local.log" 2>&1 &
      echo "Started $n sample jobs locally, $par at a time (pid $!). Check: tfpm v2 status";;
    *) die "use --runner wr or --runner local";;
  esac
}

cmd_status(){ need_py; (cd "$ROOT/mining" && "$PY" -m motif_cooccurrence.cli status --manifest "$MANIFEST" --out "$OUT"); }
cmd_aggregate(){ need_py; (cd "$ROOT/mining" && "$PY" -m motif_cooccurrence.cli aggregate --out "$OUT" --site "$OUT/site" --meme "$MEME"); }
cmd_publish(){
  local push=1; [ "${1:-}" = "--no-push" ] && push=0
  cmd_aggregate
  rm -rf "$ROOT/docs/v2/data"; mkdir -p "$ROOT/docs/v2"; cp -r "$OUT/site" "$ROOT/docs/v2/data"
  cd "$ROOT" && git add docs/v2/data && git commit -qm "v2: publish motif co-occurrence results ($(date -u +%F))" || echo "nothing to commit"
  [ $push = 1 ] && git push || true
}

case "${1:-}" in
  setup) shift; cmd_setup "$@";; discover) cmd_discover;; run) shift; cmd_run "$@";; status) cmd_status;;
  aggregate) cmd_aggregate;; publish) shift; cmd_publish "$@";; *) usage;;
esac
