#!/usr/bin/env python3
"""TFpairminer v2 cluster pipeline: raw public data -> four evidence levels -> portal JSON.

  ./tfpm v2 pipeline run     --config C [--stage N|all] [--runner local|wr] [--only ID,..] [--force]
  ./tfpm v2 pipeline status  --config C            sample x stage table (done/failed/skipped/pending)
  ./tfpm v2 pipeline resume  --config C [--runner]  clear .failed markers and re-run everything not done
  ./tfpm v2 pipeline dry-run --config C [--stage]   print the exact wr commands (and local unit commands)
  ./tfpm v2 pipeline exec    --config C --stage S --unit U   run ONE unit (what every job calls)
  ./tfpm v2 pipeline expand  --config C --stage S   (wr) submit the units of S and wait for them

Stages (sub-stages run in this order; '--stage 5' = 05 + 05b):
  00   env check               01   fetch + verify inputs, refs, models   02   fragments/cCRE -> windows (+ Corpus pseudobulk)
  03   motif scan + pairs/sets (full permutations)   03b  cross-profile aggregate (site index for 05)
  04   ChromBPNet model per profile (ENCODE pretrained, else bias+chrombpnet training on GPU)
  05   contribution shards     05b  TF-MoDISco             06   ablation synergy shards   06b merge shards
  07a/07/07b  CELLxGENE RNA (strict plan / targeted per profile + strict shards / strict aggregate)
  08a/08/08b  Tabula Sapiens RNA (same, independent replication screen)
  09   integration -> atac_rna (targeted, proxy-labelled) and atac_rna_strict (exact contexts only)
  10   aggregate + publish -> publish.dir                  11   QC report
Every unit writes OUT/state/<stage>/<unit>.{done,failed,skipped} + OUT/logs/<stage>/<unit>.log, validates
its outputs before writing .done, and is skipped when .done exists (idempotent; --force re-runs).
"""
from __future__ import annotations

import argparse, csv, fcntl, socket, uuid, gzip, hashlib, json, os, re, shlex, shutil, subprocess, sys, tarfile, threading, time, traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MINING = ROOT / "mining"
STAGES = ["00", "01", "02", "03", "03b", "04", "05", "05b", "06", "06b", "07a", "07", "07b",
          "08a", "08", "08b", "09", "10", "11"]
STAGE_NAMES = {"00": "env", "01": "fetch", "02": "windows", "03": "motifs", "03b": "aggregate_atac", "04": "cbp_model",
               "05": "contrib", "05b": "modisco", "06": "ablation", "06b": "ablation_merge", "07a": "cxg_strict_plan",
               "07": "cxg_rna", "07b": "cxg_strict_agg", "08a": "ts_strict_plan", "08": "ts_rna", "08b": "ts_strict_agg",
               "09": "integrate", "10": "publish", "11": "qc"}
DATASET = {"CCLA": "CCLA (scatlas HCL cell lines)", "Corpus": "Human scATAC Corpus"}
UA = "Mozilla/5.0 TFpairminer-v2-pipeline"


class Skip(Exception):
    """Raise inside a unit to mark it skipped (with a reason) instead of failed."""


# ----------------------------------------------------------------------------- config
def _expand(s):
    if not isinstance(s, str):
        return s
    def rep(m):
        name, _, default = m.group(1).partition(":-")
        return os.environ.get(name) or default
    s = re.sub(r"\$\{([^}]+)\}", rep, s)
    return os.path.expanduser(os.path.expandvars(s))


def _walk(x):
    if isinstance(x, dict):
        return {k: _walk(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_walk(v) for v in x]
    return _expand(x)


class Cfg:
    def __init__(self, path):
        import yaml
        self.path = Path(path).resolve()
        self.c = _walk(yaml.safe_load(self.path.read_text()))
        self.out = Path(self.c["out"]).resolve()

    def __getitem__(self, k):
        return self.c[k]

    def get(self, *keys, default=None):
        x = self.c
        for k in keys:
            if not isinstance(x, dict) or k not in x:
                return default
            x = x[k]
        return x

    def path_(self, p):
        if not p:
            return None
        p = Path(p)
        return p if p.is_absolute() else (ROOT / p)

    def cbp_container(self):
        c = self.get("chrombpnet", "container")
        return self.path_(c) if c else None

    def gpu_py(self, gpu=True):
        """Python for ChromBPNet stages: the GPU conda env, or `apptainer exec [--nv] <sif> python` when
        chrombpnet.container is set (tensorflow 2.8 / numpy 1.23 / py3.8 inside the kundajelab image)."""
        sif = self.cbp_container()
        if not sif:
            return self.py("gpu")
        C = self["chrombpnet"]
        binds = []
        for b in [str(ROOT), str(self.out)] + list(C.get("container_binds") or []):
            b = os.path.expandvars(str(b))
            if b and b not in binds:
                binds.append(b)
        pre = ""
        if C.get("module_load"):
            pre = ("{ type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; }; "
                   f"{C['module_load']} >/dev/null 2>&1; ")
        nv = "--nv " if gpu else ""
        extra = [os.path.expandvars(str(x)) for x in (C.get("extra_pythonpath") or []) ] if isinstance(C.get("extra_pythonpath"), list) \
            else [os.path.expandvars(str(C["extra_pythonpath"]))] if C.get("extra_pythonpath") else []
        pypath = ":".join(extra + [str(MINING)])      # e.g. pip --target dir with anndata (missing from the image)
        return (f"{pre}{C.get('container_exec', 'apptainer')} exec {nv}-B {shlex.quote(','.join(binds))} "
                f"--env PYTHONPATH={shlex.quote(pypath)} {shlex.quote(str(sif))} {C.get('container_python', 'python')}")

    def py(self, which):
        p = self.get("python", which) or sys.executable
        return p if (os.path.exists(p) or shutil.which(p)) else sys.executable if which == "main" else p

    @property
    def genome_dir(self):
        return self.out / "ref" / "hg38_codes"

    @property
    def fasta(self):
        return self.out / "ref" / "hg38.fa"

    @property
    def meme(self):
        return self.path_(self["refs"]["jaspar_meme"])

    @property
    def blacklist(self):
        return self.path_(self["refs"]["blacklist"])

    def res(self, stage):
        r = dict(cpus=1, mem=4, time="1h", gpu=0); r.update(self.get("resources", stage, default={}) or {}); return r


def read_tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(p, rows, fields):
    # unique tmp per writer: parallel units (threads or wr jobs on other hosts) may write the same file
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(f".{p.name}.{socket.gethostname()}.{os.getpid()}.{threading.get_ident()}.{uuid.uuid4().hex[:8]}.tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    os.replace(tmp, p)


def sha256(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(buf), b""):
            h.update(b)
    return h.hexdigest()


def md5(p, buf=1 << 22):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(buf), b""):
            h.update(b)
    return h.hexdigest()


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


# ----------------------------------------------------------------------------- samples, profiles, contexts
def samples(cfg):
    rows = read_tsv(cfg.path_(cfg["manifest"]))
    only = cfg.get("only") or []
    if only:
        rows = [r for r in rows if any(r["profile_id"] == o or re.fullmatch(o, r["profile_id"]) for o in only)]
    return rows


def sample_dir(cfg, sid):
    return cfg.out / "samples" / sid


def local_inputs(cfg, r, col="inputs"):
    return [str(cfg.out / "inputs" / r["profile_id"] / url_name(u)) for u in r[col].split(",") if u]


def url_name(u):
    return u.rstrip("/").split("/")[-1].split("?")[0]


def profiles(cfg, need_windows=True):
    """ATAC profiles = CCLA lines (1 each) + Corpus sample x cell type (known after stage 02)."""
    out = []
    for r in samples(cfg):
        sd = sample_dir(cfg, r["profile_id"])
        if r["kind"] == "fragments":
            out.append({"profile": r["profile_id"], "sample": r["profile_id"], "dir": sd, "atlas": r["atlas"],
                        "context": r["context"], "celltype": "", "fragments": local_inputs(cfg, r, "fragments"),
                        "pretrained": r.get("cbp_pretrained", ""), "pretrained_sha256": r.get("cbp_pretrained_sha256", ""), "row": r})
            continue
        meta = sd / "celltypes" / "corpus_windows.json"
        if not meta.exists():
            continue
        split = {}
        sj = sd / "celltypes" / "fragment_split.json"
        if sj.exists():
            split = json.loads(sj.read_text()).get("files", {})
        for ct in json.loads(meta.read_text())["cell_types"]:
            out.append({"profile": f"{r['profile_id']}__{ct['slug']}", "sample": r["profile_id"],
                        "dir": sd / "celltypes" / ct["slug"], "atlas": r["atlas"], "context": f"{r['context']} | {ct['cell_type']}",
                        "celltype": ct["cell_type"], "fragments": [split[ct["slug"]]] if ct["slug"] in split else [],
                        "pretrained": "", "pretrained_sha256": "", "row": r})
    return out


def context_rules(cfg):
    p = cfg.path_(cfg.get("context_map"))
    return read_tsv(p) if p and p.exists() else []


def rna_context(cfg, prof):
    """First matching context_map rule. CCLA rules match profile_id; Corpus rules match the cell-type label
    (optionally also a tissue/sample regex). Returns None when unmapped (-> RNA levels 'unassessed')."""
    for rule in context_rules(cfg):
        key = prof["profile"] if rule["match_on"] == "profile" else prof["celltype"]
        if not key or not re.fullmatch(rule["pattern"], key, flags=re.I):
            continue
        if rule.get("sample_pattern") and not re.search(rule["sample_pattern"], prof["sample"], flags=re.I):
            continue
        if rule["mapping_level"] == "unmapped":
            return None
        return rule
    return None


# ----------------------------------------------------------------------------- units
@dataclass
class Unit:
    stage: str
    uid: str
    fn: object
    outputs: list = field(default_factory=list)   # files that must exist + be non-empty for .done
    sample: str = ""


def run(cmd, log, env=None, cwd=None):
    """Run a shell command, tee-ing into the unit log; raise on failure."""
    with open(log, "a") as lf:
        lf.write(f"\n$ {cmd}\n"); lf.flush()
        e = dict(os.environ); e.update(env or {})
        p = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], stdout=lf, stderr=subprocess.STDOUT, env=e, cwd=cwd or MINING)
    if p.returncode:
        raise RuntimeError(f"command failed ({p.returncode}): {cmd[:300]}")


def curl(url, dest, log):
    dest = Path(dest); dest.parent.mkdir(parents=True, exist_ok=True)
    if "://" not in url:
        src = Path(url) if Path(url).is_absolute() else ROOT / url
        if dest.resolve() != src.resolve():
            shutil.copy(src, dest)
        return dest
    tmp = Path(str(dest) + ".part")
    run(f"curl -fL --retry 8 --retry-delay 20 -C - -sS -A {shlex.quote(UA)} -o {shlex.quote(str(tmp))} {shlex.quote(url)}", log)
    os.replace(tmp, dest)
    return dest


def fetch_verified(url, dest, log, size="", md5sum="", sha="", lock=None):
    """Download once; verify size/md5/sha256; pin sha256 in the lock so re-downloads must be identical."""
    dest = Path(dest)
    if not (dest.exists() and dest.stat().st_size > 0):
        curl(url, dest, log)
    got = dest.stat().st_size
    if size and str(size).isdigit() and int(size) != got:
        dest.rename(str(dest) + ".badsize")
        raise RuntimeError(f"{dest.name}: size {got} != expected {size} (moved to .badsize; re-run to re-download)")
    if md5sum and md5(dest) != md5sum:
        raise RuntimeError(f"{dest.name}: md5 mismatch")
    h = sha256(dest)
    if sha and h != sha:
        raise RuntimeError(f"{dest.name}: sha256 {h} != expected {sha}")
    if lock is not None:
        prev = {r["file"]: r for r in read_tsv(lock)} if Path(lock).exists() else {}
        if dest.name in prev and prev[dest.name]["sha256"] != h:
            raise RuntimeError(f"{dest.name}: sha256 differs from pinned lock {lock}")
        prev[dest.name] = {"file": dest.name, "url": url, "bytes": got, "sha256": h}
        write_tsv(lock, list(prev.values()), ["file", "url", "bytes", "sha256"])
    return dest


# ---- stage implementations: each returns Units
def st00(cfg):
    def check(u, log):
        rep = {}
        for k in ("main", "rna", "gpu"):
            if k == "gpu" and cfg.cbp_container():
                sif = cfg.cbp_container()
                rep[k] = {"python": cfg.gpu_py(), "ok": sif.exists(), "err": "" if sif.exists() else f"container missing: {sif}",
                          "note": "container mode: GPU/TF are checked inside the image at stage 04 (not on this node)"}
                continue
            py = cfg.py(k); mods = {"main": "numpy,numba,anndata,scipy,motif_cooccurrence.engine",
                                    "rna": "numpy,anndata,tf_census", "gpu": "tensorflow,modiscolite"}[k]
            p = subprocess.run([py, "-c", f"import {mods}; print('ok')"], cwd=MINING, capture_output=True, text=True) \
                if (os.path.exists(py) or shutil.which(py)) else None
            rep[k] = {"python": py, "ok": bool(p and p.returncode == 0), "err": (p.stderr[-400:] if p else "missing")}
        tools = {t: bool(shutil.which(t)) for t in ("curl", "wr", "macs3", "bedtools", "bgzip", "chrombpnet")}
        gb = Path(cfg.get("python", "gpu_bin") or "/nonexistent")
        tools.update({f"gpu_env:{t}": (gb / t).exists() for t in ("chrombpnet", "macs3", "bedtools", "bgzip")})
        (cfg.out / "env_report.json").write_text(json.dumps({"python": rep, "tools": tools}, indent=1))
        if not rep["main"]["ok"]:
            raise RuntimeError("main python env broken: " + rep["main"]["err"])
    return [Unit("00", "env", check, [cfg.out / "env_report.json"])]


def st01(cfg):
    units = []
    def refs(u, log):
        ref = cfg.out / "ref"; ref.mkdir(parents=True, exist_ok=True); R = cfg["refs"]
        if not (cfg.genome_dir / "chr1.npy").exists() and not (cfg.genome_dir / "DONE").exists():
            gz = ref / url_name(R["genome_fasta_url"])
            want = ""
            if R.get("genome_md5_url"):
                tbl = ref / "md5sum.txt"; curl(R["genome_md5_url"], tbl, log)
                want = next((l.split()[0] for l in tbl.read_text().splitlines() if l.strip().endswith(gz.name)), "")
            fetch_verified(R["genome_fasta_url"], gz, log, md5sum=want, lock=ref / "refs.lock.tsv")
            if not cfg.fasta.exists():
                run(f"{'zcat' if gz.suffix == '.gz' else 'cat'} {gz} > {cfg.fasta}.tmp && mv {cfg.fasta}.tmp {cfg.fasta}", log)
            run(f"{cfg.py('main')} -m motif_cooccurrence.cli prep-genome --fasta {cfg.fasta} --out {cfg.genome_dir}", log)
            (cfg.genome_dir / "DONE").write_text("ok")
        sizes = ref / "hg38.chrom.sizes"
        if R.get("chrom_sizes_url"):
            curl(R["chrom_sizes_url"], ref / "chrom.sizes.all", log)
        else:  # derive from fasta
            run(f"awk '/^>/{{if(n)print n\"\\t\"l; n=substr($1,2); l=0; next}}{{l+=length($0)}}END{{print n\"\\t\"l}}' {cfg.fasta} > {ref}/chrom.sizes.all", log)
        keep = set(R["chroms"])
        sizes.write_text("".join(l + "\n" for l in (ref / "chrom.sizes.all").read_text().splitlines() if l.split("\t")[0] in keep))
        f = R["fold"]; test, valid = f["test"], f["valid"]
        train = [c for c in R["chroms"] if c not in test + valid]
        (ref / "fold_0.json").write_text(json.dumps({"train": train, "valid": valid, "test": test}))
        meta = {k: {"path": str(cfg.path_(R[k])), "sha256": sha256(cfg.path_(R[k]))} for k in ("jaspar_meme", "blacklist")}
        (ref / "refs.json").write_text(json.dumps({**meta, "fasta": str(cfg.fasta), "genome_url": R["genome_fasta_url"]}, indent=1))
    units.append(Unit("01", "refs", refs, [cfg.out / "ref" / "refs.json", cfg.out / "ref" / "hg38.chrom.sizes",
                                           cfg.out / "ref" / "fold_0.json", cfg.fasta]))
    for r in samples(cfg):
        def fetch(u, log, r=r):
            d = cfg.out / "inputs" / r["profile_id"]; lock = d / "inputs.lock.tsv"
            for col, scol in (("inputs", "sizes"), ("fragments", "fragment_sizes")):
                urls = [x for x in r[col].split(",") if x]; sz = (r.get(scol) or "").split(",")
                for i, url in enumerate(urls):
                    try:
                        fetch_verified(url, d / url_name(url), log, size=sz[i] if i < len(sz) else "", lock=lock)
                    except RuntimeError:
                        if col == "fragments" and r["kind"] == "corpus":   # optional (only for ChromBPNet training)
                            (d / "fragments.MISSING").write_text(traceback.format_exc()); continue
                        raise
            if r.get("cbp_pretrained") and cfg.get("chrombpnet", "pretrained") == "manifest" and not cfg.get("chrombpnet", "force_pretrained"):
                fetch_verified(r["cbp_pretrained"], d / "pretrained" / url_name(r["cbp_pretrained"]), log,
                               sha=r.get("cbp_pretrained_sha256", ""), lock=lock)
        units.append(Unit("01", r["profile_id"], fetch, [cfg.out / "inputs" / r["profile_id"] / "inputs.lock.tsv"], r["profile_id"]))
    return units


def run_manifest(cfg):
    """Local-path manifest in the motif_cooccurrence.cli format (written atomically, deterministic)."""
    rows = [{"sample_id": r["profile_id"], "dataset": DATASET.get(r["atlas"], r["atlas"]), "context": r["context"],
             "kind": r["kind"], "inputs": ",".join(local_inputs(cfg, r)), "celltype_col": r.get("celltype_col", "")}
            for r in read_tsv(cfg.path_(cfg["manifest"]))]
    p = cfg.out / "run_manifest.tsv"
    with open(cfg.out / ".run_manifest.lock", "w") as lk:      # serialise writers (flock; Lustre client may make it node-local)
        try:
            fcntl.flock(lk, fcntl.LOCK_EX)
        except OSError:
            pass
        write_tsv(p, rows, ["sample_id", "dataset", "context", "kind", "inputs", "celltype_col"])
    return p


def split_corpus_fragments(cfg, r, log):
    """Corpus: per-cell-type pseudobulk fragment files (for ChromBPNet training) via barcode match to the h5ad."""
    import anndata as ad
    sd = sample_dir(cfg, r["profile_id"]) / "celltypes"; res = {"files": {}, "reason": ""}
    tars = [p for p in local_inputs(cfg, r, "fragments") if Path(p).exists()]
    if not tars:
        res["reason"] = "no fragment file published/fetched"; return res
    a = ad.read_h5ad(local_inputs(cfg, r)[0], backed="r")
    col = r.get("celltype_col") or "Cell_type (HSC)"
    if col not in a.obs:
        res["reason"] = f"no {col} column"; return res
    meta = json.loads((sd / "corpus_windows.json").read_text())
    keep = {ct["cell_type"]: ct["slug"] for ct in meta["cell_types"]}
    bc2ct = {str(b): keep[c] for b, c in zip(a.obs_names, a.obs[col].astype(str)) if c in keep}
    tmp = sd / "_frag_extract"; tmp.mkdir(exist_ok=True)
    srcs = []
    for t in tars:
        if tarfile.is_tarfile(t):
            with tarfile.open(t) as tf:
                for m in tf.getmembers():
                    if m.isfile() and re.search(r"(fragments|frag).*\.(tsv|bed)(\.gz)?$", m.name):
                        tf.extract(m, tmp, filter="data") if sys.version_info >= (3, 12) else tf.extract(m, tmp)
                        srcs.append(tmp / m.name)
        else:
            srcs.append(Path(t))
    norm = lambda b: re.split(r"[#:]", b)[-1]
    alt = {norm(b): c for b, c in bc2ct.items()}
    outs, n_tot, n_hit = {}, 0, 0
    try:
        for s in srcs:
            op = gzip.open if str(s).endswith(".gz") else open
            with op(s, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split("\t", 4); n_tot += 1
                    if len(f) < 4:
                        continue
                    b = f[3].strip(); c = bc2ct.get(b) or alt.get(norm(b))
                    if c is None:   # try sample-prefixed barcodes (sample#barcode / sample_barcode)
                        continue
                    n_hit += 1
                    if c not in outs:
                        outs[c] = gzip.open(sd / f"{c}.fragments.tsv.gz", "wt", compresslevel=3)
                    outs[c].write(line if line.endswith("\n") else line + "\n")
    finally:
        for h in outs.values():
            h.close()
        shutil.rmtree(tmp, ignore_errors=True)
    rate = n_hit / max(n_tot, 1); res["fragment_barcode_match_rate"] = round(rate, 4); res["n_fragments"] = n_tot
    if rate < 0.3:
        for c in outs:
            (sd / f"{c}.fragments.tsv.gz").unlink(missing_ok=True)
        res["reason"] = f"only {rate:.1%} of fragments matched h5ad barcodes; ChromBPNet de-novo training not possible"
        return res
    res["files"] = {c: str(sd / f"{c}.fragments.tsv.gz") for c in outs}
    return res


def st02(cfg):
    units = []
    M = cfg["motifs"]
    for r in samples(cfg):
        def windows(u, log, r=r):
            sys.path.insert(0, str(MINING))
            from motif_cooccurrence import engine as E
            sd = sample_dir(cfg, r["profile_id"]); sd.mkdir(parents=True, exist_ok=True)
            if r["kind"] == "fragments":
                meta = E.windows_from_fragments(local_inputs(cfg, r), str(cfg.genome_dir), str(cfg.blacklist), sd / "windows.tsv",
                                                M["window_bp"], M["max_windows"], M["min_cuts"], M["min_rep_cuts"], None, r["profile_id"])
                (sd / "windows_meta.json").write_text(json.dumps(meta, indent=1, default=str))
            else:
                meta = E.windows_from_corpus(local_inputs(cfg, r)[0], str(cfg.genome_dir), str(cfg.blacklist), sd / "celltypes",
                                             r.get("celltype_col") or "Cell_type (HSC)", M["window_bp"], M["max_windows"],
                                             M["min_cells"], label=r["profile_id"])
                (sd / "windows_meta.json").write_text(json.dumps(meta, indent=1, default=str))
                split = split_corpus_fragments(cfg, r, log)
                (sd / "celltypes" / "fragment_split.json").write_text(json.dumps(split, indent=1))
                if not meta["cell_types"]:
                    raise Skip(f"no cell type with >= {M['min_cells']} cells")
        outs = [sample_dir(cfg, r["profile_id"]) / ("windows.tsv" if r["kind"] == "fragments" else "celltypes/corpus_windows.json")]
        units.append(Unit("02", r["profile_id"], windows, outs, r["profile_id"]))
    return units


def st03(cfg):
    M = cfg["motifs"]; units = []
    for r in samples(cfg):
        def motifs(u, log, r=r):
            man = run_manifest(cfg)
            flags = " ".join(f"--{k.replace('_', '-')} {M[k]}" for k in ("window_bp", "max_windows", "min_cuts", "min_rep_cuts",
                                                                          "min_cells", "motif_pvalue", "permutations", "q", "min_lift"))
            sd = sample_dir(cfg, r["profile_id"])
            (sd / "DONE").unlink(missing_ok=True) if u.get("force") else None
            run(f"{cfg.py('main')} -m motif_cooccurrence.cli sample --manifest {man} --sample {r['profile_id']} --out {cfg.out} "
                f"--genome {cfg.genome_dir} --meme {cfg.meme} --blacklist {cfg.blacklist} {flags}", log)
            for p in profiles(cfg):
                if p["sample"] == r["profile_id"]:
                    for f in ("pairs.tsv", "summary.json"):
                        if not (p["dir"] / f).exists():
                            raise RuntimeError(f"missing {p['dir'] / f}")
        units.append(Unit("03", r["profile_id"], motifs, [sample_dir(cfg, r["profile_id"]) / "DONE"], r["profile_id"]))
    return units


def st03b(cfg):
    def agg(u, log):
        run(f"{cfg.py('main')} -m motif_cooccurrence.cli aggregate --out {cfg.out} --site {cfg.out / 'site'} --meme {cfg.meme}", log)
    return [Unit("03b", "aggregate", agg, [cfg.out / "site" / "index.json"])]


def cbp(prof):
    return prof["dir"] / "chrombpnet"


def link_windows(prof):
    """Corpus cell-type dirs keep windows in ../<slug>.windows.tsv; chrombpnet.py expects <dir>/windows.tsv."""
    w = prof["dir"] / "windows.tsv"
    if not w.exists():
        src = prof["dir"].parent / f"{prof['dir'].name}.windows.tsv"
        if src.exists():
            tmp = w.with_name(f".windows.tsv.{os.getpid()}.{threading.get_ident()}")   # atomic: shards race here
            tmp.symlink_to(os.path.relpath(src, w.parent)); os.replace(tmp, w)


def st04(cfg):
    C = cfg["chrombpnet"]; units = []; G = cfg.gpu_py(); ref = cfg.out / "ref"
    for p in profiles(cfg):
        def model(u, log, p=p):
            o = cbp(p); o.mkdir(parents=True, exist_ok=True)
            base = f"{G} -m motif_cooccurrence.chrombpnet train --out {o}"
            if C.get("force_pretrained"):
                run(f"{base} --pretrained {cfg.path_(C['force_pretrained'])}", log)
                s = json.loads((o / "model" / "SOURCE.json").read_text()); s["source"] = "smoke/forced model: " + str(C["force_pretrained"])
                (o / "model" / "SOURCE.json").write_text(json.dumps(s)); return
            if p["pretrained"] and C.get("pretrained") == "manifest":
                local = cfg.out / "inputs" / p["sample"] / "pretrained" / url_name(p["pretrained"])
                run(f"{base} --pretrained {local}", log)
                s = json.loads((o / "model" / "SOURCE.json").read_text())
                s["source"] = "ENCODE ChromBPNet ATAC " + p["pretrained"].split("/resolve/")[0].split("/")[-1] + " fold_0 bias-corrected"
                s["note"] += " Contributions: integrated gradients on the counts head."; (o / "model" / "SOURCE.json").write_text(json.dumps(s))
                return
            frs = [f for f in p["fragments"] if Path(f).exists()]
            if not frs:
                raise Skip("no pretrained model and no fragments for this profile (Corpus barcode split unavailable)")
            nfrag = sum(1 for f in frs for _ in gzip.open(f, "rt"))
            if nfrag < int(C.get("min_fragments", 0)):
                raise Skip(f"{nfrag} fragments < min_fragments {C['min_fragments']}: too shallow for de-novo ChromBPNet")
            bias = f"--bias-model {C['bias_model']}" if C.get("bias_model") else ""
            env = None if cfg.cbp_container() else {"PATH": f"{cfg.get('python', 'gpu_bin')}:{os.environ['PATH']}"}
            run(f"{base} --name {slug(p['profile'])} --fragments {' '.join(frs)} --fasta {cfg.fasta} --chrom-sizes {ref / 'hg38.chrom.sizes'} "
                f"--blacklist {cfg.blacklist} --fold {ref / 'fold_0.json'} --max-peaks {C['max_peaks']} {bias} --execute", log, env=env)
            (o / "model" / "SOURCE.json").write_text(json.dumps({"source": f"trained de novo on {p['profile']} fragments ({nfrag})",
                "note": "ChromBPNet trained by the v2 pipeline (bias model + chrombpnet pipeline, fold_0).",
                "sha256": sha256(o / "model" / "chrombpnet_nobias.h5")}))
        units.append(Unit("04", p["profile"], model, [cbp(p) / "model" / "chrombpnet_nobias.h5"], p["sample"]))
    return units


def has_model(p):
    return (cbp(p) / "model" / "chrombpnet_nobias.h5").exists()


def st05(cfg):
    C = cfg["chrombpnet"]; n = int(C["contrib_shards"]); units = []
    for p in profiles(cfg):
        for k in range(n):
            def contrib(u, log, p=p, k=k):
                if not has_model(p):
                    raise Skip("no ChromBPNet model (stage 04 skipped)")
                link_windows(p)
                run(f"{cfg.gpu_py()} -m motif_cooccurrence.chrombpnet contrib --out {cbp(p)} --sample-dir {p['dir']} "
                    f"--model {cbp(p) / 'model'} --genome {cfg.genome_dir} --site-index {cfg.out / 'site' / 'index.json'} "
                    f"--shard {k} --nshards {n} --method {C['method']} {C.get('contrib_extra', '')}", log)
            units.append(Unit("05", f"{p['profile']}:{k}", contrib, [cbp(p) / f"contrib.shard{k}.npz"], p["sample"]))
    return units


def st05b(cfg):
    C = cfg["chrombpnet"]; units = []
    for p in profiles(cfg):
        def modisco(u, log, p=p):
            if not has_model(p):
                raise Skip("no ChromBPNet model")
            try:
                run(f"{cfg.gpu_py(gpu=False)} -m motif_cooccurrence.chrombpnet modisco --out {cbp(p)} --meme {cfg.meme} --max-seqlets {C['modisco_max_seqlets']}", log)
            except RuntimeError:
                tail = Path(log).read_text()[-4000:]
                if re.search(r"Sample weights must contain at least one non-zero|No seqlets|zero-size array", tail):
                    raise Skip("TF-MoDISco found no significant seqlets (weak/uninformative contributions)")
                raise
        units.append(Unit("05b", p["profile"], modisco, [cbp(p) / "modisco_patterns.tsv"], p["sample"]))
    return units


def st06(cfg):
    C = cfg["chrombpnet"]; n = int(C["pairs_shards"]); units = []
    for p in profiles(cfg):
        for k in range(n):
            def pairs(u, log, p=p, k=k):
                if not has_model(p):
                    raise Skip("no ChromBPNet model")
                run(f"{cfg.gpu_py()} -m motif_cooccurrence.chrombpnet pairs --out {cbp(p)} --sample-dir {p['dir']} --model {cbp(p) / 'model'} "
                    f"--genome {cfg.genome_dir} --shard {k} --nshards {n} {C.get('pairs_extra', '')}", log)
            units.append(Unit("06", f"{p['profile']}:{k}", pairs, [cbp(p) / f"pairs_cbp.shard{k}.json"], p["sample"]))
    return units


def st06b(cfg):
    units = []
    for p in profiles(cfg):
        def merge(u, log, p=p):
            if not has_model(p):
                raise Skip("no ChromBPNet model")
            run(f"{cfg.py('main')} -m motif_cooccurrence.chrombpnet merge --out {cbp(p)}", log)
        units.append(Unit("06b", p["profile"], merge, [cbp(p) / "pairs_cbp.tsv", cbp(p) / "cbp_summary.json"], p["sample"]))
    return units


def strict_env(cfg, screen=None, shard=None):
    S = cfg["rna"]
    e = {"TFPM_REPO": str(ROOT), "WORK_ROOT": str(cfg.out / "rna_strict"), "REGISTER": str(ROOT / "compendium" / "gene_register.tsv"),
         "CENSUS_VERSION": S["census_version"], "TABULA_DATASET_ID": S["tabula_dataset_id"],
         "PATH": f"{Path(cfg.py('rna')).parent}:{os.environ['PATH']}"}
    if screen:
        e["SCREEN"] = screen
    if shard is not None:
        e.update(SHARD_INDEX=str(shard), NUM_SHARDS=str(S["strict"]["shards"]))
    return e


def strict_on(cfg):
    return bool(cfg.get("rna", "strict", "enabled")) and not cfg.get("rna", "local_h5ad")


def rna_stage(cfg, screen, sub):
    R = cfg["rna"]; units = []; rs = cfg.out / "rna_strict"
    if sub == "plan":
        def plan(u, log):
            if not strict_on(cfg):
                raise Skip("strict atlas screens disabled (rna.strict.enabled=false or local_h5ad set)")
            if not (rs / screen / "plan.json").exists() and not list((rs / screen).glob("plan*")):
                # run_remote.sh plan plans BOTH screens; run it once (07a), 08a only verifies.
                run(f"bash {MINING / 'tf_pair_evidence' / 'run_remote.sh'} plan", log, env=strict_env(cfg))
        return [Unit("07a" if screen == "cellxgene" else "08a", "plan", plan, [])]
    if sub == "agg":
        def agg(u, log):
            if not strict_on(cfg):
                raise Skip("strict atlas screens disabled")
            run(f"bash {MINING / 'tf_pair_evidence' / 'run_remote.sh'} aggregate", log, env=strict_env(cfg, screen))
        return [Unit("07b" if screen == "cellxgene" else "08b", "aggregate", agg, [rs / screen / "aggregate" / "summary.json"])]
    st = "07" if screen == "cellxgene" else "08"
    reg = ROOT / "compendium" / "gene_register.tsv"
    for p in profiles(cfg):
        def targeted(u, log, p=p):
            pc = cbp(p) / "pairs_cbp.tsv"
            if not pc.exists():
                raise Skip("no pairs_cbp.tsv (ChromBPNet levels not available for this profile)")
            rule = rna_context(cfg, p)
            if rule is None:
                raise Skip("profile not mapped to an RNA cell context in context_map.tsv")
            src = f"--local-h5ad {cfg.path_(R['local_h5ad'])}" if R.get("local_h5ad") else f"--census-version {R['census_version']}"
            label = f"{rule['mapping_level']}: {rule.get('note') or p['context']} -> {rule['rna_context']}"
            run(f"{cfg.py('rna')} -m tf_pair_evidence.targeted_rna --pairs-cbp {pc} --register {reg} --context {slug(rule['rna_context'])} "
                f"--cl-ids {rule['cl_ids']} --context-label {shlex.quote(label)} --screen {screen} {src} "
                f"--tabula-dataset-id {R['tabula_dataset_id']} {R['thresholds']} --out {cfg.out / 'rna' / screen / p['profile']}", log)
        units.append(Unit(st, p["profile"], targeted, [cfg.out / "rna" / screen / p["profile"] / "_SUCCESS"], p["sample"]))
    if strict_on(cfg):
        for k in range(int(R["strict"]["shards"])):
            def worker(u, log, k=k):
                run(f"bash {MINING / 'tf_pair_evidence' / 'run_remote.sh'} worker", log, env=strict_env(cfg, screen, k))
            units.append(Unit(st, f"strict:{k}", worker, []))
    return units


def st09(cfg):
    reg = ROOT / "compendium" / "gene_register.tsv"; E = cfg.py("rna")
    def targeted(u, log):
        lines = []
        for p in profiles(cfg):
            pc = cbp(p) / "pairs_cbp.tsv"; r1 = cfg.out / "rna" / "cellxgene" / p["profile"]; r2 = cfg.out / "rna" / "tabula" / p["profile"]
            rule = rna_context(cfg, p)
            if pc.exists() and rule and (r1 / "manifest.json").exists() and (r2 / "manifest.json").exists():
                lines.append(f"--line {shlex.quote(':'.join([p['profile'], slug(rule['rna_context']), str(pc), str(r1), str(r2)]))}")
        if not lines:
            raise Skip("no profile has ChromBPNet pairs + both RNA screens")
        run(f"{E} -m tf_pair_evidence.targeted_integrate {' '.join(lines)} --register {reg} --out {cfg.out / 'atac_rna'}", log)
        s = json.loads((cfg.out / "atac_rna" / "summary.json").read_text())
        s["label"] = "TARGETED proxy/exact-context screen (not atlas-wide); see mapping_level per row"
        (cfg.out / "atac_rna" / "summary.json").write_text(json.dumps(s, indent=1))
    def strict(u, log):
        rs = cfg.out / "rna_strict"; od = cfg.out / "atac_rna_strict"
        if not strict_on(cfg):
            raise Skip("strict screens disabled")
        miss = [s for s in ("cellxgene", "tabula") if not (rs / s / "aggregate" / "summary.json").exists()]
        if miss:
            raise Skip(f"strict atlas screens incomplete: {miss}")
        man, cmap = [], []
        for p in profiles(cfg):
            pc = cbp(p) / "pairs_cbp.tsv"; rule = rna_context(cfg, p)
            if not pc.exists():
                continue
            man.append({"atac_dataset": p["profile"], "atac_context": p["context"], "pairs_cbp_path": str(pc)})
            for cl in (rule["cl_ids"].split(",") if rule else [""]):
                cmap.append({"atac_dataset": p["profile"], "atac_context": p["context"], "cell_type_id": cl,
                             "tissue_id": (rule or {}).get("tissue_id", ""), "disease_id": (rule or {}).get("disease_id", "PATO:0000461"),
                             "mapping_level": (rule or {}).get("mapping_level", "unmapped")})
        od.mkdir(parents=True, exist_ok=True)
        write_tsv(od / "atac_manifest.tsv", man, ["atac_dataset", "atac_context", "pairs_cbp_path"])
        write_tsv(od / "context_map.tsv", cmap, ["atac_dataset", "atac_context", "cell_type_id", "tissue_id", "disease_id", "mapping_level"])
        run(f"{E} -m tf_pair_evidence.integrate --atac-manifest {od / 'atac_manifest.tsv'} --context-map {od / 'context_map.tsv'} "
            f"--register {reg} --rna1 {rs / 'cellxgene'} --rna2 {rs / 'tabula'} --out {od}", log)
    return [Unit("09", "targeted", targeted, [cfg.out / "atac_rna" / "master.tsv"]),
            Unit("09", "strict", strict, [cfg.out / "atac_rna_strict" / "summary.json"])]


def st10(cfg):
    def publish(u, log):
        run(f"{cfg.py('main')} -m motif_cooccurrence.cli aggregate --out {cfg.out} --site {cfg.out / 'site'} --meme {cfg.meme}", log)
        dest = cfg.path_(cfg["publish"]["dir"]).resolve()
        if dest in ((ROOT / "docs").resolve(), (ROOT / "docs" / "data").resolve()) or not str(dest).rstrip("/").endswith("data"):
            raise RuntimeError(f"refusing to publish into {dest} (must be a .../data dir, never docs/ (v1))")
        if (cfg.out / "atac_rna_strict" / "master.tsv").exists():
            shutil.copy(cfg.out / "atac_rna_strict" / "master.tsv", cfg.out / "site" / "atac_rna_strict_master.tsv")
        tmp = dest.with_name(dest.name + ".new"); shutil.rmtree(tmp, ignore_errors=True); shutil.copytree(cfg.out / "site", tmp)
        old = dest.with_name(dest.name + ".old"); shutil.rmtree(old, ignore_errors=True)
        if dest.exists():
            dest.rename(old)
        tmp.rename(dest); shutil.rmtree(old, ignore_errors=True)
        (cfg.out / "published.json").write_text(json.dumps({"dir": str(dest), "time": time.strftime("%F %T"),
                                                             "files": sorted(x.name for x in dest.iterdir())}, indent=1))
    return [Unit("10", "publish", publish, [cfg.out / "published.json"])]


def st11(cfg):
    def qc(u, log):
        rows = []
        for p in profiles(cfg):
            s = json.loads((p["dir"] / "summary.json").read_text()) if (p["dir"] / "summary.json").exists() else {}
            cb = json.loads((cbp(p) / "cbp_summary.json").read_text()) if (cbp(p) / "cbp_summary.json").exists() else {}
            src = json.loads((cbp(p) / "model" / "SOURCE.json").read_text()) if (cbp(p) / "model" / "SOURCE.json").exists() else {}
            rule = rna_context(cfg, p)
            def rn(screen):
                d = cfg.out / "rna" / screen / p["profile"]
                if not (d / "_SUCCESS").exists():
                    return "-"
                st = [json.loads(f.read_text()) for f in (d / "strata").glob("*.json")]
                return f"{len(st)} strata/{sum(x.get('hit', False) for s_ in st for x in s_['pairs'])} hits"
            rows.append({"profile": p["profile"], "atlas": p["atlas"], "n_windows": s.get("n_windows", ""),
                         "n_pairs_tested": s.get("n_pairs_tested", s.get("n_pairs", "")),
                         "n_sig_pairs": s.get("n_sig_pairs", s.get("n_priority_distinct", "")),
                         "cbp_model": (src.get("source") or "")[:60], "cbp_tested": cb.get("n_tested", ""),
                         "cbp_synergistic": cb.get("n_synergistic", ""), "rna_context": rule["rna_context"] if rule else "unmapped",
                         "mapping_level": rule["mapping_level"] if rule else "unmapped", "cellxgene": rn("cellxgene"), "tabula": rn("tabula")})
        F = list(rows[0]) if rows else ["profile"]
        write_tsv(cfg.out / "qc" / "profiles.tsv", rows, F)
        st = status_matrix(cfg)
        fails = [(s, u) for s in STAGES for u, v in st["units"].get(s, {}).items() if v == "failed"]
        skips = []
        for s in STAGES:
            for f in sorted((cfg.out / "state" / s).glob("*.skipped")) if (cfg.out / "state" / s).exists() else []:
                skips.append((s, f.stem, f.read_text().strip()[:160]))
        md = ["# TFpairminer v2 pipeline QC", "", f"config: `{cfg.path}`  out: `{cfg.out}`  generated {time.strftime('%F %T %Z')}", "",
              f"profiles: {len(rows)}; failed units: {len(fails)}; skipped units: {len(skips)}", "",
              "| " + " | ".join(F) + " |", "|" + "---|" * len(F)]
        md += ["| " + " | ".join(str(r[k]) for k in F) + " |" for r in rows]
        md += ["", "## Failed units", ""] + [f"- {s} {u}: see logs/{s}/{u.replace(':', '_')}.log" for s, u in fails]
        md += ["", "## Skipped units (reason)", ""] + [f"- {s} {u}: {why}" for s, u, why in skips]
        sj = cfg.out / "atac_rna" / "summary.json"
        if sj.exists():
            md += ["", "## ATAC+RNA (targeted) counts", "", "```", json.dumps(json.loads(sj.read_text()).get("counts", {}), indent=1), "```"]
        stj = cfg.out / "atac_rna_strict" / "summary.json"
        md += ["", "## Strict atlas integration", "", (stj.read_text()[:2000] if stj.exists() else "not run (see skipped units)")]
        (cfg.out / "qc" / "QC.md").write_text("\n".join(md) + "\n")
    return [Unit("11", "report", qc, [cfg.out / "qc" / "QC.md"])]


BUILDERS = {"00": st00, "01": st01, "02": st02, "03": st03, "03b": st03b, "04": st04, "05": st05, "05b": st05b,
            "06": st06, "06b": st06b, "07a": lambda c: rna_stage(c, "cellxgene", "plan"), "07": lambda c: rna_stage(c, "cellxgene", "units"),
            "07b": lambda c: rna_stage(c, "cellxgene", "agg"), "08a": lambda c: rna_stage(c, "tabula", "plan"),
            "08": lambda c: rna_stage(c, "tabula", "units"), "08b": lambda c: rna_stage(c, "tabula", "agg"),
            "09": st09, "10": st10, "11": st11}
DYNAMIC = {"04", "05", "05b", "06", "06b", "07", "08", "09", "11"}   # unit list depends on stage-02 output


# ----------------------------------------------------------------------------- markers / exec
def mk(cfg, stage, uid, kind):
    return cfg.out / "state" / stage / f"{uid.replace(':', '_').replace('/', '_')}.{kind}"


def unit_state(cfg, stage, uid):
    for k in ("done", "skipped", "failed", "running"):
        if mk(cfg, stage, uid, k).exists():
            return k
    return "pending"


def exec_unit(cfg, u, force=False):
    if not force and unit_state(cfg, u.stage, u.uid) in ("done", "skipped"):
        return unit_state(cfg, u.stage, u.uid)
    d = cfg.out / "state" / u.stage; d.mkdir(parents=True, exist_ok=True)
    logd = cfg.out / "logs" / u.stage; logd.mkdir(parents=True, exist_ok=True)
    log = logd / f"{u.uid.replace(':', '_').replace('/', '_')}.log"
    for k in ("done", "failed", "skipped"):
        mk(cfg, u.stage, u.uid, k).unlink(missing_ok=True)
    mk(cfg, u.stage, u.uid, "running").write_text(f"{os.uname().nodename} {os.getpid()} {time.strftime('%F %T')}")
    t0 = time.time()
    with open(log, "a") as lf:
        lf.write(f"=== {time.strftime('%F %T')} stage {u.stage} ({STAGE_NAMES[u.stage]}) unit {u.uid} on {os.uname().nodename}\n")
    try:
        u.fn({"force": force}, log)
        bad = [str(p) for p in u.outputs if not Path(p).exists() or Path(p).stat().st_size == 0]
        if bad:
            raise RuntimeError("output validation failed, missing/empty: " + ", ".join(bad))
        mk(cfg, u.stage, u.uid, "done").write_text(json.dumps({"seconds": round(time.time() - t0), "outputs": [str(p) for p in u.outputs]}))
        state = "done"
    except Skip as e:
        mk(cfg, u.stage, u.uid, "skipped").write_text(str(e)); state = "skipped"
        with open(log, "a") as lf:
            lf.write(f"SKIPPED: {e}\n")
    except BaseException as e:
        with open(log, "a") as lf:
            lf.write(traceback.format_exc())
        mk(cfg, u.stage, u.uid, "failed").write_text(f"{type(e).__name__}: {e}"[:2000]); state = "failed"
    finally:
        mk(cfg, u.stage, u.uid, "running").unlink(missing_ok=True)
    return state


def stage_units(cfg, stage):
    return BUILDERS[stage](cfg)


def expand_stages(sel):
    if sel in (None, "all"):
        return STAGES
    out = []
    for s in str(sel).split(","):
        if "-" in s:
            a, b = s.split("-"); ia = STAGES.index(norm_stage(a)); ib = max(i for i, x in enumerate(STAGES) if x.startswith(norm_stage(b)))
            out += STAGES[ia:ib + 1]
        else:
            ns = norm_stage(s); out += [x for x in STAGES if x == ns or (x[:2] == ns and len(ns) == 2)]
    return [s for s in STAGES if s in out]


def norm_stage(s):
    s = s.strip().lower()
    return s.zfill(2) if s.isdigit() else (s[:-1].zfill(2) + s[-1] if s[-1].isalpha() else s)


# ----------------------------------------------------------------------------- runners
def run_local(cfg, stages, force=False):
    L = cfg.get("runner", "local"); summary = {}
    for s in stages:
        units = stage_units(cfg, s); r = cfg.res(s)
        par = max(1, min(int(L["cpus"]) // max(1, int(r["cpus"])), int(L["mem_gb"]) // max(1, int(r["mem"]))))
        if r.get("gpu") and int(L.get("gpus", 0)) > 0:
            par = min(par, int(L["gpus"]))
        todo = [u for u in units if force or unit_state(cfg, s, u.uid) not in ("done", "skipped")]
        print(f"[{s} {STAGE_NAMES[s]}] {len(units)} units, {len(todo)} to run, parallel {par}", flush=True)
        os.environ.setdefault("NUMBA_NUM_THREADS", str(max(1, int(r["cpus"]))))
        with ThreadPoolExecutor(par) as ex:
            res = list(ex.map(lambda u: (u.uid, exec_unit(cfg, u, force)), todo))
        states = {uid: unit_state(cfg, s, uid) for uid in (u.uid for u in units)}
        cnt = {k: sum(v == k for v in states.values()) for k in ("done", "skipped", "failed", "pending")}
        summary[s] = cnt; print(f"[{s}] {cnt}", flush=True)
        if cnt["failed"]:
            for uid, v in states.items():
                if v == "failed":
                    print(f"   FAILED {uid}: {mk(cfg, s, uid, 'failed').read_text()[:300]}")
            print(f"stopping: stage {s} has failed units (fix, then `pipeline resume`)"); return 1
    return 0


def me(cfg):
    return f"{shlex.quote(sys.executable)} {shlex.quote(str(Path(__file__).resolve()))} --config {shlex.quote(str(cfg.path))}"


def wr_res(cfg, s, r=None):
    r = r or cfg.res(s); W = cfg.get("runner", "wr")
    a = ["--cpus", str(r["cpus"]), "-m", f"{r['mem']}G", "-t", str(r["time"]), "-r", str(W.get("retries", 2))]
    if W.get("override"):                       # e.g. 'higher': wr must not shrink our mem/time from learned values
        a += ["-o", str(W["override"])]
    q = r.get("queue") or (W.get("gpu_queue") if r.get("gpu") else W.get("queue"))
    if q:
        a += ["--queue", str(q)]
    misc = " ".join(x for x in [W.get("gpu_misc") if r.get("gpu") else "", r.get("misc", "")] if x)
    if misc:                                    # wr adds -M/-R mem itself; misc is only e.g. -gpu "num=1"
        a += ["--misc", misc]
    return a


def wr_unit_line(cfg, s, uid):
    return f"{me(cfg)} exec --stage {s} --unit {shlex.quote(uid)}"


def wr_commands(cfg, stages):
    """Expander chain: job 'expand S' depends on the previous expander; it submits S's units and waits."""
    W = cfg.get("runner", "wr"); P = W["rep_prefix"]; cmds = []; prev = None
    for s in stages:
        dep = ["-d", f"{P}-exp-{prev}"] if prev else []
        line = f"{me(cfg)} expand --stage {s}"
        args = ["wr", "add", "-i", f"{P}-exp-{s}", "-e", f"{P}-exp-{s}", *dep, "--cpus", "1", "-m", "2G", "-t", "72h", "-r", "0"]
        cmds.append((s, line, args)); prev = s
    return cmds


def submit_wr(cfg, stages, dry=False):
    if not dry:
        if not shutil.which("wr"):
            print("wr not on PATH"); return 1
        if subprocess.run(["wr", "manager", "status"], capture_output=True).returncode:
            print("no wr manager running (wr manager start ...; see RUNBOOK)"); return 1
    for s, line, args in wr_commands(cfg, stages):
        print(f"echo {shlex.quote(line)} | {' '.join(shlex.quote(a) for a in args)}")
        if not dry:
            subprocess.run(args, input=line + "\n", text=True, check=True)
    return 0


def _secs(t):
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([smhd]?)", str(t).strip())
    return float(m.group(1)) * {"": 1, "s": 1, "m": 60, "h": 3600, "d": 86400}[m.group(2)] if m else 86400


def expand(cfg, s):
    """Runs as a wr job: add this stage's units (rep group <prefix>-<stage>), then wait for their markers."""
    W = cfg.get("runner", "wr"); P = W["rep_prefix"]
    units = [u for u in stage_units(cfg, s) if unit_state(cfg, s, u.uid) not in ("done", "skipped")]
    if units:
        for u in units:
            mk(cfg, s, u.uid, "failed").unlink(missing_ok=True)
        lines = "\n".join(wr_unit_line(cfg, s, u.uid) for u in units) + "\n"
        args = ["wr", "add", "-i", f"{P}-{s}", "-e", f"{P}-{s}", *wr_res(cfg, s)]
        print("+", " ".join(shlex.quote(a) for a in args), f"<<< {len(units)} units", flush=True)
        subprocess.run(args, input=lines, text=True, check=True)
    while True:
        st = {u.uid: unit_state(cfg, s, u.uid) for u in units}
        limit = _secs(cfg.res(s).get("time", "24h")) * (1 + int(W.get("retries", 0))) + 3600
        for uid, v in st.items():   # a SIGKILLed job (OOM / walltime) leaves a stale .running marker
            r = mk(cfg, s, uid, "running")
            if v == "running" and r.exists() and time.time() - r.stat().st_mtime > limit:
                mk(cfg, s, uid, "failed").write_text("stale .running marker (job killed: OOM/walltime?)\n"); r.unlink(missing_ok=True)
                st[uid] = "failed"
        n = {k: sum(v == k for v in st.values()) for k in ("done", "skipped", "failed", "pending", "running")}
        if n["pending"] + n["running"] == 0:
            break
        if n["failed"]:   # a failed unit may be retried by wr (-r); give up only when wr has buried it
            q = subprocess.run(["wr", "status", "-i", f"{P}-{s}", "-b", "-o", "c"], capture_output=True, text=True)
            if q.returncode == 0 and "buried" in q.stdout.lower() and n["pending"] + n["running"] == 0:
                break
        print(time.strftime("%T"), s, n, flush=True); time.sleep(int(W.get("poll_seconds", 60)))
    bad = [uid for uid, v in st.items() if v == "failed"]
    if bad:
        print(f"stage {s}: {len(bad)} failed units: {bad[:20]}"); return 1
    return 0


# ----------------------------------------------------------------------------- status
def status_matrix(cfg):
    units = {}
    for s in STAGES:
        try:
            us = stage_units(cfg, s)
        except Exception:
            us = []
        units[s] = {u.uid: unit_state(cfg, s, u.uid) for u in us}
        units[s + "_sample"] = {u.uid: u.sample for u in us}
    return {"units": units}


def print_status(cfg):
    st = status_matrix(cfg)["units"]
    sams = [r["profile_id"] for r in samples(cfg)]
    sym = {"done": "✓", "skipped": "–", "failed": "✗", "running": "~", "pending": "."}
    print(f"{'sample':38s} " + " ".join(f"{s:>4s}" for s in STAGES))
    for sid in sams:
        cells = []
        for s in STAGES:
            v = [st[s][u] for u, sm in st[s + "_sample"].items() if sm == sid]
            if not v:
                cells.append("    "); continue
            c = {k: v.count(k) for k in sym}
            cell = "✗" if c["failed"] else "~" if c["running"] else "✓" if c["done"] == len(v) else "–" if c["skipped"] == len(v) \
                else "✓–" if c["done"] + c["skipped"] == len(v) else f"{c['done']}/{len(v)}"
            cells.append(f"{cell:>4s}")
        print(f"{sid[:38]:38s} " + " ".join(cells))
    print(f"{'(global units)':38s} " + " ".join(
        f"{''.join(sym[x] for x in [st[s][u] for u, sm in st[s + '_sample'].items() if not sm])[:4]:>4s}" for s in STAGES))
    tot = {k: sum(v == k for s in STAGES for v in st[s].values()) for k in sym}
    print(f"\nunits: {tot}   legend: ✓ done  – skipped  ✗ failed  ~ running  . pending  n/m partially done")
    print(f"logs: {cfg.out}/logs/<stage>/<unit>.log   markers: {cfg.out}/state/")


# ----------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(HERE / "config.yaml"))
    sp = ap.add_subparsers(dest="cmd", required=True)
    for n in ("run", "resume", "dry-run"):
        p = sp.add_parser(n); p.add_argument("--stage", default="all"); p.add_argument("--runner", choices=["local", "wr"], default="local")
        p.add_argument("--only", default=""); p.add_argument("--force", action="store_true")
    sp.add_parser("status").add_argument("--only", default="")
    e = sp.add_parser("exec"); e.add_argument("--stage", required=True); e.add_argument("--unit", required=True); e.add_argument("--force", action="store_true")
    x = sp.add_parser("expand"); x.add_argument("--stage", required=True)
    sp.add_parser("stages")
    a = ap.parse_args(argv)
    cfg = Cfg(a.config)
    if getattr(a, "only", ""):
        cfg.c["only"] = a.only.split(",")
    cfg.out.mkdir(parents=True, exist_ok=True)
    if a.cmd == "stages":
        for s in STAGES:
            r = cfg.res(s); print(f"{s:4s} {STAGE_NAMES[s]:16s} cpus={r['cpus']} mem={r['mem']}G time={r['time']} gpu={r.get('gpu', 0)}")
        return 0
    if a.cmd == "status":
        print_status(cfg); return 0
    if a.cmd == "exec":
        s = norm_stage(a.stage); u = [u for u in stage_units(cfg, s) if u.uid == a.unit]
        if not u:
            print(f"no unit {a.unit} in stage {s}"); return 2
        st = exec_unit(cfg, u[0], a.force); print(st); return 1 if st == "failed" else 0
    if a.cmd == "expand":
        return expand(cfg, norm_stage(a.stage))
    stages = expand_stages(a.stage)
    if a.cmd == "resume":
        for s in stages:
            for f in (cfg.out / "state" / s).glob("*.failed") if (cfg.out / "state" / s).exists() else []:
                f.unlink()
            for f in (cfg.out / "state" / s).glob("*.running") if (cfg.out / "state" / s).exists() else []:
                f.unlink()
    if a.cmd == "dry-run":
        print(f"# config {cfg.path}\n# out {cfg.out}\n# wr submission (each expander submits its stage's units, then waits):")
        submit_wr(cfg, stages, dry=True)
        print("\n# units each expander will add now (dynamic stages are re-expanded at run time):")
        for s in stages:
            us = stage_units(cfg, s); todo = [u for u in us if unit_state(cfg, s, u.uid) not in ("done", "skipped")]
            args = ["wr", "add", "-i", f"{cfg.get('runner', 'wr')['rep_prefix']}-{s}", "-e", f"{cfg.get('runner', 'wr')['rep_prefix']}-{s}", *wr_res(cfg, s)]
            note = " (+ profiles discovered after stage 02)" if s in DYNAMIC else ""
            print(f"## {s} {STAGE_NAMES[s]}: {len(us)} units, {len(todo)} pending{note}\n#   {' '.join(shlex.quote(x) for x in args)}")
            for u in todo[:3]:
                print("    " + wr_unit_line(cfg, s, u.uid))
            if len(todo) > 3:
                print(f"    ... ({len(todo) - 3} more)")
        return 0
    if a.runner == "wr":
        return submit_wr(cfg, stages)
    return run_local(cfg, stages, a.force)


if __name__ == "__main__":
    sys.exit(main())
