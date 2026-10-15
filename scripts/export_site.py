#!/usr/bin/env python3
"""Validate a finished tf_census.py run and export it as a versioned, sharded JSON release for the site.

Layout (schema_version 2):
  <site>/runs/index.json                 list of published runs; "latest" is the default on the site
  <site>/runs/<run_id>/run.json          provenance, thresholds, counts, positive controls, file list
  <site>/runs/<run_id>/strata.json       cell contexts (cell type x tissue x disease x assay)
  <site>/runs/<run_id>/tfs.json          TF index (pairs, partners, prior pairs)
  <site>/runs/<run_id>/top_pairs.json    top-N candidate rows (default view)
  <site>/runs/<run_id>/heatmap.json      TF x context detection for the top TFs/contexts
  <site>/runs/<run_id>/tf/<SYMBOL>.json  all candidate rows + profile for one TF
  <site>/runs/<run_id>/ct/<idx>.json     all candidate rows + TF profile for one context
  <site>/runs/<run_id>/*.tsv(.gz)        downloadable source tables
Every number is copied from the aggregate tables; nothing is imputed.
"""
import argparse, csv, datetime, gzip, hashlib, json, re, shutil, subprocess, sys
from collections import defaultdict
from pathlib import Path

csv.field_size_limit(sys.maxsize)
SCHEMA_VERSION = 2
STRATUM = ["cell_type_id", "cell_type", "tissue_id", "tissue", "disease_id", "disease", "assay_id", "assay"]
PCOLS = ["s", "tf1", "tf2", "n_datasets", "n_collections", "n_cells", "co_detected", "r", "lift", "donor_r",
         "primary_stat", "primary_r", "min_q", "k_meta", "meta_r", "meta_p", "meta_q", "tau2", "i2",
         "prior", "prior_evidence", "prior_classes", "pmids", "pair_id", "tf1_det", "tf2_det",
         "replication", "atac_level"]  # last two: future layers, null = not yet run
SRC = {"n_datasets": "n_datasets_with_hit", "n_collections": "n_collections_with_hit", "n_cells": "n_cells_in_hit_datasets",
       "co_detected": "co_detected_cells", "r": "cell_r_weighted_by_cells", "lift": "mean_co_detection_lift",
       "donor_r": "mean_donor_r_when_available", "primary_stat": "primary_stat", "primary_r": "best_primary_r",
       "min_q": "min_q_within_stratum", "k_meta": "n_datasets_in_meta", "meta_r": "meta_r_random", "meta_p": "meta_p_random",
       "meta_q": "meta_q_random", "tau2": "meta_tau2", "i2": "meta_i2", "prior_evidence": "prior_best_evidence",
       "prior_classes": "prior_source_classes", "pmids": "prior_pmids"}
REQUIRED = ["cell_type_id", "tf1_symbol", "tf2_symbol", "n_collections_with_hit", "cell_r_weighted_by_cells", "mean_co_detection_lift"]


def rows(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def num(v, sig=4):
    if v in ("", None):
        return None
    try:
        f = float(v)
    except ValueError:
        return v
    if f.is_integer() and abs(f) < 1e15 and "." not in str(v) and "e" not in str(v).lower():
        return int(f)
    return float(f"{f:.{sig}g}")


def fail(msg):
    sys.exit(f"export_site: VALIDATION FAILED: {msg}")


def validate(run):
    agg = run / "aggregate"
    for f in [run / "manifest.json", run / "features.tsv", run / "datasets.tsv", agg / "summary.json",
              agg / "cell_type_pair_candidates.tsv", agg / "cell_type_profiles.tsv"]:
        if not f.exists():
            fail(f"missing {f}")
    summary = json.loads((agg / "summary.json").read_text())
    if summary.get("n_completed_datasets", 0) < 1:
        fail("no completed datasets")
    with open(agg / "cell_type_pair_candidates.tsv", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    missing = [c for c in REQUIRED if c not in header]
    if missing:
        fail(f"candidate table lacks columns {missing}")
    warnings = []
    if not summary.get("complete"):
        warnings.append(f"run incomplete: {summary['n_completed_datasets']}/{summary['n_planned_datasets']} datasets")
    if "meta_r_random" not in header:
        warnings.append("v1 aggregate (no FDR / meta-analysis columns); those fields will be blank")
    return summary, warnings


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="tf_census.py --out directory (after aggregate)")
    p.add_argument("--site-data", default="docs/data", help="site data root (runs/ is created inside)")
    p.add_argument("--priors", default=None, help="interaction_priors.tsv (adds Pair_ID)")
    p.add_argument("--label", default="")
    p.add_argument("--run-id", default="", help="override the automatic run id")
    p.add_argument("--top-n", type=int, default=3000)
    p.add_argument("--min-profile-detection", type=float, default=0.02)
    p.add_argument("--no-latest", action="store_true", help="publish without making this the default run")
    a = p.parse_args()
    run = Path(a.run).resolve(); agg = run / "aggregate"
    summary, warnings = validate(run)
    for w in warnings:
        print("export_site: warning:", w)
    manifest = json.loads((run / "manifest.json").read_text())
    rc = summary.get("run_config", {})
    run_id = a.run_id or "{}_{}_{}_{}".format(manifest.get("resolved_census_version", "na"), manifest.get("tf_set", "na"),
        "full" if summary.get("atlas_wide") else "pilot",
        hashlib.sha256(json.dumps([manifest.get("register_sha256"), manifest.get("n_datasets"), manifest.get("dataset_selection"), rc], sort_keys=True).encode()).hexdigest()[:8])
    out = Path(a.site_data) / "runs" / run_id
    keep_embed = None
    if out.exists():
        if (out / "embed").exists():
            keep_embed = out.parent / f".{run_id}.embed.tmp"
            shutil.rmtree(keep_embed, ignore_errors=True); shutil.move(str(out / "embed"), keep_embed)
        old_info = json.loads((out / "run.json").read_text()) if (out / "run.json").exists() else {}
        shutil.rmtree(out)
    else:
        old_info = {}
    (out / "tf").mkdir(parents=True); (out / "ct").mkdir()

    feats = {r["ensembl_id"]: r for r in rows(run / "features.tsv")}
    tf_ens = {e for e, r in feats.items() if r["is_tf"] == "True"}
    pair_ids = defaultdict(set)
    if a.priors and Path(a.priors).exists():
        for r in rows(a.priors):
            if r.get("A_HGNC_ID") and r.get("B_HGNC_ID") and r.get("Pair_ID"):
                pair_ids[tuple(sorted((r["A_HGNC_ID"], r["B_HGNC_ID"])))].add(r["Pair_ID"])

    strata, sidx, scells = [], {}, {}
    def stratum(r):
        k = tuple(r[c] for c in STRATUM)
        if k not in sidx:
            sidx[k] = len(strata); strata.append(list(k))
        return sidx[k]
    det = {}
    tf_prof, ct_prof = defaultdict(list), defaultdict(list)
    for r in rows(agg / "cell_type_profiles.tsv"):
        if r["ensembl_id"] not in tf_ens:
            continue
        s = stratum(r); d = float(r["detection_fraction"]); mn = float(r["mean_log1p_cptt"])
        det[(s, r["ensembl_id"])] = d
        scells[s] = int(r["n_cells"])
        if d >= a.min_profile_detection:
            sym = feats[r["ensembl_id"]]["symbol"]
            tf_prof[sym].append([s, round(d, 4), round(mn, 4)])
            ct_prof[s].append([sym, round(d, 4), round(mn, 4)])

    P = []
    for r in rows(agg / "cell_type_pair_candidates.tsv"):
        s = stratum(r)
        h = sorted((feats.get(r["tf1_ensembl"], {}).get("hgnc_id", ""), feats.get(r["tf2_ensembl"], {}).get("hgnc_id", "")))
        row = {"s": s, "tf1": r["tf1_symbol"], "tf2": r["tf2_symbol"]}
        for k, src in SRC.items():
            row[k] = num(r.get(src, "")) if k not in ("primary_stat", "prior_evidence", "prior_classes", "pmids") else r.get(src, "")
        row["prior"] = 1 if r.get("physical_prior_recorded") == "Yes" else 0
        row["pair_id"] = " | ".join(sorted(pair_ids.get(tuple(h), ())))
        d1, d2 = det.get((s, r["tf1_ensembl"])), det.get((s, r["tf2_ensembl"]))
        row["tf1_det"] = None if d1 is None else round(d1, 4); row["tf2_det"] = None if d2 is None else round(d2, 4)
        row["replication"] = None; row["atac_level"] = None
        P.append([row[c] for c in PCOLS])
    # stable rank: collections desc, primary r desc
    ci, pi, ri = PCOLS.index("n_collections"), PCOLS.index("primary_r"), PCOLS.index("r")
    P.sort(key=lambda x: (-(x[ci] or 0), -((x[pi] if x[pi] is not None else x[ri]) or 0)))

    dump = lambda path, obj: Path(path).write_text(json.dumps(obj, separators=(",", ":")))
    by_tf, by_ct = defaultdict(list), defaultdict(list)
    for x in P:
        by_tf[x[1]].append(x); by_tf[x[2]].append(x); by_ct[x[0]].append(x)
    tfs = sorted(set(by_tf) | set(tf_prof))
    for t in tfs:
        dump(out / "tf" / f"{safe(t)}.json", {"tf": t, "cols": PCOLS, "rows": by_tf.get(t, []), "profile_cols": ["s", "detection", "mean_log1p_cptt"], "profile": tf_prof.get(t, [])})
    for s in range(len(strata)):
        dump(out / "ct" / f"{s}.json", {"s": s, "cols": PCOLS, "rows": by_ct.get(s, []), "profile_cols": ["tf", "detection", "mean_log1p_cptt"], "profile": sorted(ct_prof.get(s, []), key=lambda v: -v[1])})
    dump(out / "strata.json", {"cols": STRATUM + ["n_cells", "n_pair_rows"], "rows": [st + [scells.get(i), len(by_ct.get(i, []))] for i, st in enumerate(strata)]})
    pr = PCOLS.index("prior")
    dump(out / "tfs.json", {"cols": ["tf", "file", "n_pair_rows", "n_partners", "n_prior_rows", "n_contexts_detected"],
                            "rows": [[t, safe(t), len(by_tf.get(t, [])), len({x[1] if x[2] == t else x[2] for x in by_tf.get(t, [])}),
                                      sum(x[pr] for x in by_tf.get(t, [])), len(tf_prof.get(t, []))] for t in tfs]})
    dump(out / "top_pairs.json", {"cols": PCOLS, "rows": P[:a.top_n], "total": len(P)})
    deg = defaultdict(int)
    for x in P[:a.top_n]:
        deg[x[1]] += 1; deg[x[2]] += 1
    htf = sorted(deg, key=lambda t: -deg[t])[:40]
    hs = []
    for x in P[:a.top_n]:
        if x[0] not in hs:
            hs.append(x[0])
        if len(hs) >= 60:
            break
    hv = {(t, s): d for t in htf for s, d, _ in tf_prof.get(t, [])}
    dump(out / "heatmap.json", {"tfs": htf, "strata": hs, "detection": [[hv.get((t, s)) for s in hs] for t in htf]})

    downloads = {}
    for src in [agg / "cell_type_pair_candidates.tsv", agg / "positive_controls.tsv", run / "manifest.json", agg / "summary.json"]:
        if not src.exists():
            continue
        if src.stat().st_size > 40 * 1024 * 1024:
            with open(src, "rb") as fi, gzip.open(out / (src.name + ".gz"), "wb") as fo:
                shutil.copyfileobj(fi, fo)
            downloads[src.name] = src.name + ".gz"
        else:
            shutil.copyfile(src, out / src.name); downloads[src.name] = src.name
    controls = list(rows(agg / "positive_controls.tsv")) if (agg / "positive_controls.tsv").exists() else []
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=Path(__file__).parent).stdout.strip()
    except Exception:
        sha = ""
    info = {"schema_version": SCHEMA_VERSION, "run_id": run_id, "label": a.label, "pilot": not summary.get("atlas_wide", False),
            "exported_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "code_commit_at_export": sha, "source_run_dir": str(run), "warnings": warnings,
            "manifest": manifest, "summary": summary, "positive_controls": controls,
            "counts": {"pair_rows": len(P), "unique_pairs": len({(x[1], x[2]) for x in P}), "strata": len(strata),
                       "tfs": len(tfs), "pairs_with_prior": sum(x[pr] for x in P)},
            "layers": {"discovery_cellxgene": {"status": "done", "title": "CELLxGENE Census discovery"},
                       "replication": {"status": "not_run", "title": "Replication (SCEA / Tabula Sapiens)", "column": "replication"},
                       "scatac": {"status": "not_run", "title": "scATAC evidence (levels A-D)", "column": "atac_level", "levels": ["A", "B", "C", "D"]}},
            "profile_detection_floor": a.min_profile_detection, "top_n": a.top_n, "downloads": downloads}
    if keep_embed is not None:
        shutil.move(str(keep_embed), out / "embed")
        if "embedding" in old_info:
            info["embedding"] = old_info["embedding"]
    info.setdefault("embedding", {"status": "not_run", "datasets": []})
    (out / "run.json").write_text(json.dumps(info, indent=1))

    idx_path = Path(a.site_data) / "runs" / "index.json"
    idx = json.loads(idx_path.read_text()) if idx_path.exists() else {"schema_version": SCHEMA_VERSION, "latest": None, "runs": []}
    idx["runs"] = [r for r in idx["runs"] if r["run_id"] != run_id]
    idx["runs"].insert(0, {"run_id": run_id, "label": a.label, "pilot": info["pilot"], "exported_at": info["exported_at"],
                           "census": manifest.get("resolved_census_version"), "n_datasets": summary.get("n_completed_datasets"),
                           "pair_rows": len(P), "strata": len(strata)})
    if not a.no_latest or not idx.get("latest"):
        idx["latest"] = run_id
    idx["schema_version"] = SCHEMA_VERSION
    idx_path.write_text(json.dumps(idx, indent=1))
    (agg / "published_run_id").write_text(run_id + "\n")
    print(json.dumps({"run_id": run_id, **info["counts"]}), "->", out)


if __name__ == "__main__":
    main()
