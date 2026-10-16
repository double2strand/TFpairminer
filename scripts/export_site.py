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
         "replication", "atac_level", "atac_corpus"]  # last two: future layers, null = not yet run
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
    p.add_argument("--replication", default=None, help="Tabula Sapiens validation dir (default: <run>/replication_ts if aggregated)")
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

    # ---- replication layer (Tabula Sapiens validation, mining/tabula_sapiens_validate.py)
    rep_dir = Path(a.replication) if a.replication else run / "replication_ts"
    rep_best, rep_assays, rep_info = {}, defaultdict(list), None
    if (rep_dir / "aggregate" / "candidate_summary.tsv.gz").exists():
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "mining"))
        import tabula_sapiens_validate as tsv
        rep_info = {"manifest": json.loads((rep_dir / "manifest.json").read_text()),
                    "summary": json.loads((rep_dir / "aggregate" / "summary.json").read_text())}
        def slim(r):
            return {"status": r["validation_status"], "assay": r["ts_assay"], "n_cells": num(r["ts_n_cells"]), "n_donors": num(r["ts_n_donors"]),
                    "r": num(r["ts_cell_r"]), "lift": num(r["ts_lift"]), "codetected": num(r["ts_codetected_cells"]),
                    "tf1_det": num(r["ts_tf1_detection"]), "tf2_det": num(r["ts_tf2_detection"]), "donor_mean_r": num(r["ts_donor_mean_r"]),
                    "donors_tested": num(r["ts_donors_tested"]), "pos_donor_frac": num(r["ts_positive_donor_fraction"]), "reason": r["reason"]}
        for r in rows(rep_dir / "aggregate" / "replication_results.tsv.gz"):
            rep_assays[r["candidate_id"]].append(slim(r))
        for r in rows(rep_dir / "aggregate" / "candidate_summary.tsv.gz"):
            rep_best[r["candidate_id"]] = {**slim(r), "n_assays_tested": num(r["n_assays_tested"]), "n_assays_reported": num(r["n_assays_reported"])}
    # ---- scATAC layer (CATlas same-cCRE motif co-enrichment, mining/catlas_pipeline.py) ~ level B
    ATAC_LABEL = {"catlas_coenriched": "co-enriched", "not_coenriched": "not co-enriched",
                  "provisional_coenriched": "co-enriched (provisional mapping)", "provisional_not_coenriched": "not co-enriched (provisional mapping)",
                  "proxy_coenriched": "co-enriched in fetal proxy", "proxy_not_coenriched": "not co-enriched in fetal proxy",
                  "no_exact_celltype": "no matching CATlas cell type", "ontology_present_cluster_unresolved": "CATlas cluster unresolved",
                  "motif_unavailable": "no JASPAR motif", "shared_motif_ambiguous": "shared/ambiguous motif",
                  "inconclusive_sparse_motif": "inconclusive (sparse motif)", "insufficient_matched_cCREs": "too few matched cCREs"}
    def atac_load(d):
        if not (d / "aggregate" / "catlas_pair_results.tsv.gz").exists():
            return None
        res = defaultdict(list)
        for r in rows(d / "aggregate" / "catlas_pair_results.tsv.gz"):
            res[r["candidate_id"]].append({"result": r["catlas_result"], "label": ATAC_LABEL.get(r["catlas_result"], r["catlas_result"]),
                "cluster": r["catlas_cluster"], "stage": r["catlas_stage"], "tier": r["mapping_tier"], "reviewed": r["mapping_reviewed"],
                "mapping_status": r["mapping_status"], "motif_status": r["motif_status"], "tf1_motifs": r["tf1_motifs"], "tf2_motifs": r["tf2_motifs"],
                **{k: num(r[k]) for k in ["n_accessible_distal", "n_matched_foreground", "n_matched_background", "tf1_fg_hits", "tf1_bg_hits",
                   "tf2_fg_hits", "tf2_bg_hits", "joint_fg_hits", "joint_bg_hits", "tf1_odds_ratio", "tf2_odds_ratio", "joint_odds_ratio",
                   "tf1_q", "tf2_q", "joint_q", "joint_lift_fg", "joint_lift_bg"]}})
        man = json.loads((d / "manifest.json").read_text()); summ = json.loads((d / "aggregate" / "summary.json").read_text())
        clusters = list(rows(d / "clusters.tsv")); review = list(rows(d / "mapping_review.tsv"))
        cw = list(rows(d.parent / (d.name + ".inputs") / "crosswalk.tsv")) if (d.parent / (d.name + ".inputs") / "crosswalk.tsv").exists() else []
        fetal = {c["catlas_cluster"] for c in clusters if c["stage"] == "fetal"}
        for rv in review:   # honest label for the coverage gap: a fetal-only counterpart exists
            fo = [n for n in rv["name_suggestions"].split(";") if n in fetal and rv["cell_type"].split()[-1].lower() in n.lower()]
            rv["fetal_only_clusters"] = ";".join(fo)
        return {"res": res, "manifest": man, "summary": summ, "review": review, "crosswalk": cw, "dir": d}
    atac_exact = atac_load(run / "atac_catlas"); atac_proxy = atac_load(run / "atac_catlas_proxy_fetal")
    fetal_only = {rv["cell_type_id"]: rv["fetal_only_clusters"] for rv in (atac_exact or {"review": []})["review"] if rv["fetal_only_clusters"]}
    def atac_obj(cid, ct_id, rep_ok):
        if not (atac_exact or atac_proxy):
            return None
        if cid is None or not rep_ok:
            return {"level": "B", "status": "not_tested", "label": "not tested (not replicated)", "tier": None}
        ex = (atac_exact or {"res": {}})["res"].get(cid, []); px = (atac_proxy or {"res": {}})["res"].get(cid, [])
        o = {"level": "B", "status": "not_tested", "label": "not tested", "tier": None, "exact": ex, "proxy_fetal": px}
        if ex:
            st = ex[0]["result"]
            o.update(status=st, label=ATAC_LABEL.get(st, st), tier="exact" if st in ("catlas_coenriched", "not_coenriched") else ex[0]["tier"] or None)
            if st == "no_exact_celltype" and ct_id in fetal_only:
                o.update(status="no_matching_context", label="no adult CATlas match (fetal only)", fetal_clusters=fetal_only[ct_id])
        if px:
            o["proxy_label"] = ATAC_LABEL.get(px[0]["result"], px[0]["result"]); o["proxy_status"] = px[0]["result"]
        return o
    # ---- Human-scATAC-Corpus layer (mining/human_scatac_corpus_integration): separate from CATlas, RNA-first
    CORPUS_LABEL = {"corpus_coenriched_independent_source": "co-enriched (independent study)",
                    "corpus_coenriched_shared_catlas_source": "co-enriched (CATlas source study)",
                    "not_coenriched": "not co-enriched", "proxy_coenriched": "co-enriched (proxy)", "proxy_not_coenriched": "not co-enriched (proxy)",
                    "no_reviewed_corpus_mapping": "no matched ATAC test", "no_matched_test": "no matched ATAC test"}
    cdir = run / "atac_corpus"; corpus = None
    if (cdir / "aggregate" / "corpus_pair_results.tsv.gz").exists():
        cres = defaultdict(list)
        for r in rows(cdir / "aggregate" / "corpus_pair_results.tsv.gz"):
            cres[r["candidate_id"]].append({"result": r["corpus_result"], "label": CORPUS_LABEL.get(r["corpus_result"], r["corpus_result"].replace("_", " ")),
                "sample_id": r["sample_id"], "sample": r["sample"], "dataset": r["dataset"], "hsc_cell_type": r["hsc_cell_type"],
                "hsc_n_cells": num(r["hsc_n_cells"]), "stage": r["source_stage"], "overlaps_catlas_source": r["overlaps_catlas_source"],
                "tier": r["mapping_tier"], "reviewed": r["mapping_reviewed"], "mapping_status": r["mapping_status"], "motif_status": r["motif_status"],
                "tf1_motifs": r["tf1_motifs"], "tf2_motifs": r["tf2_motifs"],
                **{k: num(r[k]) for k in ["n_accessible_cCREs", "n_matched_foreground", "n_matched_background", "tf1_fg_hits", "tf1_bg_hits", "tf2_fg_hits",
                   "tf2_bg_hits", "joint_fg_hits", "joint_bg_hits", "tf1_odds_ratio", "tf2_odds_ratio", "joint_odds_ratio", "tf1_q", "tf2_q", "joint_q",
                   "joint_lift_fg", "joint_lift_bg"]}})
        csum = {r["candidate_id"]: r for r in rows(cdir / "aggregate" / "candidate_summary.tsv")}
        idx = list(rows(cdir / "celltypes.tsv")) if (cdir / "celltypes.tsv").exists() else []
        recs = list(rows(cdir / "sample_recommendations.tsv")) if (cdir / "sample_recommendations.tsv").exists() else []
        samples = {}
        for r in idx:
            sm = samples.setdefault(r["sample_id"], {"sample_id": r["sample_id"], "sample": r["sample"], "dataset": r["dataset"], "tissue": r["tissue"],
                                     "organ": r["organ"], "pmid": r["pmid"], "stage": r["source_stage"], "overlaps_catlas_source": r["overlaps_catlas_source"],
                                     "n_cCREs": num(r["n_cCREs"]), "labels": []})
            sm["labels"].append([r["hsc_cell_type"], num(r["n_cells"])])
        for sm in samples.values():
            sm["labels"].sort(key=lambda x: -(x[1] or 0))
        corpus = {"res": cres, "sum": csum, "samples": samples, "recs": recs,
                  "manifest": json.loads((cdir / "manifest.json").read_text()), "summary": json.loads((cdir / "aggregate" / "summary.json").read_text()),
                  "review": list(rows(cdir / "mapping_review.tsv")),
                  "combined": json.loads((cdir / "combined_rna_atac_evidence.tsv.manifest.json").read_text()) if (cdir / "combined_rna_atac_evidence.tsv.manifest.json").exists() else None}
    def corpus_obj(cid, r, rep_ok):
        if corpus is None:
            return None
        if cid is None or not rep_ok:
            return {"status": "not_tested", "label": "not tested (not replicated)"}
        res = corpus["res"].get(cid, [])
        if not res:
            return {"status": "not_tested", "label": "not tested"}
        tested = [x for x in res if x["result"] in CORPUS_LABEL and x["result"] not in ("no_reviewed_corpus_mapping",)]
        order = ["corpus_coenriched_independent_source", "corpus_coenriched_shared_catlas_source", "proxy_coenriched", "not_coenriched", "proxy_not_coenriched"]
        best = min(tested, key=lambda x: order.index(x["result"]) if x["result"] in order else 99) if tested else res[0]
        st = best["result"]
        if st == "no_reviewed_corpus_mapping":
            st = "no_matched_test"
            recs = [x for x in corpus["recs"] if x["tissue_id"] == r["tissue_id"]]
            idxd = [corpus["samples"][x["sample_id"]] for x in recs if x["sample_id"] in corpus["samples"]]
            why = (f"Corpus samples for {r['tissue']}: " + (", ".join(f"{x['sample_id']} ({x['dataset']}, {corpus['samples'][x['sample_id']]['stage'] + ' tissue, indexed' if x['sample_id'] in corpus['samples'] else 'not indexed'})" for x in recs) or "none") +
                   ". " + ("; ".join(f"Indexed sample {s['sample_id']} has {len(s['labels'])} cell labels, none reviewed as {r['cell_type']}" for s in idxd) or "No sample indexed yet") +
                   ". No reviewed crosswalk row, so no test was run (not a negative result).")
        else:
            why = ""
        cs = corpus["sum"].get(cid, {})
        return {"status": st, "label": CORPUS_LABEL.get(st, st), "summary": cs.get("corpus_summary"), "reason": why,
                "n_mapped_samples": num(cs.get("n_corpus_mapped_samples", "0")), "n_tested_samples": num(cs.get("n_corpus_tested_samples", "0")),
                "n_independent_supporting_datasets": num(cs.get("n_independent_supporting_datasets", "0")),
                "catlas_source_overlap": any(x["overlaps_catlas_source"] == "Yes" for x in res if x["sample_id"]) or
                                         any(x["overlaps_catlas_source"] == "Yes" for x in corpus["recs"] if x["tissue_id"] == r["tissue_id"]),
                "results": [x for x in res if x["sample_id"]]}
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
        row["replication"] = None; row["atac_level"] = None; row["atac_corpus"] = None
        if rep_info is not None:
            cid = tsv.candidate_id(r)
            if cid in rep_best:
                row["replication"] = {**rep_best[cid], "candidate_id": cid, "assays": rep_assays.get(cid, [])}
            row["atac_level"] = atac_obj(cid, r["cell_type_id"], row["replication"] is not None and row["replication"]["status"] == "replicated_donor_supported")
            row["atac_corpus"] = corpus_obj(cid, r, row["replication"] is not None and row["replication"]["status"] == "replicated_donor_supported")
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
    def corpus_layer():
        if corpus is None:
            return {"status": "not_run", "title": "scATAC: Human-scATAC-Corpus", "column": "atac_corpus"}
        return {"status": "done", "title": "scATAC: Human-scATAC-Corpus (RNA-first, separate from CATlas)", "column": "atac_corpus",
                "reference": "Human-scATAC-Corpus, Nucleic Acids Res 2026;54(D1):D175 (per-sample cell x cCRE matrices)",
                "url": "https://health.tsinghua.edu.cn/human-scatac-corpus/",
                "manifest": corpus["manifest"], "summary": corpus["summary"], "combined": corpus["combined"], "recommendations": corpus["recs"],
                "samples": list(corpus["samples"].values()), "mapping_review": corpus["review"],
                "caveat": ("Same-cCRE motif co-enrichment in accessible vs matched inaccessible corpus cCREs for a REVIEWED matching cell type. "
                           "Not co-binding. Pairs are selected by Tabula Sapiens RNA replication only; CATlas results never select pairs. "
                           "Corpus samples from Zhang2021/Domcke2020 share source studies with CATlas and are not independent. "
                           "'No matched ATAC test' is not a negative result.")}
    def coverage():
        cat_rev = {r["cell_type_id"]: r for r in (atac_exact or atac_proxy or {"review": []})["review"]}
        out_ = []
        for i, st in enumerate(strata):
            st = dict(zip(STRATUM, st))
            rr = [p for p in P if p[0] == i]
            rep_i = PCOLS.index("replication")
            n_rep = sum(1 for p in rr if p[rep_i] and p[rep_i]["status"] == "replicated_donor_supported")
            cr = cat_rev.get(st["cell_type_id"])
            cmap = [m for m in ((atac_proxy or {}).get("crosswalk") or []) if m["cell_type_id"] == st["cell_type_id"]]
            catlas = (None if not (atac_exact or atac_proxy) else
                      {"status": "exact" if cr and cr["exact_ontology_clusters"] else ("fetal_only" if cr and cr.get("fetal_only_clusters") else ("not_checked" if not cr else "none")),
                       "exact": cr["exact_ontology_clusters"] if cr else "", "fetal_only": cr.get("fetal_only_clusters", "") if cr else "",
                       "proxy": ";".join(f"{m['catlas_cluster']} ({m['mapping_tier']})" for m in cmap)})
            cor = None
            if corpus is not None:
                recs = [x for x in corpus["recs"] if x["tissue_id"] == st["tissue_id"]]
                tested_ct = any(x["sample_id"] for res in corpus["res"].values() for x in res if x.get("hsc_cell_type") and x["mapping_status"] == "mapped")
                cor = {"status": "mapped" if tested_ct else ("indexed_no_match" if any(x["sample_id"] in corpus["samples"] for x in recs) else
                                                              ("samples_available" if recs else "no_samples")),
                       "samples": [{"sample_id": x["sample_id"], "dataset": x["dataset"], "indexed": x["sample_id"] in corpus["samples"],
                                    "overlaps_catlas_source": x["overlaps_catlas_source"], "n_cells": num(x["n_cells_catalog"])} for x in recs]}
            out_.append({"s": i, "n_rows": len(rr), "n_replicated": n_rep, "catlas": catlas, "corpus": cor})
        return out_
    def atac_layer():
        L = {"status": "not_run", "title": "scATAC evidence (levels A-D)", "column": "atac_level",
             "levels": {"A": {"name": "cell-type accessibility of both TF loci", "status": "not_run"},
                        "B": {"name": "same-cCRE motif co-enrichment in accessible distal cCREs (CATlas)", "status": "not_run"},
                        "C": {"name": "footprinting / co-accessibility", "status": "not_run"},
                        "D": {"name": "ChIP / perturbation support", "status": "not_run"}}}
        if not (atac_exact or atac_proxy):
            return L
        L["status"] = "done"; L["levels"]["B"]["status"] = "done"
        L["caveat"] = ("Same-cCRE sequence-motif co-enrichment in accessible vs GC/width-matched inaccessible distal CATlas cCREs. "
                       "This is not co-binding, not a TF-TF interaction test and not target-gene assignment; AP-1 family motifs are "
                       "near-identical, so positive AP-1 pairs need motif-similarity filtering, footprinting/ChIP or perturbation. "
                       "Coverage gaps (no matching CATlas cell type) are not negative results.")
        for k, A in (("exact", atac_exact), ("proxy_fetal", atac_proxy)):
            if A:
                L[k] = {"manifest": A["manifest"], "summary": A["summary"], "mapping_review": A["review"], "crosswalk": A["crosswalk"]}
        if atac_proxy:
            L["proxy_fetal"]["note"] = ("Exploratory developmental proxy: adult cell types mapped to a FETAL CATlas cluster by an explicit "
                                        "proxy_fetal crosswalk. Never adult validation.")
        return L
    srcs = [agg / "cell_type_pair_candidates.tsv", agg / "positive_controls.tsv", run / "manifest.json", agg / "summary.json"]
    if rep_info:
        srcs += [rep_dir / "aggregate" / "replication_results.tsv.gz", rep_dir / "aggregate" / "candidate_summary.tsv.gz"]
    atac_srcs = [(A["dir"] / "aggregate" / "catlas_pair_results.tsv.gz", f"atac_{tag}_catlas_pair_results.tsv.gz") for tag, A in (("exact", atac_exact), ("proxy_fetal", atac_proxy)) if A]
    if corpus is not None:
        atac_srcs += [(cdir / "aggregate" / "corpus_pair_results.tsv.gz", "corpus_pair_results.tsv.gz")]
        if (cdir / "combined_rna_atac_evidence.tsv").exists():
            atac_srcs += [(cdir / "combined_rna_atac_evidence.tsv", "combined_rna_atac_evidence.tsv")]
    for src, name in atac_srcs:
        shutil.copyfile(src, out / name); downloads[name] = name
    for src in srcs:
        if not src.exists():
            continue
        if src.suffix == ".gz":
            name = src.name if src.name.startswith("replication_") else "replication_" + src.name
            shutil.copyfile(src, out / name); downloads[name] = name
        elif src.stat().st_size > 40 * 1024 * 1024:
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
                       "replication": {"status": "done" if rep_info else "not_run", "title": "Replication (Tabula Sapiens; SCEA pending)", "column": "replication",
                                       **({"dataset": rep_info["manifest"]["ts_dataset_title"], "dataset_id": rep_info["manifest"]["ts_dataset_id"],
                                           "thresholds": rep_info["manifest"]["thresholds"], "summary": rep_info["summary"],
                                           "independence_note": rep_info["manifest"].get("independence_note", ""),
                                           "complete": rep_info["summary"].get("complete")} if rep_info else {})},
                       "scatac": atac_layer(), "atac_corpus": corpus_layer()},
             "coverage": coverage(),
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
