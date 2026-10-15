#!/usr/bin/env python3
"""Validate preselected CELLxGENE TF pairs in the nonduplicated Tabula Sapiens atlas.

The Tabula Sapiens All Cells dataset is discovered from collection metadata in the
same pinned Census release. This is independent-study validation, not an
independent data portal. Discovery overlap is checked before any analysis.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sqlite3
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import tf_census as base


RESULT_FIELDS = [
    "candidate_id", "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
    "cell_type_id", "cell_type", "tissue_id", "tissue", "disease_id", "disease",
    "discovery_assay_id", "discovery_assay", "discovery_r", "discovery_lift",
    "discovery_n_collections", "physical_prior_recorded", "ts_dataset_id",
    "ts_assay_id", "ts_assay", "ts_n_cells", "ts_n_donors",
    "ts_tf1_detection", "ts_tf2_detection", "ts_codetected_cells", "ts_lift",
    "ts_cell_r", "ts_donor_mean_r", "ts_donors_tested", "ts_positive_donor_fraction",
    "validation_status", "reason",
]
CONTEXT_FIELDS = ["cell_type_id", "tissue_id", "disease_id"]
ONTOLOGY_RE = re.compile(r"^[A-Za-z0-9_:.-]+$")


def digest(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:20]


def context_id(row):
    return digest("|".join(row[k] for k in CONTEXT_FIELDS))


def candidate_id(row):
    genes = sorted((row["tf1_ensembl"], row["tf2_ensembl"]))
    return digest("|".join([row["cell_type_id"], row["tissue_id"], row["disease_id"], row["assay_id"], *genes]))


def assert_no_discovery_overlap(discovery_run, ts_dataset_ids):
    discovery_run = Path(discovery_run)
    planned = list(base.read_tsv(discovery_run / "datasets.tsv"))
    overlap = []
    for dataset in planned:
        if dataset["dataset_id"] not in ts_dataset_ids:
            continue
        marker = discovery_run / "datasets" / dataset["dataset_id"] / "_SUCCESS"
        if marker.exists() and json.loads(marker.read_text(encoding="utf-8")).get("n_profile_rows", 0) > 0:
            overlap.append(dataset["dataset_id"])
    if overlap:
        raise RuntimeError("Discovery already used Tabula Sapiens cells in dataset(s) " + ", ".join(overlap) + ". Re-run discovery excluding the Tabula Sapiens collection before claiming independent replication.")
    return len(planned)


def plan(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "manifest.json").exists():
        raise SystemExit("Validation plan already exists; use a new --out directory")
    discovery = Path(args.discovery_run)
    source = discovery / "aggregate" / "cell_type_pair_candidates.tsv"
    if not source.exists():
        raise FileNotFoundError(f"Missing mined candidates: {source}")
    source_manifest = json.loads((discovery / "manifest.json").read_text(encoding="utf-8"))
    source_summary = json.loads((discovery / "aggregate" / "summary.json").read_text(encoding="utf-8"))
    release = source_manifest["resolved_census_version"]
    api = base.census_module()
    with api.open_soma(census_version=release) as census:
        datasets = census["census_info"]["datasets"].read(column_names=["soma_joinid", "dataset_id", "collection_id", "collection_name", "dataset_title", "dataset_total_cell_count"]).concat().to_pandas()
        collection = datasets[datasets["collection_name"] == args.collection_name]
        all_cells = collection[collection["dataset_title"] == args.all_cells_title]
        if len(all_cells) != 1:
            raise RuntimeError(f"Expected exactly one {args.all_cells_title!r} dataset in {args.collection_name!r}; found {len(all_cells)}")
        target = {k: base.clean(v) for k, v in all_cells.iloc[0].to_dict().items()}
        ts_ids = set(collection["dataset_id"])
        assert_no_discovery_overlap(discovery, ts_ids)
        db = out / "candidates.sqlite"
        conn = sqlite3.connect(db)
        conn.execute("CREATE TABLE candidates (candidate_id TEXT PRIMARY KEY, context_id TEXT, cell_type_id TEXT, tissue_id TEXT, disease_id TEXT, tf1_ensembl TEXT, tf2_ensembl TEXT, source_json TEXT)")
        conn.execute("CREATE INDEX candidates_context ON candidates (context_id)")
        needed = set()
        n = 0
        for row in base.read_tsv(source):
            required = ["cell_type_id", "tissue_id", "disease_id", "assay_id", "tf1_ensembl", "tf2_ensembl", "cell_r_weighted_by_cells"]
            if any(not row.get(k) for k in required):
                raise ValueError(f"Candidate row lacks a required field: {row}")
            cid = candidate_id(row)
            xid = context_id(row)
            try:
                conn.execute("INSERT INTO candidates VALUES (?,?,?,?,?,?,?,?)", (cid, xid, row["cell_type_id"], row["tissue_id"], row["disease_id"], row["tf1_ensembl"], row["tf2_ensembl"], json.dumps(row)))
            except sqlite3.IntegrityError as exc:
                raise ValueError(f"Duplicate candidate key: {cid}") from exc
            needed.update((row["tf1_ensembl"], row["tf2_ensembl"]))
            n += 1
        if not n:
            raise ValueError("No candidate pairs supplied")
        conn.commit()
        contexts = [dict(zip(["context_id", *CONTEXT_FIELDS, "n_candidates"], r)) for r in conn.execute("SELECT context_id,cell_type_id,tissue_id,disease_id,COUNT(*) FROM candidates GROUP BY context_id ORDER BY context_id")]
        conn.close()
        base.write_tsv(out / "contexts.tsv", contexts, ["context_id", *CONTEXT_FIELDS, "n_candidates"])
        var = census["census_data"]["homo_sapiens"].ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()
        var = var[var["feature_id"].isin(needed)]
        features = [{"ensembl_id": r.feature_id, "symbol": base.clean(r.feature_name), "soma_joinid": int(r.soma_joinid), "is_tf": True} for r in var.itertuples(index=False)]
        present_joinids = base.get_present_feature_ids(census, int(target["soma_joinid"]), features)
        feature_rows = [{**f, "ts_feature_present": f["soma_joinid"] in present_joinids} for f in features]
        found = {f["ensembl_id"] for f in features}
        feature_rows.extend({"ensembl_id": gid, "symbol": "", "soma_joinid": "", "is_tf": True, "ts_feature_present": False} for gid in sorted(needed - found))
        base.write_tsv(out / "features.tsv", sorted(feature_rows, key=lambda r: r["ensembl_id"]), ["ensembl_id", "symbol", "soma_joinid", "is_tf", "ts_feature_present"])
    source_config = source_summary.get("run_config", {})
    thresholds = {"min_cells": args.min_cells, "min_detection": args.min_detection if args.min_detection is not None else source_config.get("min_detection", 0.05),
                  "min_rho": args.min_rho if args.min_rho is not None else source_config.get("min_rho", 0.25),
                  "min_lift": args.min_lift if args.min_lift is not None else source_config.get("min_lift", 1.1),
                  "min_codetected": args.min_codetected if args.min_codetected is not None else max(25, source_config.get("min_codetected", 25)),
                  "min_donor_cells": args.min_donor_cells, "min_donors": args.min_donors,
                  "min_positive_donor_fraction": args.min_positive_donor_fraction,
                  "assay_ids": sorted(set(args.assay_id))}
    if (thresholds["min_cells"] < 3 or not 0 <= thresholds["min_detection"] <= 1 or
            not -1 <= thresholds["min_rho"] <= 1 or thresholds["min_lift"] < 0 or
            thresholds["min_codetected"] < 1 or thresholds["min_donor_cells"] < 3 or
            thresholds["min_donors"] < 1 or not 0 <= thresholds["min_positive_donor_fraction"] <= 1):
        raise ValueError("Invalid validation thresholds")
    manifest = {"census_version": release, "source_discovery_run": str(discovery.resolve()),
                "source_candidates_sha256": base.sha256(source), "source_atlas_wide": source_summary.get("atlas_wide", False),
                "source_complete": source_summary.get("complete", False),
                "ts_collection_name": args.collection_name, "ts_collection_id": target["collection_id"],
                "ts_dataset_id": target["dataset_id"], "ts_dataset_title": target["dataset_title"],
                "ts_dataset_joinid": int(target["soma_joinid"]), "n_candidates": n,
                "n_contexts": len(contexts), "n_unique_genes": len(needed),
                "n_available_genes": sum(bool(r["ts_feature_present"]) for r in feature_rows),
                "thresholds": thresholds,
                "independence_note": "Biologically distinct Tabula Sapiens study accessed through the same Census portal; source dataset overlap was prohibited"}
    base.write_json_atomic(out / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


def metadata_to_sqlite(census, dataset_id, cell_types, tempdir):
    conn = sqlite3.connect(str(Path(tempdir) / "ts_cells.sqlite"))
    conn.execute("CREATE TABLE cells (jid INTEGER PRIMARY KEY, cell_type_id TEXT, cell_type TEXT, tissue_id TEXT, tissue TEXT, disease_id TEXT, disease TEXT, assay_id TEXT, assay TEXT, donor_id TEXT, raw_sum REAL)")
    allowed = set(cell_types)
    filt = f"is_primary_data == True and dataset_id == '{dataset_id}'"
    if len(allowed) <= 80 and all(ONTOLOGY_RE.fullmatch(x) for x in allowed):
        values = ", ".join(repr(x) for x in sorted(allowed))
        filt += f" and cell_type_ontology_term_id in [{values}]"
    insert = "INSERT INTO cells VALUES (?,?,?,?,?,?,?,?,?,?,?)"
    obs = census["census_data"]["homo_sapiens"].obs
    for table in obs.read(value_filter=filt, column_names=base.OBS_COLUMNS):
        for arrow_batch in table.to_batches(max_chunksize=10000):
            batch = []
            for r in arrow_batch.to_pylist():
                if r["cell_type_ontology_term_id"] not in allowed:
                    continue
                raw_sum = float(r["raw_sum"] or 0)
                if raw_sum <= 0 or not math.isfinite(raw_sum):
                    continue
                batch.append((int(r["soma_joinid"]), base.clean(r["cell_type_ontology_term_id"]), base.clean(r["cell_type"]),
                              base.clean(r["tissue_general_ontology_term_id"]), base.clean(r["tissue_general"]),
                              base.clean(r["disease_ontology_term_id"]), base.clean(r["disease"]),
                              base.clean(r["assay_ontology_term_id"]), base.clean(r["assay"]),
                              base.clean(r["donor_id"]), raw_sum))
            conn.executemany(insert, batch)
    conn.commit()
    conn.execute("CREATE INDEX cells_context ON cells (cell_type_id,tissue_id,disease_id,assay_id)")
    return conn


def expression_for_stratum(census, conn, key, features, batch_size, track_donors):
    acc = base.Accumulator(features, track_donor_covariance=track_donors)
    feature_ids = [f["soma_joinid"] for f in features]
    api = base.census_module()
    for batch in base.cell_batches(conn, key, batch_size):
        ids = np.asarray([r[0] for r in batch], dtype=np.int64)
        cells = {r[0]: r for r in batch}
        adata = api.get_anndata(census, organism="Homo sapiens", measurement_name="RNA", X_name="raw",
                               obs_coords=ids, var_coords=np.asarray(feature_ids, dtype=np.int64),
                               obs_column_names=["soma_joinid"], var_column_names=["soma_joinid", "feature_id"])
        obs_ids = [int(v) for v in adata.obs["soma_joinid"]]
        var_ids = [int(v) for v in adata.var["soma_joinid"]]
        if len(obs_ids) != len(batch) or set(obs_ids) != set(ids) or len(var_ids) != len(feature_ids) or set(var_ids) != set(feature_ids):
            raise RuntimeError("Tabula Sapiens Census coordinates do not match request")
        positions = {jid: j for j, jid in enumerate(var_ids)}
        x = adata.X[:, [positions[jid] for jid in feature_ids]]
        acc.add(x, [cells[j][1] for j in obs_ids], [cells[j][2] for j in obs_ids])
    return acc


def pair_metrics(acc, i, j, thresholds):
    n = acc.n
    d1, d2 = int(acc.n_detected[i]), int(acc.n_detected[j])
    s1, s2 = acc.sum_log_tf[i], acc.sum_log_tf[j]
    v1 = acc.sum_sq_tf[i] - s1 * s1 / n
    v2 = acc.sum_sq_tf[j] - s2 * s2 / n
    covariance = acc.gram[i, j] - s1 * s2 / n
    r = covariance / math.sqrt(v1 * v2) if v1 > 0 and v2 > 0 else 0.0
    codetected = int(acc.codetect[i, j])
    lift = codetected * n / (d1 * d2) if d1 and d2 else 0.0
    donor_r = []
    donor_means = []
    for donor, entry in acc.donors.items():
        dn, dsums = entry[0], entry[1]   # v2 Accumulator stores [n, sum_y, sum_log_lib]
        if dn < thresholds["min_donor_cells"]:
            continue
        donor_means.append((dsums[i] / dn, dsums[j] / dn))
        if donor not in acc.donor_gram:
            continue
        dv1 = acc.donor_sum_sq[donor][i] - dsums[i] * dsums[i] / dn
        dv2 = acc.donor_sum_sq[donor][j] - dsums[j] * dsums[j] / dn
        dcov = acc.donor_gram[donor][i, j] - dsums[i] * dsums[j] / dn
        if dv1 > 0 and dv2 > 0:
            donor_r.append(float(np.clip(dcov / math.sqrt(dv1 * dv2), -1, 1)))
    mean_r = ""
    if len(donor_means) >= 4:
        arr = np.asarray(donor_means)
        if arr[:, 0].std() > 0 and arr[:, 1].std() > 0:
            mean_r = float(np.corrcoef(arr.T)[0, 1])
    positive_fraction = sum(x > 0 for x in donor_r) / len(donor_r) if donor_r else ""
    metrics = {"ts_n_cells": n, "ts_n_donors": len(donor_means), "ts_tf1_detection": d1 / n,
               "ts_tf2_detection": d2 / n, "ts_codetected_cells": codetected, "ts_lift": lift,
               "ts_cell_r": float(np.clip(r, -1, 1)), "ts_donor_mean_r": mean_r,
               "ts_donors_tested": len(donor_r), "ts_positive_donor_fraction": positive_fraction}
    if d1 / n < thresholds["min_detection"] or d2 / n < thresholds["min_detection"] or codetected < thresholds["min_codetected"]:
        return metrics, "inconclusive_low_expression", "Insufficient detection or co-detection"
    if r < thresholds["min_rho"] or lift < thresholds["min_lift"]:
        return metrics, "not_replicated", "Cell-level correlation or co-detection lift below prespecified threshold"
    if len(donor_r) >= thresholds["min_donors"] and positive_fraction >= thresholds["min_positive_donor_fraction"]:
        return metrics, "replicated_donor_supported", "Cell-level thresholds pass and within-donor direction is consistent"
    return metrics, "cell_level_support_only", "Cell-level thresholds pass; donor support insufficient or unavailable"


def result_base(candidate, ts_dataset_id, assay_id="", assay="", n_cells=""):
    return {"candidate_id": candidate_id(candidate), "tf1_ensembl": candidate["tf1_ensembl"],
            "tf1_symbol": candidate["tf1_symbol"], "tf2_ensembl": candidate["tf2_ensembl"],
            "tf2_symbol": candidate["tf2_symbol"], "cell_type_id": candidate["cell_type_id"],
            "cell_type": candidate["cell_type"], "tissue_id": candidate["tissue_id"],
            "tissue": candidate["tissue"], "disease_id": candidate["disease_id"],
            "disease": candidate["disease"], "discovery_assay_id": candidate["assay_id"],
            "discovery_assay": candidate["assay"], "discovery_r": candidate["cell_r_weighted_by_cells"],
            "discovery_lift": candidate["mean_co_detection_lift"],
            "discovery_n_collections": candidate["n_collections_with_hit"],
            "physical_prior_recorded": candidate.get("physical_prior_recorded", ""),
            "ts_dataset_id": ts_dataset_id, "ts_assay_id": assay_id, "ts_assay": assay,
            "ts_n_cells": n_cells}


def run_context(census, conn, candidates, manifest, feature_map, args, writer):
    context = {k: candidates[0][k] for k in CONTEXT_FIELDS}
    cell_type_id, tissue_id, disease_id = (context[k] for k in CONTEXT_FIELDS)
    type_n = conn.execute("SELECT COUNT(*) FROM cells WHERE cell_type_id=?", (cell_type_id,)).fetchone()[0]
    tissue_n = conn.execute("SELECT COUNT(*) FROM cells WHERE cell_type_id=? AND tissue_id=?", (cell_type_id, tissue_id)).fetchone()[0]
    disease_n = conn.execute("SELECT COUNT(*) FROM cells WHERE cell_type_id=? AND tissue_id=? AND disease_id=?", (cell_type_id, tissue_id, disease_id)).fetchone()[0]
    if not type_n or not tissue_n or not disease_n:
        status, reason = (("no_cell_type_match", "No matching Tabula Sapiens cell-type ontology ID") if not type_n else
                          ("no_tissue_match", "Cell type exists but not in the discovery general tissue") if not tissue_n else
                          ("no_disease_match", "No matching disease/condition ontology ID"))
        for candidate in candidates:
            writer.writerow({**result_base(candidate, manifest["ts_dataset_id"]), "validation_status": status, "reason": reason})
        return
    assays = list(conn.execute("SELECT assay_id,MIN(assay),COUNT(*) FROM cells WHERE cell_type_id=? AND tissue_id=? AND disease_id=? GROUP BY assay_id ORDER BY COUNT(*) DESC", (cell_type_id, tissue_id, disease_id)))
    if manifest["thresholds"]["assay_ids"]:
        assays = [a for a in assays if a[0] in manifest["thresholds"]["assay_ids"]]
    if not assays:
        for candidate in candidates:
            writer.writerow({**result_base(candidate, manifest["ts_dataset_id"]), "validation_status": "no_allowed_assay", "reason": "No matching validation assay after filter"})
        return
    for assay_id, assay, n_cells in assays:
        if n_cells < manifest["thresholds"]["min_cells"]:
            for candidate in candidates:
                writer.writerow({**result_base(candidate, manifest["ts_dataset_id"], assay_id, assay, n_cells), "validation_status": "insufficient_cells", "reason": "Fewer than minimum validation cells"})
            continue
        needed = set()
        for candidate in candidates:
            needed.update((candidate["tf1_ensembl"], candidate["tf2_ensembl"]))
        available = [feature_map[g] for g in sorted(needed) if g in feature_map and feature_map[g]["ts_feature_present"]]
        if len(available) < 2:
            for candidate in candidates:
                writer.writerow({**result_base(candidate, manifest["ts_dataset_id"], assay_id, assay, n_cells), "validation_status": "gene_unavailable", "reason": "One or both genes absent from Tabula Sapiens feature set"})
            continue
        key = {"cell_type_id": cell_type_id, "tissue_id": tissue_id, "disease_id": disease_id, "assay_id": assay_id}
        acc = expression_for_stratum(census, conn, key, available, args.batch_size, True)
        if acc.n != n_cells:
            raise RuntimeError(f"Incomplete Tabula Sapiens stratum: {acc.n} versus {n_cells}")
        gene_pos = {f["ensembl_id"]: i for i, f in enumerate(available)}
        for candidate in candidates:
            row = result_base(candidate, manifest["ts_dataset_id"], assay_id, assay, n_cells)
            a, b = candidate["tf1_ensembl"], candidate["tf2_ensembl"]
            if a not in gene_pos or b not in gene_pos:
                writer.writerow({**row, "validation_status": "gene_unavailable", "reason": "One or both genes absent from Tabula Sapiens feature set"})
                continue
            metrics, status, reason = pair_metrics(acc, gene_pos[a], gene_pos[b], manifest["thresholds"])
            writer.writerow({**row, **metrics, "validation_status": status, "reason": reason})


def run(args):
    out = Path(args.out)
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    if args.num_shards < 1 or not 0 <= args.shard_index < args.num_shards or args.batch_size < 1:
        raise ValueError("Invalid shard index/count or batch size")
    contexts = [r for r in base.read_tsv(out / "contexts.tsv") if int(r["context_id"], 16) % args.num_shards == args.shard_index]
    feature_map = {}
    for r in base.read_tsv(out / "features.tsv"):
        if r["soma_joinid"]:
            feature_map[r["ensembl_id"]] = {"ensembl_id": r["ensembl_id"], "symbol": r["symbol"], "soma_joinid": int(r["soma_joinid"]), "is_tf": True, "ts_feature_present": r["ts_feature_present"] == "True"}
    pending = [r for r in contexts if not (out / "groups" / f"{r['context_id']}.done").exists()]
    print(f"Tabula Sapiens {manifest['census_version']}; shard {args.shard_index}/{args.num_shards}; {len(pending)} pending contexts", flush=True)
    if not pending:
        return
    with tempfile.TemporaryDirectory(prefix="ts_validate_", dir=args.tmpdir) as tempdir:
        api = base.census_module()
        with api.open_soma(census_version=manifest["census_version"]) as census:
            conn = metadata_to_sqlite(census, manifest["ts_dataset_id"], {r["cell_type_id"] for r in pending}, tempdir)
            cand_db = sqlite3.connect(f"file:{out / 'candidates.sqlite'}?mode=ro", uri=True)
            for context in pending:
                xid = context["context_id"]
                candidates = [json.loads(r[0]) for r in cand_db.execute("SELECT source_json FROM candidates WHERE context_id=?", (xid,))]
                path = out / "groups" / f"{xid}.tsv.gz"
                with base.tsv_writer(path, RESULT_FIELDS) as writer:
                    run_context(census, conn, candidates, manifest, feature_map, args, writer)
                base.write_json_atomic(out / "groups" / f"{xid}.done", {"context_id": xid, "n_candidates": len(candidates), "census_version": manifest["census_version"]})
                print(f"Completed context {xid}: {len(candidates)} candidates", flush=True)
            cand_db.close()
            conn.close()


def aggregate(args):
    out = Path(args.out)
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    contexts = list(base.read_tsv(out / "contexts.tsv"))
    completed = [r for r in contexts if (out / "groups" / f"{r['context_id']}.done").exists()]
    if not completed:
        raise SystemExit("No completed validation contexts")
    aggregate_dir = out / "aggregate"
    if (aggregate_dir / "replication_results.tsv.gz").exists():
        raise SystemExit("Aggregate already exists; use a new output directory or move the existing aggregate aside")
    rank = {"replicated_donor_supported": 5, "cell_level_support_only": 4, "not_replicated": 3,
            "inconclusive_low_expression": 2, "gene_unavailable": 1, "insufficient_cells": 1,
            "no_allowed_assay": 0, "no_disease_match": 0, "no_tissue_match": 0, "no_cell_type_match": 0}
    conn = sqlite3.connect(str(out / "candidate_summary.sqlite"))
    conn.execute("CREATE TABLE IF NOT EXISTS summary (candidate_id TEXT PRIMARY KEY, best_rank INTEGER, best_json TEXT, n_assays_tested INTEGER, n_assays_reported INTEGER)")
    if conn.execute("SELECT COUNT(*) FROM summary").fetchone()[0]:
        raise SystemExit("Summary database already populated; use a new output directory")
    status_counts = Counter()
    with base.tsv_writer(aggregate_dir / "replication_results.tsv.gz", RESULT_FIELDS) as writer:
        for context in completed:
            path = out / "groups" / f"{context['context_id']}.tsv.gz"
            with base.open_tsv(path) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    writer.writerow(row)
                    status = row["validation_status"]
                    status_counts[status] += 1
                    cid = row["candidate_id"]
                    score = rank.get(status, -1)
                    tested = int(status in ("replicated_donor_supported", "cell_level_support_only", "not_replicated"))
                    existing = conn.execute("SELECT best_rank,best_json,n_assays_tested,n_assays_reported FROM summary WHERE candidate_id=?", (cid,)).fetchone()
                    if existing is None:
                        conn.execute("INSERT INTO summary VALUES (?,?,?,?,?)", (cid, score, json.dumps(row), tested, 1))
                    else:
                        best_rank, best_json, n_tested, n_reported = existing
                        if score > best_rank or (score == best_rank and float(row.get("ts_n_cells") or 0) > float(json.loads(best_json).get("ts_n_cells") or 0)):
                            best_rank, best_json = score, json.dumps(row)
                        conn.execute("UPDATE summary SET best_rank=?,best_json=?,n_assays_tested=?,n_assays_reported=? WHERE candidate_id=?", (best_rank, best_json, n_tested + tested, n_reported + 1, cid))
    summary_fields = [*RESULT_FIELDS, "n_assays_tested", "n_assays_reported"]
    candidate_counts = Counter()
    with base.tsv_writer(aggregate_dir / "candidate_summary.tsv.gz", summary_fields) as writer:
        for raw, n_tested, n_reported in conn.execute("SELECT best_json,n_assays_tested,n_assays_reported FROM summary ORDER BY best_rank DESC,candidate_id"):
            row = json.loads(raw)
            candidate_counts[row["validation_status"]] += 1
            writer.writerow({**row, "n_assays_tested": n_tested, "n_assays_reported": n_reported})
    conn.close()
    base.write_json_atomic(aggregate_dir / "summary.json", {"n_candidates_planned": manifest["n_candidates"],
        "n_contexts_planned": len(contexts), "n_contexts_completed": len(completed),
        "complete": len(completed) == len(contexts), "n_candidates_summarized": sum(candidate_counts.values()),
        "per_assay_status_counts": dict(status_counts), "best_candidate_status_counts": dict(candidate_counts),
        "thresholds": manifest["thresholds"], "ts_dataset_id": manifest["ts_dataset_id"],
        "source_atlas_wide": manifest["source_atlas_wide"]})
    print(f"Aggregated {sum(candidate_counts.values())} candidates across {len(completed)} contexts")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    pp = sub.add_parser("plan")
    pp.add_argument("--discovery-run", required=True)
    pp.add_argument("--out", required=True)
    pp.add_argument("--collection-name", default="Tabula Sapiens")
    pp.add_argument("--all-cells-title", default="Tabula Sapiens - All Cells")
    pp.add_argument("--min-cells", type=int, default=100)
    pp.add_argument("--min-detection", type=float)
    pp.add_argument("--min-rho", type=float)
    pp.add_argument("--min-lift", type=float)
    pp.add_argument("--min-codetected", type=int)
    pp.add_argument("--min-donor-cells", type=int, default=20)
    pp.add_argument("--min-donors", type=int, default=2)
    pp.add_argument("--min-positive-donor-fraction", type=float, default=0.75)
    pp.add_argument("--assay-id", action="append", default=[])
    pp.set_defaults(func=plan)
    rp = sub.add_parser("run")
    rp.add_argument("--out", required=True)
    rp.add_argument("--num-shards", type=int, default=1)
    rp.add_argument("--shard-index", type=int, default=0)
    rp.add_argument("--batch-size", type=int, default=1000)
    rp.add_argument("--tmpdir", default=None)
    rp.set_defaults(func=run)
    ap = sub.add_parser("aggregate")
    ap.add_argument("--out", required=True)
    ap.set_defaults(func=aggregate)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
