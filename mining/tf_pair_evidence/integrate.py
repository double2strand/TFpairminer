#!/usr/bin/env python3
"""Join TFpairminer ChromBPNet pair tests to two independent RNA screens.

The four named outputs are evidence views, not a partition: a pair supported in
both RNA screens appears in both ATAC+RNA files. Untested pairs stay unassessed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path


CONTEXT = ("cell_type_id", "tissue_id", "disease_id")
FIELDS = ("atac_dataset", "atac_context", *CONTEXT, "mapping_level",
          "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
          "atac_status", "cbp_rows", "cbp_best_synergy_q", "cbp_best_synergy",
          "cbp_cooperative", "rna1_status", "rna1_best_r", "rna1_min_q",
          "rna2_status", "rna2_best_r", "rna2_min_q", "evidence_view", "note")


def read_tsv(path: Path):
    handle = (gzip.open(path, "rt", newline="", encoding="utf-8") if path.suffix == ".gz"
              else path.open(newline="", encoding="utf-8"))
    with handle:
        yield from csv.DictReader(handle, delimiter="\t")


def required(row: dict, names: tuple[str, ...], source: Path):
    absent = [name for name in names if name not in row]
    if absent:
        raise ValueError(f"{source}: missing columns {absent}")


def clean_id(value: str) -> str:
    return (value or "").split(".", 1)[0].strip()


def pair(a: str, b: str) -> tuple[str, str]:
    if not a or not b or a == b:
        raise ValueError(f"Invalid gene pair: {a!r}, {b!r}")
    return tuple(sorted((clean_id(a), clean_id(b))))


def number(value: str | None, default: float | None = None):
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def context(row: dict) -> tuple[str, str, str]:
    return tuple((row.get(k) or "").strip() for k in CONTEXT)


def load_register(path: Path):
    by_symbol = defaultdict(set)
    by_id = {}
    for row in read_tsv(path):
        required(row, ("Symbol", "Ensembl_gene_ID"), path)
        symbol, gene = row["Symbol"].strip(), clean_id(row["Ensembl_gene_ID"])
        if symbol and gene:
            by_symbol[symbol.upper()].add(gene)
            by_id[gene] = symbol
    return by_symbol, by_id


def load_rna(root: Path, label: str):
    aggregate = root / "aggregate"
    summary_path = aggregate / "summary.json"
    with summary_path.open(encoding="utf-8") as handle:
        summary = json.load(handle)
    config = summary.get("run_config", {})
    if not summary.get("complete") or int(config.get("max_strata", -1)) != 0:
        raise ValueError(f"{label}: incomplete or capped screen; cannot infer absence of hits")
    if int(config.get("top_k", -1)) != 0:
        raise ValueError(f"{label}: top-k truncation prevents absence-of-hit classification")
    if config.get("primary_stat") not in {"metacell", "pseudobulk"}:
        raise ValueError(f"{label}: donor-aware primary statistic is required")
    if number(config.get("max_q"), 1) > 0.05:
        raise ValueError(f"{label}: q-value cutoff must be <= 0.05")
    with (root / "manifest.json").open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    dataset_path = root / "datasets.tsv"
    dataset_rows = list(read_tsv(dataset_path))
    datasets = {row["dataset_id"] for row in dataset_rows}
    if not datasets:
        raise ValueError(f"{label}: no source datasets")
    min_cells = int(config.get("min_cells", 100))
    min_detection = float(config.get("min_detection", 0.05))
    candidates_path = aggregate / "cell_type_pair_candidates.tsv"
    candidates = defaultdict(list)
    for row in read_tsv(candidates_path):
        required(row, (*CONTEXT, "assay_id", "tf1_ensembl", "tf2_ensembl",
                       "primary_stat", "best_primary_r", "min_q_within_stratum",
                       "mean_co_detection_lift", "co_detected_cells"), candidates_path)
        key = (context(row), pair(row["tf1_ensembl"], row["tf2_ensembl"]))
        candidates[key].append(row)
    return {"label": label, "root": root, "summary": summary, "config": config,
            "manifest": manifest, "datasets": datasets, "dataset_rows": dataset_rows,
            "min_cells": min_cells, "min_detection": min_detection,
            "covered_pairs": set(), "candidates": candidates}


def load_pair_coverage(screen: dict, targets: set):
    """Find same-dataset/assay coverage; aggregate profiles cannot establish this."""
    wanted = defaultdict(set)
    for ctx, genes in targets:
        wanted[ctx].add(genes)
    wanted_genes = {ctx: {g for genes in pairs for g in genes} for ctx, pairs in wanted.items()}
    covered = set()
    for dataset_id in sorted(screen["datasets"]):
        path = screen["root"] / "datasets" / dataset_id / "profiles.tsv.gz"
        if not path.is_file():
            raise ValueError(f"{screen['label']}: missing per-dataset profiles {path}")
        current_key = None
        group = {}

        def flush():
            if current_key is None:
                return
            ctx, _assay = current_key
            for genes in wanted.get(ctx, ()):
                if not all(g in group for g in genes):
                    continue
                a, b = (group[g] for g in genes)
                if (min(int(a["n_cells"]), int(b["n_cells"])) >= screen["min_cells"]
                        and min(float(a["detection_fraction"]), float(b["detection_fraction"]))
                        >= screen["min_detection"]):
                    covered.add((ctx, genes))

        for row in read_tsv(path):
            required(row, (*CONTEXT, "assay_id", "ensembl_id", "n_cells", "detection_fraction"), path)
            key = (context(row), row["assay_id"])
            if key != current_key:
                flush()
                current_key, group = key, {}
            gene = clean_id(row["ensembl_id"])
            if gene in wanted_genes.get(key[0], ()):
                group[gene] = row
        flush()
    screen["covered_pairs"] = covered


def rna_evidence(screen: dict, ctx: tuple[str, str, str], genes: tuple[str, str]):
    if not all(ctx):
        return "unassessed", "", ""
    if (ctx, genes) not in screen["covered_pairs"]:
        return "unassessed", "", ""
    qualifying = []
    for row in screen["candidates"].get((ctx, genes), []):
        method = row["primary_stat"].split("|")
        q = number(row["min_q_within_stratum"])
        r = number(row["best_primary_r"])
        lift = number(row["mean_co_detection_lift"])
        codetected = number(row["co_detected_cells"])
        if (method and set(method) <= {"metacell", "pseudobulk"}
                and q is not None and q <= 0.05 and r is not None and r >= 0.25
                and lift is not None and lift >= 1.1 and codetected is not None
                and codetected >= 25):
            qualifying.append((r, q))
    if qualifying:
        return "hit", max(x[0] for x in qualifying), min(x[1] for x in qualifying)
    if (ctx, genes) in screen["candidates"]:
        return "candidate_nonqualifying", "", ""
    return "screened_no_hit", "", ""


def load_atac(manifest: Path, mapping: Path, by_symbol: dict):
    maps = {}
    for row in read_tsv(mapping):
        required(row, ("atac_dataset", "atac_context", *CONTEXT, "mapping_level"), mapping)
        key = row["atac_dataset"], row["atac_context"]
        if key in maps:
            raise ValueError(f"Duplicate context mapping: {key}")
        if row["mapping_level"] not in {"exact", "proxy", "unmapped"}:
            raise ValueError(f"Invalid mapping_level for {key}")
        maps[key] = row
    events = defaultdict(list)
    provenance = []
    rejects = []
    for source in read_tsv(manifest):
        required(source, ("atac_dataset", "atac_context", "pairs_cbp_path"), manifest)
        key = source["atac_dataset"], source["atac_context"]
        if key not in maps:
            raise ValueError(f"No context mapping for {key}")
        path = Path(source["pairs_cbp_path"])
        if not path.is_absolute():
            path = (manifest.parent / path).resolve()
        if path.name != "pairs_cbp.tsv":
            raise ValueError(f"{path}: require TFpairminer ChromBPNet pairs_cbp.tsv, not a motif proxy")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        provenance.append({"atac_dataset": key[0], "atac_context": key[1],
                           "pairs_cbp_path": str(path), "sha256": digest})
        for row in read_tsv(path):
            required(row, ("tfs", "k", "status", "contribution_supported", "synergy_q"), path)
            names = row["tfs"].split("+")
            if len(names) != 2 or any("::" in name for name in names):
                rejects.append({"source": str(path), "tfs": row["tfs"], "reason": "not_two_unambiguous_TFs"})
                continue
            found = [by_symbol.get(name.strip().upper(), set()) for name in names]
            if any(len(x) != 1 for x in found):
                rejects.append({"source": str(path), "tfs": row["tfs"], "reason": "unknown_or_ambiguous_symbol"})
                continue
            genes = pair(next(iter(found[0])), next(iter(found[1])))
            events[(key, genes)].append(row)
    return maps, events, provenance, rejects


def atac_evidence(rows: list[dict]):
    if not rows:
        return "unassessed", 0, "", "", ""
    tested = [r for r in rows if r["status"] == "tested" and r["k"] == "2"]
    supported = [r for r in tested if r["contribution_supported"].lower() == "yes"]
    if supported:
        best = min(supported, key=lambda r: number(r.get("synergy_q"), 1))
        return "supported", len(rows), best.get("synergy_q", ""), best.get("synergy", ""), best.get("cooperative", "")
    if tested and len(tested) == len(rows) and all(r["contribution_supported"].lower() == "no" for r in tested):
        return "tested_no_support", len(rows), "", "", ""
    return "inconclusive", len(rows), "", "", ""


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict], fields=FIELDS):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def module_rows(rows: list[dict]):
    """Connected TF components, kept separate by evidence view and ATAC context."""
    grouped = defaultdict(list)
    for row in rows:
        for view in row["evidence_view"].split("|"):
            if view == "unassessed":
                continue
            key = (view, row["atac_dataset"], row["atac_context"], context(row))
            grouped[key].append(row)
    modules = []
    for (view, dataset, atac_ctx, rna_ctx), edges in sorted(grouped.items()):
        graph = defaultdict(set)
        symbols = {}
        for edge in edges:
            a, b = edge["tf1_ensembl"], edge["tf2_ensembl"]
            graph[a].add(b)
            graph[b].add(a)
            symbols[a], symbols[b] = edge["tf1_symbol"], edge["tf2_symbol"]
        unseen = set(graph)
        while unseen:
            start = min(unseen)
            stack, component = [start], set()
            while stack:
                gene = stack.pop()
                if gene in component:
                    continue
                component.add(gene)
                stack.extend(graph[gene] - component)
            unseen -= component
            if len(component) < 3:
                continue
            members = sorted(component)
            count = sum(edge["tf1_ensembl"] in component and edge["tf2_ensembl"] in component for edge in edges)
            identity = "|".join((view, dataset, atac_ctx, *rna_ctx, *members))
            modules.append({"module_id": hashlib.sha256(identity.encode()).hexdigest()[:12],
                            "view": view, "atac_dataset": dataset, "atac_context": atac_ctx,
                            **{field: rna_ctx[i] for i, field in enumerate(CONTEXT)},
                            "n_tfs": len(members), "n_pairs": count,
                            "tf_ensembl": ";".join(members),
                            "tf_symbols": ";".join(symbols[g] for g in members)})
    return modules


def integrate(args):
    by_symbol, by_id = load_register(args.register)
    rna1, rna2 = load_rna(args.rna1, "rna1"), load_rna(args.rna2, "rna2")
    if rna1["datasets"] & rna2["datasets"]:
        raise ValueError("RNA1 and RNA2 share source dataset IDs; replication is not independent")
    if {r.get("collection_id") for r in rna1["dataset_rows"] if r.get("collection_id")} & {r.get("collection_id") for r in rna2["dataset_rows"] if r.get("collection_id")}:
        raise ValueError("RNA1 and RNA2 share source collections; replication is not independent")
    if rna1["manifest"].get("dataset_selection") != "all" or not rna1["summary"].get("atlas_wide_except_heldout_collections"):
        raise ValueError("RNA1 must be a complete, uncapped Census screen excluding Tabula Sapiens")
    if not any("tabula sapiens" in str(x).lower() for x in rna1["manifest"].get("excluded_collection_names", [])):
        raise ValueError("RNA1 did not explicitly exclude the Tabula Sapiens collection")
    if any("tabula sapiens" in str(r.get("collection_name", "")).lower() for r in rna1["dataset_rows"]):
        raise ValueError("RNA1 still contains Tabula Sapiens source datasets")
    if (len(rna2["dataset_rows"]) != 1 or rna2["manifest"].get("dataset_selection") != "explicit_ids"
            or "tabula sapiens" not in str(rna2["dataset_rows"][0].get("collection_name", "")).lower()
            or "all cells" not in str(rna2["dataset_rows"][0].get("dataset_title", "")).lower()):
        raise ValueError("RNA2 must be the single Tabula Sapiens All Cells dataset")
    if rna1["manifest"].get("resolved_census_version") != rna2["manifest"].get("resolved_census_version"):
        raise ValueError("RNA screens use different Census releases; pin both to the same build")
    register_digest = sha256(args.register)
    if any(screen["manifest"].get("register_sha256") != register_digest for screen in (rna1, rna2)):
        raise ValueError("RNA screens must use the supplied gene register")
    maps, events, provenance, rejects = load_atac(args.atac_manifest, args.context_map, by_symbol)
    contexts_to_atac = defaultdict(list)
    for atac_key, mapping in maps.items():
        if mapping["mapping_level"] == "exact":
            contexts_to_atac[context(mapping)].append(atac_key)
    tasks = {(atac_key, context(maps[atac_key]) if maps[atac_key]["mapping_level"] == "exact" else ("", "", ""), genes)
             for atac_key, genes in events}
    rna_intersection = set(rna1["candidates"]) & set(rna2["candidates"])
    for ctx, genes in rna_intersection:
        matches = contexts_to_atac.get(ctx) or [("", "")]
        tasks.update((key, ctx, genes) for key in matches)
    targets = {(ctx, genes) for _key, ctx, genes in tasks if all(ctx)}
    load_pair_coverage(rna1, targets)
    load_pair_coverage(rna2, targets)
    output = []
    for atac_key, ctx, genes in sorted(tasks):
        map_row = maps.get(atac_key, {})
        level = map_row.get("mapping_level", "unmapped")
        a_status, count, sy_q, synergy, coop = atac_evidence(events.get((atac_key, genes), []))
        r1_status, r1_r, r1_q = rna_evidence(rna1, ctx, genes)
        r2_status, r2_r, r2_q = rna_evidence(rna2, ctx, genes)
        if level != "exact" and atac_key[0]:
            view, note = "unassessed", "ATAC context is not an exact RNA match"
        elif a_status == "supported" and r1_status == "hit" and r2_status == "hit":
            view, note = "atac_rna1|atac_rna2", "Supported in both independent RNA screens"
        elif a_status == "supported" and r1_status == "hit":
            view, note = "atac_rna1", "RNA2 status retained separately"
        elif a_status == "supported" and r2_status == "hit":
            view, note = "atac_rna2", "RNA1 status retained separately"
        elif a_status == "supported" and r1_status == r2_status == "screened_no_hit":
            view, note = "atac_only", "No significant hit in either complete RNA screen"
        elif a_status == "tested_no_support" and r1_status == r2_status == "hit":
            view, note = "rna_only", "ChromBPNet tested, but no contribution support under this criterion"
        else:
            view, note = "unassessed", "Insufficient evidence for one of the four primary views"
        output.append({"atac_dataset": atac_key[0], "atac_context": atac_key[1],
                       **{field: ctx[i] for i, field in enumerate(CONTEXT)},
                       "mapping_level": level, "tf1_ensembl": genes[0],
                       "tf1_symbol": by_id.get(genes[0], ""), "tf2_ensembl": genes[1],
                       "tf2_symbol": by_id.get(genes[1], ""), "atac_status": a_status,
                       "cbp_rows": count, "cbp_best_synergy_q": sy_q,
                       "cbp_best_synergy": synergy, "cbp_cooperative": coop,
                       "rna1_status": r1_status, "rna1_best_r": r1_r, "rna1_min_q": r1_q,
                       "rna2_status": r2_status, "rna2_best_r": r2_r, "rna2_min_q": r2_q,
                       "evidence_view": view, "note": note})
    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "master.tsv", output)
    views = {"atac_only": "01_atac_only.tsv", "atac_rna1": "02_atac_plus_rna1.tsv",
             "atac_rna2": "03_atac_plus_rna2.tsv", "rna_only": "04_rna_replicated_cbp_no_support.tsv"}
    counts = {}
    for view, filename in views.items():
        rows = [r for r in output if view in r["evidence_view"].split("|")]
        write_tsv(args.out / filename, rows)
        counts[view] = len(rows)
    write_tsv(args.out / "05_unassessed.tsv", [r for r in output if r["evidence_view"] == "unassessed"])
    write_tsv(args.out / "06_atac_plus_both_rna.tsv", [r for r in output if r["evidence_view"] == "atac_rna1|atac_rna2"])
    write_tsv(args.out / "07_rna_replicated_atac_unassessed.tsv",
              [r for r in output if r["rna1_status"] == r["rna2_status"] == "hit"
               and r["atac_status"] in {"unassessed", "inconclusive"}])
    write_tsv(args.out / "modules.tsv", module_rows(output),
              ("module_id", "view", "atac_dataset", "atac_context", *CONTEXT,
               "n_tfs", "n_pairs", "tf_ensembl", "tf_symbols"))
    write_tsv(args.out / "rejected_cbp_labels.tsv", rejects, ("source", "tfs", "reason"))
    with (args.out / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump({"counts": counts, "master_rows": len(output), "rejected_cbp_labels": len(rejects),
                   "chrombpnet_sources": provenance,
                   "rna1": {"path": str(args.rna1), "summary_sha256": sha256(args.rna1 / "aggregate" / "summary.json")},
                   "rna2": {"path": str(args.rna2), "summary_sha256": sha256(args.rna2 / "aggregate" / "summary.json")},
                   "context_map_sha256": sha256(args.context_map),
                   "register_sha256": register_digest,
                   "rna_dataset_overlap": sorted(rna1["datasets"] & rna2["datasets"]),
                   "resolved_census_version": rna1["manifest"].get("resolved_census_version")}, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"master_rows": len(output), "counts": counts}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("atac-manifest", "context-map", "register", "rna1", "rna2", "out"):
        parser.add_argument("--" + flag, type=Path, required=True)
    integrate(parser.parse_args())


if __name__ == "__main__":
    main()
