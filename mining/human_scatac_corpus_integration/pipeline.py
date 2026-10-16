#!/usr/bin/env python3
"""Screen RNA-replicated TF pairs in Human-scATAC-Corpus, independently of CATlas.

Stages: catalog -> fetch -> index -> plan -> scan -> compile -> analyze -> aggregate.
Run as ``python -m human_scatac_corpus_integration.pipeline`` from the repository root.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
import os
import re
import sqlite3
import ssl
import subprocess
import tempfile
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import numpy as np

from catlas_tf_integration.catlas_pipeline import (
    bh, fasta_gc, matched_indices, motif_matches, motif_presence, odds_p,
    parse_meme, read_tsv, sha256, write_json, write_tsv,
)


BASE = "https://health.tsinghua.edu.cn/human-scatac-corpus/"
CATLAS_SOURCE_DATASETS = {"Zhang2021", "Domcke2020"}
COORD_RE = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")
INDEX_FIELDS = ["sample_id", "sample", "dataset", "tissue", "organ", "pmid", "hsc_cell_type",
                "n_cells", "n_donors", "n_batches", "n_cCREs", "universe_sha256",
                "source_stage", "overlaps_catlas_source", "h5ad_sha256"]
PLAN_FIELDS = ["candidate_id", "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
               "cell_type_id", "cell_type", "tissue_id", "tissue", "disease_id", "disease",
               "validation_status", "ts_cell_r", "ts_lift", "ts_n_donors", "sample_id", "sample",
               "dataset", "hsc_cell_type", "hsc_n_cells", "hsc_n_donors", "hsc_n_batches",
               "source_stage", "overlaps_catlas_source", "mapping_tier", "mapping_reviewed",
               "mapping_status", "tf1_motifs", "tf2_motifs", "motif_status"]
RAW_FIELDS = [*PLAN_FIELDS, "n_accessible_cCREs", "n_matched_foreground", "n_matched_background",
              "matching_retention", "tf1_fg_hits", "tf1_bg_hits", "tf2_fg_hits", "tf2_bg_hits",
              "joint_fg_hits", "joint_bg_hits", "tf1_odds_ratio", "tf2_odds_ratio",
              "joint_odds_ratio", "tf1_p", "tf2_p", "joint_p", "joint_lift_fg",
              "joint_lift_bg", "analysis_status"]
FINAL_FIELDS = [*RAW_FIELDS, "tf1_q", "tf2_q", "joint_q", "corpus_result"]
COMBINED_FIELDS = ["candidate_id", "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
                   "cell_type_id", "cell_type", "tissue_id", "tissue", "validation_status",
                   "ts_cell_r", "ts_lift", "ts_n_donors", "catlas_statuses", "n_catlas_tested",
                   "catlas_strict_support", "catlas_provisional_support", "catlas_proxy_support",
                   "corpus_summary",
                   "n_corpus_mapped_samples", "n_corpus_tested_samples",
                   "n_corpus_supporting_datasets", "n_corpus_independent_supporting_datasets",
                   "n_corpus_proxy_supporting_datasets",
                   "combined_atac_support"]


def url_context(insecure):
    return ssl._create_unverified_context() if insecure else None


def first_href(value):
    match = re.search(r"href=['\"]([^'\"]+)", value or "")
    return urllib.parse.urljoin(BASE, match.group(1)) if match else ""


def catalog_rows(insecure=False):
    payload = urllib.parse.urlencode({"filter": "all", "item": "", "start": 0,
                                      "length": 1000, "draw": 1}).encode()
    request = urllib.request.Request(urllib.parse.urljoin(BASE, "getResult_browse.php"), data=payload)
    with urllib.request.urlopen(request, context=url_context(insecure), timeout=120) as response:
        document = json.load(response)
    if len(document.get("Data", [])) != int(document.get("recordsTotal", -1)):
        raise RuntimeError("Corpus catalog response is incomplete; website endpoint may have changed")
    rows = []
    for item in document["Data"]:
        match = re.search(r"Sample\s+(\d+)", item.get("id", ""))
        if not match:
            raise RuntimeError("Unrecognized corpus sample ID")
        rows.append({"sample_id": match.group(1), "sample": item["sample"],
                     "dataset": item["dataset"], "tissue": item["tissue"],
                     "organ": item.get("organ", ""), "n_cells_catalog": item.get("number_of_cells", ""),
                     "pmid": item.get("pmid", ""),
                     "cCRE_url": first_href(item.get("cell_by_cre_matrix", "")),
                     "peak_url": first_href(item.get("cell_by_peak_matrix", "")),
                     "overlaps_catlas_source": "Yes" if item["dataset"] in CATLAS_SOURCE_DATASETS else "No"})
    return rows


def command_catalog(args):
    rows = catalog_rows(args.insecure)
    write_tsv(args.out, rows, ["sample_id", "sample", "dataset", "tissue", "organ",
                               "n_cells_catalog", "pmid", "cCRE_url", "peak_url",
                               "overlaps_catlas_source"])
    print(f"Catalogued {len(rows)} corpus samples: {args.out}")


def normalized_tissue(value):
    return re.sub(r"[^a-z0-9]", "", value.lower())


def command_recommend(args):
    """Suggest tissue-matched corpus samples from RNA candidates, without ATAC gating."""
    candidates = [r for r in read_tsv(args.candidate_summary)
                  if r.get("validation_status") in set(args.include_status)]
    if not candidates:
        raise ValueError("No RNA-replicated candidates selected")
    counts = Counter((r["tissue_id"], r["tissue"]) for r in candidates)
    catalog = list(read_tsv(args.catalog))
    rows = []
    for (tissue_id, tissue), n in sorted(counts.items()):
        desired = normalized_tissue(tissue)
        matches = [(r, "organ_exact" if normalized_tissue(r["organ"]) == desired else "tissue_name_match")
                   for r in catalog if desired and (normalized_tissue(r["organ"]) == desired or
                   desired in normalized_tissue(r["tissue"]))]
        if not matches:
            rows.append({"tissue_id": tissue_id, "tissue": tissue, "n_rna_candidates": n,
                         "match_basis": "no_portal_tissue_match"})
        for r, basis in matches:
            rows.append({"tissue_id": tissue_id, "tissue": tissue, "n_rna_candidates": n,
                         "sample_id": r["sample_id"], "dataset": r["dataset"],
                         "corpus_tissue": r["tissue"], "organ": r["organ"],
                         "n_cells_catalog": r["n_cells_catalog"],
                         "overlaps_catlas_source": r["overlaps_catlas_source"],
                         "match_basis": basis})
    write_tsv(args.out, rows, ["tissue_id", "tissue", "n_rna_candidates", "sample_id",
                               "dataset", "corpus_tissue", "organ", "n_cells_catalog",
                               "overlaps_catlas_source", "match_basis"])
    print(f"Suggested {sum(bool(r.get('sample_id')) for r in rows)} tissue-matched sample rows")


def command_fetch(args):
    rows = {r["sample_id"]: r for r in read_tsv(args.catalog)}
    if not args.sample_id:
        raise ValueError("Select one or more sample IDs from the catalog")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "download_manifest.json"
    previous = {r["sample_id"]: r for r in json.loads(manifest_path.read_text()).get("files", [])} \
        if manifest_path.exists() else {}
    for sample_id in sorted(set(args.sample_id)):
        if sample_id not in rows or not rows[sample_id]["cCRE_url"]:
            raise ValueError(f"No cCRE download for sample {sample_id}")
        dest = out / f"{sample_id}.h5ad"
        if not dest.exists():
            fd, tmpname = tempfile.mkstemp(prefix=f"{sample_id}.", suffix=".part", dir=out)
            try:
                request = urllib.request.Request(rows[sample_id]["cCRE_url"],
                                                 headers={"User-Agent": "RNA-first-TF-pair-integration/1.0"})
                with os.fdopen(fd, "wb") as target, urllib.request.urlopen(
                        request, context=url_context(args.insecure), timeout=120) as source:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        target.write(block)
                os.replace(tmpname, dest)
            except BaseException:
                Path(tmpname).unlink(missing_ok=True)
                raise
        digest = sha256(dest)
        if sample_id in previous and previous[sample_id]["sha256"] != digest:
            raise RuntimeError(f"Previously downloaded sample {sample_id} changed SHA-256")
        previous[sample_id] = {"sample_id": sample_id, "path": str(dest.resolve()),
                               "sha256": digest, "bytes": dest.stat().st_size,
                               "url": rows[sample_id]["cCRE_url"]}
        print(f"Verified local sample {sample_id}: {dest} ({digest})")
    write_json(manifest_path, {"files": [previous[k] for k in sorted(previous)],
                               "catalog_sha256": sha256(args.catalog)})


def coordinate_hash(var_names):
    digest = hashlib.sha256()
    for value in var_names:
        digest.update(str(value).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def cell_column(obs):
    for key in ("Cell_type (HSC)", "cell_type", "celltype", "CellType"):
        if key in obs:
            return key
    raise ValueError("No recognized corpus cell-type column in .obs")


def known_count(values):
    return len({str(v) for v in values if str(v).lower() not in ("", "none", "nan", "na", "unknown")})


def command_index(args):
    catalog = {r["sample_id"]: r for r in read_tsv(args.catalog)}
    paths = sorted(Path(args.data_dir).glob("*.h5ad"))
    if not paths:
        raise ValueError("No sample-id-named H5AD files found")
    rows = []
    for path in paths:
        sample_id = path.stem
        if sample_id not in catalog:
            raise ValueError(f"H5AD filename is not a catalog sample ID: {path}")
        a = ad.read_h5ad(path, backed="r")
        try:
            col = cell_column(a.obs)
            labels = a.obs[col].astype(str)
            donor_col = "Donor (HSC)" if "Donor (HSC)" in a.obs else "donor" if "donor" in a.obs else None
            batch_col = "Batch (HSC)" if "Batch (HSC)" in a.obs else "batch" if "batch" in a.obs else None
            universe = coordinate_hash(a.var_names)
            file_digest = sha256(path)
            if not COORD_RE.fullmatch(str(a.var_names[0])):
                raise ValueError(f"Unexpected cCRE coordinate format in {path}")
            for label, n_cells in labels.value_counts().items():
                mask = labels == label
                rows.append({**{k: catalog[sample_id].get(k, "") for k in
                                ("sample_id", "sample", "dataset", "tissue", "organ", "pmid",
                                 "overlaps_catlas_source")},
                             "hsc_cell_type": label, "n_cells": int(n_cells),
                             "n_donors": known_count(a.obs.loc[mask, donor_col]) if donor_col else 0,
                             "n_batches": known_count(a.obs.loc[mask, batch_col]) if batch_col else 0,
                             "n_cCREs": a.n_vars, "universe_sha256": universe,
                             "source_stage": "fetal" if catalog[sample_id]["dataset"] == "Domcke2020" else
                                             "adult" if catalog[sample_id]["dataset"] == "Zhang2021" else "unknown",
                             "h5ad_sha256": file_digest})
        finally:
            a.file.close()
        print(f"Indexed sample {sample_id}: {len(set(r['hsc_cell_type'] for r in rows if r['sample_id'] == sample_id))} labels")
    write_tsv(args.out, rows, INDEX_FIELDS)


def read_crosswalk(path, index):
    if not path:
        return []
    valid = {(r["sample_id"], r["hsc_cell_type"]): r for r in index}
    rows = list(read_tsv(path))
    seen = set()
    for r in rows:
        key = (r["sample_id"], r["hsc_cell_type"])
        if key not in valid:
            raise ValueError(f"Crosswalk references an unindexed sample/cell type: {key}")
        if r["mapping_tier"] not in ("exact_reviewed", "proxy_fetal", "proxy_related_celltype", "proxy_disease_context"):
            raise ValueError("Invalid mapping_tier")
        if r["reviewed"].lower() not in ("yes", "no"):
            raise ValueError("Crosswalk reviewed must be Yes or No")
        if r["mapping_tier"] == "exact_reviewed" and (not r["tissue_id"] or r["reviewed"].lower() != "yes"):
            raise ValueError("Exact mappings require a tissue ID and explicit review")
        if r["mapping_tier"] == "proxy_fetal" and valid[key]["source_stage"] != "fetal":
            raise ValueError("proxy_fetal points to a non-fetal sample")
        identity = (r["cell_type_id"], r["tissue_id"], *key)
        if identity in seen:
            raise ValueError(f"Duplicate crosswalk row: {identity}")
        seen.add(identity)
    return rows


def command_plan(args):
    out = Path(args.out)
    if (out / "manifest.json").exists():
        raise RuntimeError("Plan already exists; choose a new run directory")
    out.mkdir(parents=True, exist_ok=True)
    candidates = [r for r in read_tsv(args.candidate_summary) if r.get("validation_status") in set(args.include_status)]
    if not candidates:
        raise ValueError("No RNA-replicated candidates selected")
    if len({r["candidate_id"] for r in candidates}) != len(candidates):
        raise ValueError("Duplicate RNA candidate IDs")
    index = list(read_tsv(args.index))
    if not index:
        raise ValueError("Empty corpus cell-type index")
    by_sample_label = {(r["sample_id"], r["hsc_cell_type"]): r for r in index}
    crosswalk = read_crosswalk(args.crosswalk, index)
    header, blocks = parse_meme(args.meme)
    genes = {(r["tf1_ensembl"], r["tf1_symbol"]) for r in candidates} | {
        (r["tf2_ensembl"], r["tf2_symbol"]) for r in candidates}
    gene_motifs = {}
    motif_rows = []
    selected = set()
    for gene_id, symbol in sorted(genes):
        direct, composite = motif_matches(symbol, blocks)
        gene_motifs[gene_id] = set(direct)
        selected.update(direct)
        motif_rows.append({"ensembl_id": gene_id, "symbol": symbol,
                           "motif_ids": ",".join(direct), "composite_ids": ",".join(composite),
                           "motif_status": "exact_single_tf" if direct else "no_exact_single_tf_motif"})
    write_tsv(out / "motif_map.tsv", motif_rows,
              ["ensembl_id", "symbol", "motif_ids", "composite_ids", "motif_status"])
    review = []
    context_set = sorted({(r["cell_type_id"], r["cell_type"], r["tissue_id"], r["tissue"]) for r in candidates})
    for cell_id, cell_name, tissue_id, tissue in context_set:
        scored = [(difflib.SequenceMatcher(None, cell_name.lower(), r["hsc_cell_type"].lower()).ratio(), r)
                  for r in index]
        nearby = [r for score, r in sorted(scored, key=lambda item: (
            item[1]["organ"].lower() != tissue.lower(), -item[0])) if score >= 0.6][:5]
        review.append({"cell_type_id": cell_id, "cell_type": cell_name, "tissue_id": tissue_id,
                       "tissue": tissue, "suggested_sample_labels": ";".join(
                           f"{r['sample_id']}:{r['hsc_cell_type']}:{r['source_stage']}" for r in nearby)})
    write_tsv(out / "mapping_review.tsv", review,
              ["cell_type_id", "cell_type", "tissue_id", "tissue", "suggested_sample_labels"])
    plan = []
    for c in candidates:
        matching = [r for r in crosswalk if r["cell_type_id"] == c["cell_type_id"] and
                    r["tissue_id"] == c["tissue_id"] and r["reviewed"].lower() == "yes"]
        healthy = (c.get("disease_id", "") in ("", "PATO:0000461") and
                   c.get("disease", "").lower() in ("", "normal", "healthy"))
        for mapping in matching or [None]:
            if mapping and not healthy and mapping["mapping_tier"] == "exact_reviewed":
                raise ValueError("Disease-context RNA candidate cannot be an exact match to healthy ATAC without review of disease status")
            indexed = by_sample_label[(mapping["sample_id"], mapping["hsc_cell_type"])] if mapping else None
            a, b = gene_motifs[c["tf1_ensembl"]], gene_motifs[c["tf2_ensembl"]]
            plan.append({**c, "sample_id": mapping["sample_id"] if mapping else "",
                         "sample": indexed["sample"] if indexed else "",
                         "dataset": indexed["dataset"] if indexed else "",
                         "hsc_cell_type": mapping["hsc_cell_type"] if mapping else "",
                         "hsc_n_cells": indexed["n_cells"] if indexed else "",
                         "hsc_n_donors": indexed["n_donors"] if indexed else "",
                         "hsc_n_batches": indexed["n_batches"] if indexed else "",
                         "source_stage": indexed["source_stage"] if indexed else "",
                         "overlaps_catlas_source": indexed["overlaps_catlas_source"] if indexed else "",
                         "mapping_tier": mapping["mapping_tier"] if mapping else "",
                         "mapping_reviewed": mapping["reviewed"] if mapping else "No",
                         "mapping_status": "mapped" if mapping else "no_reviewed_corpus_mapping",
                         "tf1_motifs": ",".join(sorted(a)), "tf2_motifs": ",".join(sorted(b)),
                         "motif_status": "both_exact_motifs" if a and b and not (a & b) else
                                         "shared_motif_ambiguous" if a & b else "motif_unavailable"})
    mapped = [r for r in plan if r["sample_id"]]
    universes = {by_sample_label[(r["sample_id"], r["hsc_cell_type"])]["universe_sha256"] for r in mapped}
    sizes = {int(by_sample_label[(r["sample_id"], r["hsc_cell_type"])]["n_cCREs"]) for r in mapped}
    if len(universes) > 1 or len(sizes) > 1:
        raise RuntimeError("Mapped samples do not share one cCRE coordinate universe; split into separate runs")
    write_tsv(out / "plan.tsv", plan, PLAN_FIELDS)
    with open(out / "selected_motifs.meme", "w", encoding="utf-8") as fh:
        fh.write(header)
        for motif_id in sorted(selected):
            fh.write(blocks[motif_id][1])
    manifest = {"selection_source": "Tabula Sapiens RNA replication only; CATlas is not an input or filter",
                "candidate_summary_sha256": sha256(args.candidate_summary), "corpus_index_sha256": sha256(args.index),
                "crosswalk_sha256": sha256(args.crosswalk) if args.crosswalk else "",
                "jaspar_meme_sha256": sha256(args.meme), "n_candidates": len(candidates),
                "n_plan_rows": len(plan), "n_mapped_rows": len(mapped),
                "n_unmapped_rows": len(plan) - len(mapped), "n_selected_motifs": len(selected),
                "n_peaks": next(iter(sizes)) if sizes else 0,
                "n_cCREs": next(iter(sizes)) if sizes else 0,
                "n_distal_peaks": next(iter(sizes)) if sizes else 0,
                "feature_space": "all corpus cCREs; n_distal_peaks is a compile-helper compatibility field",
                "universe_sha256": next(iter(universes)) if universes else "",
                "mapped_sample_ids": sorted({r["sample_id"] for r in mapped}),
                "n_catlas_source_overlap_rows": sum(r["overlaps_catlas_source"] == "Yes" for r in mapped)}
    write_json(out / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


def load_plan(out):
    out = Path(out)
    return json.loads((out / "manifest.json").read_text()), list(read_tsv(out / "plan.tsv"))


def parse_coordinate(value):
    match = COORD_RE.fullmatch(str(value))
    if not match:
        raise ValueError(f"Unrecognized cCRE coordinate: {value}")
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    if end <= start:
        raise ValueError(f"Invalid cCRE interval: {value}")
    return chrom, start, end


def coordinate_width(value):
    _, start, end = parse_coordinate(value)
    return end - start


def command_scan(args):
    out = Path(args.out)
    manifest, _ = load_plan(out)
    if not manifest["n_mapped_rows"]:
        raise RuntimeError("No reviewed corpus mappings; aggregate a coverage-only result")
    if not 0 <= args.shard_index < args.num_shards or args.num_shards < 1:
        raise ValueError("Invalid scan shard")
    if not 0 < args.fimo_p <= 1 or args.max_stored_scores < 100000:
        raise ValueError("Invalid FIMO parameters")
    if manifest["n_selected_motifs"] == 0:
        raise ValueError("No direct motifs in the screened pair set")
    marker = out / "scan" / f"{args.shard_index:04d}.done.json"
    if marker.exists():
        prior = json.loads(marker.read_text())
        if (prior["num_shards"], prior["fimo_p"], prior["max_stored_scores"]) != (
                args.num_shards, args.fimo_p, args.max_stored_scores):
            raise RuntimeError("Existing scan shard has different parameters")
        print(f"Scan shard {args.shard_index} already complete")
        return
    sample_id = manifest["mapped_sample_ids"][0]
    a = ad.read_h5ad(Path(args.data_dir) / f"{sample_id}.h5ad", backed="r")
    try:
        if a.n_vars != manifest["n_peaks"] or coordinate_hash(a.var_names) != manifest["universe_sha256"]:
            raise RuntimeError("Corpus cCRE universe differs from the plan")
        indices = np.arange(args.shard_index, a.n_vars, args.num_shards)
        with tempfile.TemporaryDirectory(prefix="hsc_scan_", dir=args.tmpdir) as tmp:
            bed, fasta = Path(tmp) / "cCREs.bed", Path(tmp) / "cCREs.fa"
            with open(bed, "w", encoding="utf-8") as fh:
                for peak in indices:
                    chrom, start, end = parse_coordinate(a.var_names[peak])
                    fh.write(f"{chrom}\t{start}\t{end}\tpeak_{peak}\n")
            with open(fasta, "w", encoding="utf-8") as fh:
                subprocess.run([args.bedtools, "getfasta", "-fi", args.genome_fasta,
                                "-bed", str(bed), "-nameOnly"], stdout=fh, check=True)
            gc_rows = []
            for name, seq in fasta_gc(fasta):
                valid = sum(base in "ACGT" for base in seq)
                gc_rows.append({"peak_index": int(name.removeprefix("peak_")),
                                "gc": sum(base in "GC" for base in seq) / valid
                                if valid and valid / len(seq) >= 0.8 else ""})
            if len(gc_rows) != len(indices):
                raise RuntimeError("Genome FASTA did not yield all cCRE sequences; check hg38 contig names")
            conn = sqlite3.connect(str(Path(tmp) / "hits.sqlite"))
            conn.execute("CREATE TABLE hits (motif_id TEXT, peak_index INTEGER, PRIMARY KEY(motif_id,peak_index)) WITHOUT ROWID")
            stderr_path = Path(tmp) / "fimo.stderr"
            cmd = [args.fimo, "--text", "--no-qvalue", "--thresh", str(args.fimo_p),
                   "--max-stored-scores", str(args.max_stored_scores),
                   str(out / "selected_motifs.meme"), str(fasta)]
            with open(stderr_path, "w", encoding="utf-8") as stderr_fh:
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr_fh, text=True)
                n_sites = 0
                try:
                    for line in proc.stdout:
                        if not line.strip() or line.startswith("#") or line.startswith("motif_id"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 3 or not parts[2].startswith("peak_"):
                            raise RuntimeError(f"Unexpected FIMO output: {line[:200]}")
                        conn.execute("INSERT OR IGNORE INTO hits VALUES (?,?)", (parts[0], int(parts[2][5:])))
                        n_sites += 1
                    code = proc.wait()
                except BaseException:
                    proc.kill()
                    proc.wait()
                    raise
            err = stderr_path.read_text()
            if code != 0 or re.search(r"max(?:imum)?[- ]stored[- ]scores|discard(?:ing|ed).*scores", err, re.I):
                raise RuntimeError(f"FIMO failed or truncated its output: {err[-2000:]}")
            write_tsv(out / "scan" / f"{args.shard_index:04d}.hits.tsv.gz",
                      ({"motif_id": mid, "peak_index": peak} for mid, peak in
                       conn.execute("SELECT motif_id,peak_index FROM hits")),
                      ["motif_id", "peak_index"])
            conn.close()
            write_tsv(out / "scan" / f"{args.shard_index:04d}.gc.tsv.gz", gc_rows, ["peak_index", "gc"])
    finally:
        a.file.close()
    write_json(marker, {"shard_index": args.shard_index, "num_shards": args.num_shards,
                        "fimo_p": args.fimo_p, "max_stored_scores": args.max_stored_scores,
                        "n_distal_peaks": len(indices), "n_fimo_sites": n_sites})
    print(f"Scanned {len(indices)} cCREs in shard {args.shard_index}")


def command_compile(args):
    # Same packed motif-presence representation as CATlas; input selection remains RNA-only.
    from catlas_tf_integration.catlas_pipeline import command_compile as shared_compile
    shared_compile(args)


def accessibility_counts(a, row_indices, batch_size):
    counts = np.zeros(a.n_vars, dtype=np.int32)
    for offset in range(0, len(row_indices), batch_size):
        x = a.X[row_indices[offset:offset + batch_size], :]
        counts += np.asarray((x > 0).sum(axis=0)).ravel().astype(np.int32)
    return counts


def command_analyze(args):
    out = Path(args.out)
    manifest, plan = load_plan(out)
    if not 0 <= args.shard_index < args.num_shards or args.num_shards < 1:
        raise ValueError("Invalid analysis shard")
    if (args.min_cells < 1 or args.min_open_cells < 1 or
            not 0 < args.min_open_fraction <= 1 or args.min_matched < 1 or
            args.min_motif_hits < 1 or args.cell_batch_size < 1):
        raise ValueError("Invalid analysis thresholds")
    if not (out / "compiled" / "manifest.json").exists():
        raise RuntimeError("Compile motif scans first")
    groups = defaultdict(list)
    for row in plan:
        if row["sample_id"]:
            groups[(row["sample_id"], row["hsc_cell_type"])].append(row)
    assigned = sorted(key for key in groups if
                      int(hashlib.sha256("|".join(key).encode()).hexdigest(), 16) % args.num_shards == args.shard_index)
    marker = out / "analysis" / f"{args.shard_index:04d}.done.json"
    config = {"num_shards": args.num_shards, "min_cells": args.min_cells,
              "min_open_cells": args.min_open_cells, "min_open_fraction": args.min_open_fraction,
              "min_matched": args.min_matched, "min_motif_hits": args.min_motif_hits,
              "cell_batch_size": args.cell_batch_size}
    if marker.exists():
        if any(json.loads(marker.read_text()).get(k) != v for k, v in config.items()):
            raise RuntimeError("Existing analysis shard has different parameters")
        print(f"Analysis shard {args.shard_index} already complete")
        return
    packed = np.load(out / "compiled" / "hits.npy", mmap_mode="r")
    gc = np.load(out / "compiled" / "gc.npy", mmap_mode="r")
    motif_index = {r["motif_id"]: int(r["row_index"]) for r in read_tsv(out / "compiled" / "motif_ids.tsv")}
    by_sample = defaultdict(list)
    for key in assigned:
        by_sample[key[0]].append(key)
    results = []
    for sample_id, keys in by_sample.items():
        a = ad.read_h5ad(Path(args.data_dir) / f"{sample_id}.h5ad", backed="r")
        try:
            if a.n_vars != manifest["n_peaks"] or coordinate_hash(a.var_names) != manifest["universe_sha256"]:
                raise RuntimeError(f"cCRE universe changed for sample {sample_id}")
            labels = a.obs[cell_column(a.obs)].astype(str).to_numpy()
            width = np.fromiter((coordinate_width(v) for v in a.var_names),
                                dtype=np.int32, count=a.n_vars)
            valid = np.isfinite(gc) & (width > 0)
            for key in keys:
                label = key[1]
                cell_ids = np.flatnonzero(labels == label)
                if len(cell_ids) < args.min_cells:
                    for row in groups[key]:
                        results.append({**row, "analysis_status": "insufficient_cells"})
                    continue
                counts = accessibility_counts(a, cell_ids, args.cell_batch_size)
                cutoff = max(args.min_open_cells, math.ceil(args.min_open_fraction * len(cell_ids)))
                foreground = counts >= cutoff
                eligible = valid & (foreground | (counts == 0))
                n_accessible = int(np.count_nonzero(foreground & valid))
                if n_accessible == 0 or not np.any(eligible & ~foreground):
                    for row in groups[key]:
                        results.append({**row, "n_accessible_cCREs": n_accessible,
                                        "analysis_status": "insufficient_matched_cCREs"})
                    continue
                seed = int(hashlib.sha256(f"{sample_id}|{label}|HSC".encode()).hexdigest()[:16], 16)
                fg, bg = matched_indices(foreground, eligible, gc, width, seed)
                retention = len(fg) / n_accessible if n_accessible else 0.0
                cache = {}
                for row in groups[key]:
                    base = {**row, "n_accessible_cCREs": n_accessible,
                            "n_matched_foreground": len(fg), "n_matched_background": len(bg),
                            "matching_retention": retention}
                    if len(fg) < args.min_matched:
                        results.append({**base, "analysis_status": "insufficient_matched_cCREs"})
                        continue
                    if row["motif_status"] != "both_exact_motifs":
                        results.append({**base, "analysis_status": row["motif_status"]})
                        continue
                    for gene, motifs in ((row["tf1_ensembl"], row["tf1_motifs"]),
                                         (row["tf2_ensembl"], row["tf2_motifs"])):
                        if gene not in cache:
                            cache[gene] = motif_presence(motifs, packed, motif_index, manifest["n_peaks"])
                    left, right = cache[row["tf1_ensembl"]], cache[row["tf2_ensembl"]]
                    af, ab = int(left[fg].sum()), int(left[bg].sum())
                    bf, bb = int(right[fg].sum()), int(right[bg].sum())
                    joint = left & right
                    jf, jb = int(joint[fg].sum()), int(joint[bg].sum())
                    counts_out = {"tf1_fg_hits": af, "tf1_bg_hits": ab, "tf2_fg_hits": bf,
                                  "tf2_bg_hits": bb, "joint_fg_hits": jf, "joint_bg_hits": jb}
                    if min(af, bf) < args.min_motif_hits:
                        results.append({**base, **counts_out, "analysis_status": "inconclusive_sparse_motif"})
                        continue
                    odds_a, p_a = odds_p(af, len(fg), ab, len(bg))
                    odds_b, p_b = odds_p(bf, len(fg), bb, len(bg))
                    odds_j, p_j = odds_p(jf, len(fg), jb, len(bg))
                    results.append({**base, **counts_out,
                                    "tf1_odds_ratio": odds_a, "tf2_odds_ratio": odds_b,
                                    "joint_odds_ratio": odds_j, "tf1_p": p_a, "tf2_p": p_b,
                                    "joint_p": p_j,
                                    "joint_lift_fg": jf * len(fg) / (af * bf) if af and bf else 0,
                                    "joint_lift_bg": jb * len(bg) / (ab * bb) if ab and bb else 0,
                                    "analysis_status": "tested"})
        finally:
            a.file.close()
    write_tsv(out / "analysis" / f"{args.shard_index:04d}.tsv.gz", results, RAW_FIELDS)
    write_json(marker, {**config, "shard_index": args.shard_index, "n_rows": len(results)})
    print(f"Analyzed {len(assigned)} sample–cell-type groups; {len(results)} pair rows")


def command_aggregate(args):
    out = Path(args.out)
    manifest, plan = load_plan(out)
    if not 0 < args.q <= 1 or args.min_joint < 1:
        raise ValueError("Invalid aggregate thresholds")
    results = []
    if manifest["n_mapped_rows"]:
        markers = [json.loads(p.read_text()) for p in sorted((out / "analysis").glob("*.done.json"))]
        if not markers:
            raise RuntimeError("No analysis shards complete")
        n_shards = markers[0]["num_shards"]
        if len(markers) != n_shards or {m["shard_index"] for m in markers} != set(range(n_shards)):
            raise RuntimeError("Analysis shards incomplete")
        config_fields = ("num_shards", "min_cells", "min_open_cells", "min_open_fraction",
                         "min_matched", "min_motif_hits", "cell_batch_size")
        if len({tuple(m[f] for f in config_fields) for m in markers}) != 1:
            raise RuntimeError("Analysis shards used different parameters")
        for i in range(n_shards):
            results.extend(read_tsv(out / "analysis" / f"{i:04d}.tsv.gz"))
    else:
        n_shards = 0
    results.extend({**r, "analysis_status": r["mapping_status"]} for r in plan if not r["sample_id"])
    expected = {(r["candidate_id"], r["sample_id"], r["hsc_cell_type"]) for r in plan}
    observed = [(r["candidate_id"], r["sample_id"], r["hsc_cell_type"]) for r in results]
    if len(expected) != len(plan) or len(observed) != len(set(observed)) or set(observed) != expected:
        raise RuntimeError("Result rows do not match RNA-first plan")
    tested = [r for r in results if r["analysis_status"] == "tested"]
    gene_p = {}
    for r in tested:
        for gene, p in ((r["tf1_ensembl"], r["tf1_p"]), (r["tf2_ensembl"], r["tf2_p"])):
            key = (r["sample_id"], r["hsc_cell_type"], gene)
            if key in gene_p and abs(gene_p[key] - float(p)) > 1e-10:
                raise RuntimeError("Inconsistent gene p-value in one sample/cell type")
            gene_p[key] = float(p)
    gene_keys = sorted(gene_p)
    gene_q = dict(zip(gene_keys, bh([gene_p[k] for k in gene_keys])))
    for r, q_joint in zip(tested, bh([float(r["joint_p"]) for r in tested])):
        r["tf1_q"] = gene_q[(r["sample_id"], r["hsc_cell_type"], r["tf1_ensembl"])]
        r["tf2_q"] = gene_q[(r["sample_id"], r["hsc_cell_type"], r["tf2_ensembl"])]
        r["joint_q"] = q_joint
    counts = Counter()
    for r in results:
        if r["analysis_status"] != "tested":
            result = r["analysis_status"]
        else:
            passed = (all(float(r[k]) <= args.q for k in ("tf1_q", "tf2_q", "joint_q")) and
                      all(float(r[k]) > 1 for k in ("tf1_odds_ratio", "tf2_odds_ratio", "joint_odds_ratio")) and
                      int(r["joint_fg_hits"]) >= args.min_joint and
                      float(r["joint_lift_fg"]) > float(r["joint_lift_bg"]))
            if r["mapping_tier"] == "exact_reviewed":
                result = ("corpus_coenriched_shared_catlas_source" if r["overlaps_catlas_source"] == "Yes"
                          else "corpus_coenriched_independent_source") if passed else "not_coenriched"
            else:
                result = "proxy_coenriched" if passed else "proxy_not_coenriched"
        r["corpus_result"] = result
        counts[result] += 1
    agg = out / "aggregate"
    if (agg / "corpus_pair_results.tsv.gz").exists():
        raise RuntimeError("Aggregate already exists; use a new run directory")
    write_tsv(agg / "corpus_pair_results.tsv.gz", results, FINAL_FIELDS)
    rollup = []
    for candidate_id in sorted({r["candidate_id"] for r in plan}):
        rows = [r for r in results if r["candidate_id"] == candidate_id]
        supported = [r for r in rows if r["corpus_result"].startswith("corpus_coenriched_")]
        proxy_supported = [r for r in rows if r["corpus_result"] == "proxy_coenriched"]
        rollup.append({"candidate_id": candidate_id, "tf1_symbol": rows[0]["tf1_symbol"],
                       "tf2_symbol": rows[0]["tf2_symbol"], "rna_validation_status": rows[0]["validation_status"],
                       "n_corpus_mapped_samples": len({r["sample_id"] for r in rows if r["sample_id"]}),
                       "n_corpus_tested_samples": len({r["sample_id"] for r in rows if r["analysis_status"] == "tested"}),
                       "n_supporting_datasets": len({r["dataset"] for r in supported}),
                       "n_independent_supporting_datasets": len({r["dataset"] for r in supported
                                                                 if r["overlaps_catlas_source"] == "No"}),
                       "n_proxy_supporting_datasets": len({r["dataset"] for r in proxy_supported}),
                       "corpus_summary": "independent_source_support" if any(
                           r["corpus_result"] == "corpus_coenriched_independent_source" for r in supported)
                           else "shared_source_support_only" if supported else
                           "proxy_support_only" if proxy_supported else
                           "no_support_or_unresolved"})
    write_tsv(agg / "candidate_summary.tsv", rollup,
              ["candidate_id", "tf1_symbol", "tf2_symbol", "rna_validation_status",
               "n_corpus_mapped_samples", "n_corpus_tested_samples", "n_supporting_datasets",
               "n_independent_supporting_datasets", "n_proxy_supporting_datasets",
               "corpus_summary"])
    write_json(agg / "summary.json", {"selection_source": manifest["selection_source"],
                "n_input_rna_candidates": manifest["n_candidates"], "n_plan_rows": len(plan),
                "n_result_rows": len(results), "n_tested_rows": len(tested),
                "n_analysis_shards": n_shards, "result_counts": dict(counts),
                "q_threshold": args.q, "min_joint_cCREs": args.min_joint,
                "interpretation": "Motif co-enrichment in accessible corpus cCREs; not direct co-occupancy or TF interaction"})
    print(json.dumps(dict(counts), indent=2))


def command_combine(args):
    """Left-join independently generated ATAC evidence onto every replicated RNA hit."""
    candidates = [r for r in read_tsv(args.candidate_summary)
                  if r.get("validation_status") in set(args.include_status)]
    selected_ids = {r["candidate_id"] for r in candidates}
    if len(selected_ids) != len(candidates) or not candidates:
        raise ValueError("RNA candidate IDs must be unique and nonempty")
    catlas = defaultdict(list)
    if args.catlas_results:
        for r in read_tsv(args.catlas_results):
            if r["candidate_id"] not in selected_ids:
                raise ValueError(f"CATlas candidate not in selected RNA set: {r['candidate_id']}")
            catlas[r["candidate_id"]].append(r)
    corpus = {}
    if args.corpus_summary:
        for r in read_tsv(args.corpus_summary):
            if r["candidate_id"] not in selected_ids:
                raise ValueError(f"Corpus candidate not in selected RNA set: {r['candidate_id']}")
            if r["candidate_id"] in corpus:
                raise ValueError(f"Duplicate corpus summary candidate: {r['candidate_id']}")
            corpus[r["candidate_id"]] = r
    rows = []
    for candidate in candidates:
        candidate_id = candidate["candidate_id"]
        crows = catlas[candidate_id]
        statuses = sorted({r.get("catlas_result", "") for r in crows if r.get("catlas_result")})
        cinfo = corpus.get(candidate_id, {})
        strict = "catlas_coenriched" in statuses
        proxy = "proxy_coenriched" in statuses
        provisional = "provisional_coenriched" in statuses
        independent = cinfo.get("corpus_summary") == "independent_source_support"
        shared = cinfo.get("corpus_summary") == "shared_source_support_only"
        corpus_proxy = cinfo.get("corpus_summary") == "proxy_support_only"
        if strict and independent:
            combined = "two_atac_resources_distinct_studies"
        elif independent:
            combined = "corpus_independent_source_only"
        elif strict:
            combined = "catlas_strict_only"
        elif proxy or provisional or shared or corpus_proxy:
            combined = "proxy_provisional_or_shared_source_support_only"
        elif any(r.get("analysis_status") == "tested" for r in crows) or \
                int(cinfo.get("n_corpus_tested_samples") or 0) > 0:
            combined = "atac_tested_no_support"
        else:
            combined = "no_matched_atac_test"
        rows.append({**candidate, "catlas_statuses": ";".join(statuses) if statuses else "not_assessed",
                     "n_catlas_tested": sum(r.get("analysis_status") == "tested" for r in crows),
                     "catlas_strict_support": "Yes" if strict else "No",
                     "catlas_provisional_support": "Yes" if provisional else "No",
                     "catlas_proxy_support": "Yes" if proxy else "No",
                     "corpus_summary": cinfo.get("corpus_summary", "not_assessed"),
                     "n_corpus_mapped_samples": cinfo.get("n_corpus_mapped_samples", ""),
                     "n_corpus_tested_samples": cinfo.get("n_corpus_tested_samples", ""),
                     "n_corpus_supporting_datasets": cinfo.get("n_supporting_datasets", ""),
                     "n_corpus_independent_supporting_datasets": cinfo.get(
                         "n_independent_supporting_datasets", ""),
                     "n_corpus_proxy_supporting_datasets": cinfo.get(
                         "n_proxy_supporting_datasets", ""),
                     "combined_atac_support": combined})
    out = Path(args.out)
    if out.exists():
        raise RuntimeError("Combined evidence table already exists; choose a new output path")
    write_tsv(out, rows, COMBINED_FIELDS)
    write_json(str(out) + ".manifest.json", {
        "selection_source": "RNA replication only; no CATlas or corpus result used to select pairs",
        "candidate_summary_sha256": sha256(args.candidate_summary),
        "catlas_results_sha256": sha256(args.catlas_results) if args.catlas_results else "",
        "corpus_summary_sha256": sha256(args.corpus_summary) if args.corpus_summary else "",
        "n_rna_candidates": len(rows),
        "combined_counts": dict(Counter(r["combined_atac_support"] for r in rows))})
    print(f"Combined independent evidence columns for {len(rows)} RNA-selected pairs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("catalog")
    p.add_argument("--out", required=True)
    p.add_argument("--insecure", action="store_true", help="Opt-in TLS bypass for a known site certificate issue")
    p.set_defaults(func=command_catalog)
    p = sub.add_parser("fetch")
    p.add_argument("--catalog", required=True)
    p.add_argument("--sample-id", action="append", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--insecure", action="store_true")
    p.set_defaults(func=command_fetch)
    p = sub.add_parser("recommend")
    p.add_argument("--candidate-summary", required=True)
    p.add_argument("--catalog", required=True)
    p.add_argument("--include-status", action="append", default=["replicated_donor_supported"])
    p.add_argument("--out", required=True)
    p.set_defaults(func=command_recommend)
    p = sub.add_parser("index")
    p.add_argument("--catalog", required=True)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=command_index)
    p = sub.add_parser("plan")
    p.add_argument("--candidate-summary", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--crosswalk")
    p.add_argument("--meme", required=True)
    p.add_argument("--include-status", action="append", default=["replicated_donor_supported"])
    p.add_argument("--out", required=True)
    p.set_defaults(func=command_plan)
    p = sub.add_parser("scan")
    p.add_argument("--out", required=True)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--fimo", default="fimo")
    p.add_argument("--bedtools", default="bedtools")
    p.add_argument("--fimo-p", type=float, default=1e-4)
    p.add_argument("--max-stored-scores", type=int, default=5000000)
    p.add_argument("--num-shards", type=int, default=32)
    p.add_argument("--shard-index", type=int, default=0)
    p.add_argument("--tmpdir")
    p.set_defaults(func=command_scan)
    p = sub.add_parser("compile")
    p.add_argument("--out", required=True)
    p.set_defaults(func=command_compile)
    p = sub.add_parser("analyze")
    p.add_argument("--out", required=True)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--num-shards", type=int, default=16)
    p.add_argument("--shard-index", type=int, default=0)
    p.add_argument("--min-cells", type=int, default=50)
    p.add_argument("--min-open-cells", type=int, default=5)
    p.add_argument("--min-open-fraction", type=float, default=0.05)
    p.add_argument("--min-matched", type=int, default=1000)
    p.add_argument("--min-motif-hits", type=int, default=20)
    p.add_argument("--cell-batch-size", type=int, default=128)
    p.set_defaults(func=command_analyze)
    p = sub.add_parser("aggregate")
    p.add_argument("--out", required=True)
    p.add_argument("--q", type=float, default=0.05)
    p.add_argument("--min-joint", type=int, default=20)
    p.set_defaults(func=command_aggregate)
    p = sub.add_parser("combine")
    p.add_argument("--candidate-summary", required=True)
    p.add_argument("--catlas-results", help="Optional CATlas pair results; never a selection filter")
    p.add_argument("--corpus-summary", help="Optional corpus candidate summary; never a selection filter")
    p.add_argument("--include-status", action="append", default=["replicated_donor_supported"])
    p.add_argument("--out", required=True)
    p.set_defaults(func=command_combine)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
