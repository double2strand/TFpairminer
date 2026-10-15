#!/usr/bin/env python3
"""Integrate replicated TF pairs with CATlas cCRE accessibility and JASPAR motifs.

Stages: plan -> scan (parallel) -> compile -> analyze (parallel) -> aggregate.
The result is motif co-enrichment in cell-type-accessible cCREs, not TF occupancy.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import gzip
import hashlib
import json
import math
import os
import re
import sqlite3
import subprocess
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.stats import fisher_exact


MIRROR_REPO = "Genentech/binary-atac-tutorial-data"
MIRROR_REVISION = "6b6d749d445605f2cdb98e5ca83018b63271b3c4"
MIRROR_SHA256 = "bbe4343c308c83efee0c4ec085229c82a6a7077df1bd66cb4096801a4578aa60"
ONTOLOGY_URL = "https://decoder-genetics.wustl.edu/catlasv1/humanenhancer/data/Cell_ontology.tsv"
JASPAR_URL = "https://jaspar2026.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
ONTOLOGY_SHA256 = "7ef6945a19298bd24c12f42768df64d5d1374029b28755f551f7cd0f14ad1062"
JASPAR_SHA256 = "ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15"
PLAN_FIELDS = ["candidate_id", "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
               "cell_type_id", "cell_type", "tissue_id", "tissue", "disease_id", "disease",
               "validation_status", "ts_cell_r", "ts_lift", "ts_n_donors", "catlas_cluster",
               "catlas_stage", "catlas_cell_type_id", "mapping_tier", "mapping_reviewed",
               "mapping_status", "tf1_motifs", "tf2_motifs", "motif_status"]
RAW_FIELDS = [*PLAN_FIELDS, "n_accessible_distal", "n_matched_foreground", "n_matched_background",
              "matching_retention", "tf1_fg_hits", "tf1_bg_hits", "tf2_fg_hits", "tf2_bg_hits",
              "joint_fg_hits", "joint_bg_hits", "tf1_odds_ratio", "tf2_odds_ratio",
              "joint_odds_ratio", "tf1_p", "tf2_p", "joint_p", "joint_lift_fg",
              "joint_lift_bg", "analysis_status"]
FINAL_FIELDS = [*RAW_FIELDS, "tf1_q", "tf2_q", "joint_q", "catlas_result"]
CLUSTER_MAP_FIELDS = ["catlas_cluster", "cell_type_id", "tissue_id", "stage", "mapping_tier", "reviewed"]


def open_text(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8", newline="")


def read_tsv(path):
    with open_text(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".tmp.gz" if str(path).endswith(".gz") else ".tmp"
    fd, name = tempfile.mkstemp(prefix=path.stem + ".", suffix=suffix, dir=path.parent)
    os.close(fd)
    try:
        with open_text(name, "wt") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(name, path)
    except BaseException:
        Path(name).unlink(missing_ok=True)
        raise


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(value, fh, indent=2)
            fh.write("\n")
        os.replace(name, path)
    except BaseException:
        Path(name).unlink(missing_ok=True)
        raise


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_fetch(args):
    """Fetch and verify the pinned bulk inputs; the genome FASTA is supplied separately."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    sources = [
        ("catlas_binary_atac.h5ad",
         f"https://huggingface.co/datasets/{MIRROR_REPO}/resolve/{MIRROR_REVISION}/data.h5ad",
         MIRROR_SHA256),
        ("catlas_cell_ontology.tsv", ONTOLOGY_URL, ONTOLOGY_SHA256),
        ("jaspar2026_core_vertebrates.meme", JASPAR_URL, JASPAR_SHA256),
    ]
    fetched = []
    for filename, url, expected in sources:
        dest = out / filename
        if dest.exists():
            observed = sha256(dest)
            if observed != expected:
                raise RuntimeError(f"Existing file has incorrect SHA256: {dest} ({observed})")
        else:
            fd, tmpname = tempfile.mkstemp(prefix=filename + ".", suffix=".part", dir=out)
            try:
                with os.fdopen(fd, "wb") as target, urllib.request.urlopen(url, timeout=120) as source:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        target.write(block)
                observed = sha256(tmpname)
                if observed != expected:
                    raise RuntimeError(f"Checksum mismatch for {filename}: {observed}")
                os.replace(tmpname, dest)
            except BaseException:
                Path(tmpname).unlink(missing_ok=True)
                raise
        fetched.append({"file": filename, "url": url, "sha256": expected})
        print(f"Verified {dest}")
    write_json(out / "input_manifest.json", {"sources": fetched, "genome_build": "GRCh38/hg38"})


def norm(name):
    return re.sub(r"[^a-z0-9]", "", name.lower().replace("cholangiocytes", "cholangiocyte"))


def parse_meme(path):
    """Return MEME header and {motif_id: (display_name, complete_block)}."""
    header = []
    blocks = {}
    current = []
    motif_id = motif_name = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("MOTIF "):
                if motif_id:
                    blocks[motif_id] = (motif_name, "".join(current))
                parts = line.strip().split(maxsplit=2)
                motif_id = parts[1]
                motif_name = parts[2] if len(parts) > 2 else motif_id
                current = [line]
            elif motif_id:
                current.append(line)
            else:
                header.append(line)
    if motif_id:
        blocks[motif_id] = (motif_name, "".join(current))
    if not blocks or not "MEME version" in "".join(header):
        raise ValueError("Input is not a MEME motif file")
    return "".join(header), blocks


def motif_matches(symbol, blocks):
    exact = sorted(mid for mid, (name, _) in blocks.items() if name.upper() == symbol.upper())
    composite = sorted(mid for mid, (name, _) in blocks.items()
                       if "::" in name and symbol.upper() in {p.upper() for p in name.split("::")})
    return exact, composite


def ontology_index(path):
    rows = list(read_tsv(path))
    index = defaultdict(list)
    for row in rows:
        name = row.get("Cell type", "")
        for cell_id in row.get("Cell Ontology ID", "").split(","):
            cell_id = cell_id.strip()
            if cell_id.startswith("CL:"):
                index[norm(name)].append((cell_id, name))
    return index, rows


def cluster_catalog(h5ad, ontology):
    if "cell type" not in h5ad.obs or not {"chrom", "start", "end", "Class"}.issubset(h5ad.var.columns):
        raise ValueError("CATlas H5AD must have obs['cell type'] and var chrom/start/end/Class")
    if h5ad.obs["cell type"].astype(str).duplicated().any():
        raise ValueError("CATlas cluster names are not unique")
    index, _ = ontology_index(ontology)
    catalog = []
    for name in h5ad.obs["cell type"].astype(str):
        matches = sorted(set(cid for cid, _ in index.get(norm(name), [])))
        catalog.append({"catlas_cluster": name,
                        "cell_type_id": matches[0] if len(matches) == 1 else "",
                        "stage": "fetal" if name.lower().startswith("fetal ") else "adult",
                        "name_match_status": "exact_name" if len(matches) == 1 else "needs_review"})
    return catalog


def load_reviewed_map(path, catalog):
    if not path:
        return []
    rows = list(read_tsv(path))
    cluster_info = {r["catlas_cluster"]: r for r in catalog}
    seen = set()
    for row in rows:
        if row["catlas_cluster"] not in cluster_info:
            raise ValueError(f"Unknown CATlas cluster in reviewed map: {row['catlas_cluster']}")
        if row["reviewed"].lower() not in ("yes", "no"):
            raise ValueError("Mapping reviewed must be Yes or No")
        if row["mapping_tier"] not in ("exact_reviewed", "proxy_fetal", "proxy_related_celltype"):
            raise ValueError("Invalid mapping_tier")
        if row["stage"] not in ("adult", "fetal"):
            raise ValueError("Stage must be adult or fetal")
        if row["stage"] != cluster_info[row["catlas_cluster"]]["stage"]:
            raise ValueError(f"Stage differs from CATlas cluster: {row['catlas_cluster']}")
        if row["mapping_tier"] == "proxy_fetal" and row["stage"] != "fetal":
            raise ValueError("proxy_fetal mappings must point to fetal clusters")
        if row["mapping_tier"] == "exact_reviewed" and not row["tissue_id"]:
            raise ValueError("exact_reviewed mappings require an explicit tissue_id")
        if (row["mapping_tier"] == "exact_reviewed" and
                cluster_info[row["catlas_cluster"]]["cell_type_id"] and
                cluster_info[row["catlas_cluster"]]["cell_type_id"] != row["cell_type_id"]):
            raise ValueError(f"Reviewed exact cell-type ID conflicts with CATlas ontology: {row['catlas_cluster']}")
        key = (row["catlas_cluster"], row["cell_type_id"], row["tissue_id"])
        if key in seen:
            raise ValueError(f"Duplicate crosswalk mapping: {key}")
        seen.add(key)
    return rows


def make_candidate_mappings(candidate, catalog, reviewed_map):
    wanted = candidate["cell_type_id"]
    tissue = candidate["tissue_id"]
    healthy = (candidate.get("disease_id", "") in ("", "PATO:0000461") and
               candidate.get("disease", "").lower() in ("", "normal", "healthy"))
    reviewed = [r for r in reviewed_map if r["cell_type_id"] == wanted and
                (not r["tissue_id"] or r["tissue_id"] == tissue) and r["reviewed"].lower() == "yes"]
    if reviewed:
        if not healthy and any(r["mapping_tier"] == "exact_reviewed" for r in reviewed):
            raise ValueError("A disease-context candidate cannot have an exact_reviewed mapping to healthy CATlas")
        return [{"catlas_cluster": r["catlas_cluster"], "catlas_stage": r["stage"],
                 "catlas_cell_type_id": r["cell_type_id"], "mapping_tier": r["mapping_tier"],
                 "mapping_reviewed": "Yes", "mapping_status": "mapped"} for r in reviewed]
    auto = [r for r in catalog if r["cell_type_id"] == wanted]
    if auto:
        return [{"catlas_cluster": r["catlas_cluster"], "catlas_stage": r["stage"],
                 "catlas_cell_type_id": wanted,
                 "mapping_tier": "ontology_exact_unreviewed" if healthy else "proxy_disease_context",
                 "mapping_reviewed": "No", "mapping_status": "mapped_provisional"} for r in auto]
    return [{"catlas_cluster": "", "catlas_stage": "", "catlas_cell_type_id": "",
             "mapping_tier": "", "mapping_reviewed": "No", "mapping_status": "no_exact_celltype"}]


def command_plan(args):
    out = Path(args.out)
    if (out / "manifest.json").exists():
        raise SystemExit("Plan already exists; choose a new output directory")
    out.mkdir(parents=True, exist_ok=True)
    candidates = [r for r in read_tsv(args.candidate_summary) if r.get("validation_status") in set(args.include_status)]
    if not candidates:
        raise ValueError("No selected candidates in Tabula Sapiens summary")
    ids = [r["candidate_id"] for r in candidates]
    if len(ids) != len(set(ids)):
        raise ValueError("Candidate IDs are not unique")
    header, blocks = parse_meme(args.meme)
    h5ad = ad.read_h5ad(args.catlas_h5ad, backed="r")
    try:
        catalog = cluster_catalog(h5ad, args.ontology)
        n_peaks = h5ad.n_vars
        n_distal_peaks = int(h5ad.var["Class"].astype(str).str.lower().eq("distal").sum())
        if not all(str(x).startswith("chr") for x in h5ad.var["chrom"].iloc[:1000]):
            raise ValueError("Expected hg38 chr-prefixed CATlas coordinates")
    finally:
        h5ad.file.close()
    reviewed = load_reviewed_map(args.cluster_map, catalog)
    _, ontology_rows = ontology_index(args.ontology)
    official_ids = {cell_id.strip() for r in ontology_rows for cell_id in r.get("Cell Ontology ID", "").split(",")}
    genes = {(r["tf1_ensembl"], r["tf1_symbol"]) for r in candidates} | {(r["tf2_ensembl"], r["tf2_symbol"]) for r in candidates}
    motif_rows = []
    gene_motifs = {}
    selected = set()
    for gid, symbol in sorted(genes):
        exact, composite = motif_matches(symbol, blocks)
        selected.update(exact)
        gene_motifs[gid] = set(exact)
        motif_rows.append({"ensembl_id": gid, "symbol": symbol, "motif_ids": ",".join(exact),
                           "composite_ids": ",".join(composite),
                           "mapping_status": "exact_single_tf" if exact else ("composite_only" if composite else "no_exact_motif")})
    write_tsv(out / "motif_map.tsv", motif_rows, ["ensembl_id", "symbol", "motif_ids", "composite_ids", "mapping_status"])
    write_tsv(out / "clusters.tsv", catalog, ["catlas_cluster", "cell_type_id", "stage", "name_match_status"])
    suggestions = []
    names = [r["catlas_cluster"] for r in catalog]
    for context in sorted({(r["cell_type_id"], r["cell_type"], r["tissue_id"], r["tissue"]) for r in candidates}):
        cell_id, cell_name, tissue_id, tissue_name = context
        exact = [r["catlas_cluster"] for r in catalog if r["cell_type_id"] == cell_id]
        near = difflib.get_close_matches(cell_name, names, n=3, cutoff=0.3)
        suggestions.append({"cell_type_id": cell_id, "cell_type": cell_name, "tissue_id": tissue_id,
                            "tissue": tissue_name, "exact_ontology_clusters": ";".join(exact),
                            "name_suggestions": ";".join(near)})
    write_tsv(out / "mapping_review.tsv", suggestions,
              ["cell_type_id", "cell_type", "tissue_id", "tissue", "exact_ontology_clusters", "name_suggestions"])
    rows = []
    for candidate in candidates:
        a = gene_motifs[candidate["tf1_ensembl"]]
        b = gene_motifs[candidate["tf2_ensembl"]]
        motif_status = "both_exact_motifs" if a and b and not (a & b) else (
            "shared_motif_ambiguous" if a & b else "motif_unavailable")
        for mapping in make_candidate_mappings(candidate, catalog, reviewed):
            if mapping["mapping_status"] == "no_exact_celltype" and candidate["cell_type_id"] in official_ids:
                mapping["mapping_status"] = "ontology_present_cluster_unresolved"
            rows.append({**candidate, **mapping, "tf1_motifs": ",".join(sorted(a)),
                         "tf2_motifs": ",".join(sorted(b)), "motif_status": motif_status})
    write_tsv(out / "plan.tsv", rows, PLAN_FIELDS)
    with open(out / "selected_motifs.meme", "w", encoding="utf-8") as fh:
        fh.write(header)
        for mid in sorted(selected):
            fh.write(blocks[mid][1])
    manifest = {"catlas_source": "CATlas-derived binary-atac-tutorial-data mirror", "catlas_mirror_repo": MIRROR_REPO,
                "catlas_mirror_revision": MIRROR_REVISION, "genome_build": "GRCh38/hg38",
                "candidate_summary_sha256": sha256(args.candidate_summary), "catlas_h5ad_sha256": sha256(args.catlas_h5ad),
                "ontology_sha256": sha256(args.ontology), "jaspar_meme_sha256": sha256(args.meme),
                "cluster_map_sha256": sha256(args.cluster_map) if args.cluster_map else "",
                "n_candidates": len(candidates), "n_plan_rows": len(rows), "n_clusters": len(catalog),
                "n_peaks": n_peaks, "n_distal_peaks": n_distal_peaks,
                "n_selected_motifs": len(selected),
                "n_mapped_rows": sum(bool(r["catlas_cluster"]) for r in rows),
                "n_unmapped_rows": sum(not r["catlas_cluster"] for r in rows),
                "included_validation_statuses": sorted(set(args.include_status)),
                "motif_interpretation": "JASPAR sequence-model matches; not direct TF occupancy or biochemical pairing"}
    write_json(out / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


def load_plan(out):
    out = Path(out)
    return json.loads((out / "manifest.json").read_text()), list(read_tsv(out / "plan.tsv"))


def distal_indices(h5ad):
    return np.flatnonzero(h5ad.var["Class"].astype(str).str.lower().eq("distal").to_numpy())


def fasta_gc(path):
    name = None
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        yield name, "".join(seq)


def command_scan(args):
    out = Path(args.out)
    manifest, _ = load_plan(out)
    if manifest["n_mapped_rows"] == 0:
        raise RuntimeError("No mapped CATlas clusters; aggregate the coverage-only result without scanning")
    if args.num_shards < 1 or not 0 <= args.shard_index < args.num_shards:
        raise ValueError("Invalid scan shard index/count")
    if manifest["n_selected_motifs"] == 0:
        raise ValueError("No direct single-TF motifs were mapped")
    if not 0 < args.fimo_p <= 1:
        raise ValueError("Invalid FIMO p-value threshold")
    if args.max_stored_scores < 100000:
        raise ValueError("FIMO max-stored-scores must be at least 100000")
    marker = out / "scan" / f"{args.shard_index:04d}.done.json"
    if marker.exists():
        old = json.loads(marker.read_text())
        if (old["num_shards"] != args.num_shards or old["fimo_p"] != args.fimo_p or
                old["max_stored_scores"] != args.max_stored_scores):
            raise RuntimeError("Completed scan shard has different parameters")
        print(f"Scan shard {args.shard_index} already complete")
        return
    h5ad = ad.read_h5ad(args.catlas_h5ad, backed="r")
    try:
        if h5ad.n_vars != manifest["n_peaks"]:
            raise RuntimeError("CATlas peak count differs from plan")
        idx = distal_indices(h5ad)
        idx = idx[idx % args.num_shards == args.shard_index]
        rows = h5ad.var.iloc[idx][["chrom", "start", "end"]].itertuples(index=False, name=None)
        with tempfile.TemporaryDirectory(prefix="catlas_scan_", dir=args.tmpdir) as tmp:
            bed = Path(tmp) / "peaks.bed"
            fa = Path(tmp) / "peaks.fa"
            with open(bed, "w", encoding="utf-8") as fh:
                for peak_id, (chrom, start, end) in zip(idx, rows):
                    fh.write(f"{chrom}\t{int(start)}\t{int(end)}\tpeak_{peak_id}\n")
            with open(fa, "w", encoding="utf-8") as fh:
                subprocess.run([args.bedtools, "getfasta", "-fi", args.genome_fasta, "-bed", str(bed), "-nameOnly"],
                               stdout=fh, check=True)
            gc_rows = []
            for name, seq in fasta_gc(fa):
                if not name.startswith("peak_"):
                    raise RuntimeError(f"Unexpected FASTA name: {name}")
                valid = sum(x in "ACGT" for x in seq)
                gc = sum(x in "GC" for x in seq) / valid if valid / max(len(seq), 1) >= 0.8 else ""
                gc_rows.append({"peak_index": int(name[5:]), "gc": gc})
            if len(gc_rows) != len(idx):
                raise RuntimeError(f"FASTA extraction incomplete: {len(gc_rows)} of {len(idx)} peaks")
            # SQLite deduplicates multiple hits of the same motif in one cCRE while bounding RAM.
            conn = sqlite3.connect(str(Path(tmp) / "hits.sqlite"))
            conn.execute("CREATE TABLE hits (motif_id TEXT, peak_index INTEGER, PRIMARY KEY(motif_id,peak_index)) WITHOUT ROWID")
            cmd = [args.fimo, "--text", "--no-qvalue", "--thresh", str(args.fimo_p),
                   "--max-stored-scores", str(args.max_stored_scores),
                   str(out / "selected_motifs.meme"), str(fa)]
            stderr_path = Path(tmp) / "fimo.stderr"
            stderr_fh = open(stderr_path, "w", encoding="utf-8")
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
                return_code = proc.wait()
                err = stderr_path.read_text()
                if return_code != 0:
                    stderr_fh.close()
                    raise RuntimeError(f"FIMO failed: {err[-2000:]}")
                if re.search(r"max(?:imum)?[- ]stored[- ]scores|discard(?:ing|ed).*scores", err, re.I):
                    raise RuntimeError("FIMO reported possible score truncation; increase --max-stored-scores")
            except BaseException:
                proc.kill()
                proc.wait()
                raise
            finally:
                stderr_fh.close()
            write_tsv(out / "scan" / f"{args.shard_index:04d}.hits.tsv.gz",
                      ({"motif_id": mid, "peak_index": pid} for mid, pid in conn.execute("SELECT motif_id,peak_index FROM hits")),
                      ["motif_id", "peak_index"])
            conn.close()
            write_tsv(out / "scan" / f"{args.shard_index:04d}.gc.tsv.gz", gc_rows, ["peak_index", "gc"])
    finally:
        h5ad.file.close()
    write_json(marker, {"shard_index": args.shard_index, "num_shards": args.num_shards,
                        "fimo_p": args.fimo_p, "max_stored_scores": args.max_stored_scores,
                        "n_distal_peaks": len(idx), "n_fimo_sites": n_sites})
    print(f"Scan shard {args.shard_index}: {len(idx)} distal peaks; {n_sites} FIMO sites")


def command_compile(args):
    out = Path(args.out)
    manifest, _ = load_plan(out)
    markers = sorted((out / "scan").glob("*.done.json"))
    if not markers:
        raise RuntimeError("No completed scan shards")
    configs = [json.loads(p.read_text()) for p in markers]
    n_shards = configs[0]["num_shards"]
    if len(configs) != n_shards or {x["shard_index"] for x in configs} != set(range(n_shards)):
        raise RuntimeError("Scan shards incomplete")
    if len({(x["fimo_p"], x["max_stored_scores"]) for x in configs}) != 1:
        raise RuntimeError("Scan shards used different FIMO parameters")
    if sum(x["n_distal_peaks"] for x in configs) != manifest["n_distal_peaks"]:
        raise RuntimeError("Scanned distal peak count differs from the plan")
    if (out / "compiled" / "hits.npy").exists():
        raise RuntimeError("Compiled hits already exist")
    _, blocks = parse_meme(out / "selected_motifs.meme")
    motif_ids = sorted(blocks)
    motif_index = {m: i for i, m in enumerate(motif_ids)}
    n_peaks = manifest["n_peaks"]
    compiled = out / "compiled"
    compiled.mkdir(parents=True, exist_ok=True)
    hits = np.lib.format.open_memmap(compiled / "hits.npy", mode="w+", dtype=np.uint8,
                                     shape=(len(motif_ids), (n_peaks + 7) // 8))
    hits[:] = 0
    gc = np.lib.format.open_memmap(compiled / "gc.npy", mode="w+", dtype=np.float32, shape=(n_peaks,))
    gc[:] = np.nan
    n_hits = 0
    n_gc = 0
    for s in range(n_shards):
        for row in read_tsv(out / "scan" / f"{s:04d}.hits.tsv.gz"):
            mid, peak = row["motif_id"], int(row["peak_index"])
            if mid not in motif_index or not 0 <= peak < n_peaks:
                raise RuntimeError("Unexpected motif or peak in scan output")
            hits[motif_index[mid], peak // 8] |= np.uint8(1 << (peak % 8))
            n_hits += 1
        for row in read_tsv(out / "scan" / f"{s:04d}.gc.tsv.gz"):
            peak = int(row["peak_index"])
            if row["gc"]:
                gc[peak] = float(row["gc"])
                n_gc += 1
    hits.flush()
    gc.flush()
    write_tsv(compiled / "motif_ids.tsv", ({"motif_id": mid, "row_index": i} for i, mid in enumerate(motif_ids)),
              ["motif_id", "row_index"])
    write_json(compiled / "manifest.json", {"n_peaks": n_peaks, "n_motifs": len(motif_ids),
                "n_unique_motif_peak_hits": n_hits, "n_peaks_with_valid_gc": n_gc,
                "scan_shards": n_shards, "fimo_p": configs[0]["fimo_p"],
                "max_stored_scores": configs[0]["max_stored_scores"]})
    print(f"Compiled {n_hits} motif-peak hits and GC for {n_gc} peaks")


def matched_indices(accessible, eligible, gc, width, seed):
    rng = np.random.default_rng(seed)
    width_bins = np.searchsorted(np.quantile(width[eligible], [0.25, 0.5, 0.75]), width, side="right")
    gc_bins = np.minimum(np.floor(np.nan_to_num(gc, nan=0) * 10).astype(np.int16), 9)
    fg_parts, bg_parts = [], []
    for gb in range(10):
        for wb in range(4):
            cell = eligible & (gc_bins == gb) & (width_bins == wb)
            fg = np.flatnonzero(cell & accessible)
            bg = np.flatnonzero(cell & ~accessible)
            n = min(len(fg), len(bg))
            if n:
                fg_parts.append(rng.choice(fg, n, replace=False))
                bg_parts.append(rng.choice(bg, n, replace=False))
    return (np.concatenate(fg_parts) if fg_parts else np.empty(0, dtype=np.int64),
            np.concatenate(bg_parts) if bg_parts else np.empty(0, dtype=np.int64))


def odds_p(hits_fg, n_fg, hits_bg, n_bg):
    odds, p = fisher_exact([[hits_fg, n_fg - hits_fg], [hits_bg, n_bg - hits_bg]], alternative="greater")
    return float(odds), float(p)


def motif_presence(motif_ids, packed_hits, motif_index, n_peaks):
    result = np.zeros(n_peaks, dtype=bool)
    for mid in motif_ids.split(","):
        if mid:
            if mid not in motif_index:
                raise RuntimeError(f"Missing compiled motif: {mid}")
            result |= np.unpackbits(packed_hits[motif_index[mid]], bitorder="little")[:n_peaks].astype(bool)
    return result


def command_analyze(args):
    out = Path(args.out)
    manifest, plan = load_plan(out)
    if args.num_shards < 1 or not 0 <= args.shard_index < args.num_shards:
        raise ValueError("Invalid analysis shard index/count")
    if not (out / "compiled" / "manifest.json").exists():
        raise RuntimeError("Compile the motif scans first")
    groups = defaultdict(list)
    for row in plan:
        if row["catlas_cluster"]:
            groups[row["catlas_cluster"]].append(row)
    clusters = sorted(c for c in groups if int(hashlib.sha256(c.encode()).hexdigest(), 16) % args.num_shards == args.shard_index)
    marker = out / "analysis" / f"{args.shard_index:04d}.done.json"
    if marker.exists():
        old = json.loads(marker.read_text())
        if (old["num_shards"] != args.num_shards or old["min_matched"] != args.min_matched or
                old["min_motif_hits"] != args.min_motif_hits):
            raise RuntimeError("Completed analysis shard has different parameters")
        print(f"Analysis shard {args.shard_index} already complete")
        return
    h5ad = ad.read_h5ad(args.catlas_h5ad, backed="r")
    try:
        if h5ad.n_vars != manifest["n_peaks"]:
            raise RuntimeError("CATlas H5AD differs from plan")
        packed = np.load(out / "compiled" / "hits.npy", mmap_mode="r")
        gc = np.load(out / "compiled" / "gc.npy", mmap_mode="r")
        motif_index = {r["motif_id"]: int(r["row_index"]) for r in read_tsv(out / "compiled" / "motif_ids.tsv")}
        names = list(h5ad.obs["cell type"].astype(str))
        name_to_index = {name: i for i, name in enumerate(names)}
        width = (h5ad.var["end"].to_numpy() - h5ad.var["start"].to_numpy()).astype(np.int32)
        distal = h5ad.var["Class"].astype(str).str.lower().eq("distal").to_numpy()
        eligible = distal & np.isfinite(gc) & (width > 0)
        results = []
        for cluster in clusters:
            x = h5ad.X[name_to_index[cluster], :]
            accessible = np.asarray(x.toarray() if hasattr(x, "toarray") else x).ravel() > 0
            n_accessible = int(np.count_nonzero(accessible & eligible))
            seed = int(hashlib.sha256((cluster + "|CATlas").encode()).hexdigest()[:16], 16)
            fg, bg = matched_indices(accessible, eligible, gc, width, seed)
            retention = len(fg) / n_accessible if n_accessible else 0.0
            cache = {}
            for row in groups[cluster]:
                base = {**row, "n_accessible_distal": n_accessible,
                        "n_matched_foreground": len(fg), "n_matched_background": len(bg),
                        "matching_retention": retention}
                if len(fg) < args.min_matched:
                    results.append({**base, "analysis_status": "insufficient_matched_cCREs"})
                    continue
                if row["motif_status"] != "both_exact_motifs":
                    results.append({**base, "analysis_status": row["motif_status"]})
                    continue
                for gene_id, motif_field in [(row["tf1_ensembl"], row["tf1_motifs"]),
                                             (row["tf2_ensembl"], row["tf2_motifs"])]:
                    if gene_id not in cache:
                        cache[gene_id] = motif_presence(motif_field, packed, motif_index, manifest["n_peaks"])
                a = cache[row["tf1_ensembl"]]
                b = cache[row["tf2_ensembl"]]
                af, ab = int(a[fg].sum()), int(a[bg].sum())
                bf, bb = int(b[fg].sum()), int(b[bg].sum())
                both = a & b
                jf, jb = int(both[fg].sum()), int(both[bg].sum())
                if min(af, bf) < args.min_motif_hits:
                    results.append({**base, "tf1_fg_hits": af, "tf1_bg_hits": ab,
                                    "tf2_fg_hits": bf, "tf2_bg_hits": bb, "joint_fg_hits": jf,
                                    "joint_bg_hits": jb, "analysis_status": "inconclusive_sparse_motif"})
                    continue
                or_a, p_a = odds_p(af, len(fg), ab, len(bg))
                or_b, p_b = odds_p(bf, len(fg), bb, len(bg))
                or_joint, p_joint = odds_p(jf, len(fg), jb, len(bg))
                lift_fg = jf * len(fg) / (af * bf) if af and bf else 0.0
                lift_bg = jb * len(bg) / (ab * bb) if ab and bb else 0.0
                results.append({**base, "tf1_fg_hits": af, "tf1_bg_hits": ab,
                                "tf2_fg_hits": bf, "tf2_bg_hits": bb, "joint_fg_hits": jf,
                                "joint_bg_hits": jb, "tf1_odds_ratio": or_a, "tf2_odds_ratio": or_b,
                                "joint_odds_ratio": or_joint, "tf1_p": p_a, "tf2_p": p_b,
                                "joint_p": p_joint, "joint_lift_fg": lift_fg, "joint_lift_bg": lift_bg,
                                "analysis_status": "tested"})
    finally:
        h5ad.file.close()
    write_tsv(out / "analysis" / f"{args.shard_index:04d}.tsv.gz", results, RAW_FIELDS)
    write_json(marker, {"shard_index": args.shard_index, "num_shards": args.num_shards,
                        "min_matched": args.min_matched, "min_motif_hits": args.min_motif_hits,
                        "n_clusters": len(clusters), "n_rows": len(results)})
    print(f"Analysis shard {args.shard_index}: {len(clusters)} CATlas clusters; {len(results)} pair rows")


def bh(pvalues):
    if not pvalues:
        return []
    n = len(pvalues)
    order = np.argsort(np.asarray(pvalues, dtype=float), kind="stable")
    out = np.ones(n, dtype=float)
    minimum = 1.0
    for rank in range(n - 1, -1, -1):
        i = int(order[rank])
        minimum = min(minimum, float(pvalues[i]) * n / (rank + 1))
        out[i] = minimum
    return out.tolist()


def command_aggregate(args):
    out = Path(args.out)
    manifest, plan = load_plan(out)
    if not 0 < args.q <= 1 or args.min_joint < 1:
        raise ValueError("Invalid FDR or minimum joint-peak threshold")
    results = []
    if manifest["n_mapped_rows"]:
        markers = sorted((out / "analysis").glob("*.done.json"))
        if not markers:
            raise RuntimeError("No analysis shards complete")
        configs = [json.loads(p.read_text()) for p in markers]
        n_shards = configs[0]["num_shards"]
        if len(configs) != n_shards or {x["shard_index"] for x in configs} != set(range(n_shards)):
            raise RuntimeError("Analysis shards incomplete")
        if len({(x["min_matched"], x["min_motif_hits"]) for x in configs}) != 1:
            raise RuntimeError("Analysis shards used different parameters")
        for s in range(n_shards):
            results.extend(read_tsv(out / "analysis" / f"{s:04d}.tsv.gz"))
    else:
        n_shards = 0
    results.extend({**r, "analysis_status": r["mapping_status"]} for r in plan if not r["catlas_cluster"])
    expected = {(r["candidate_id"], r["catlas_cluster"]) for r in plan}
    observed = [(r["candidate_id"], r["catlas_cluster"]) for r in results]
    if len(expected) != len(plan) or len(observed) != len(set(observed)) or set(observed) != expected:
        raise RuntimeError("Analysis rows do not match the planned candidate-by-cluster coverage")
    tested = [r for r in results if r["analysis_status"] == "tested"]
    gene_p = {}
    for r in tested:
        for gene, p in [(r["tf1_ensembl"], r["tf1_p"]), (r["tf2_ensembl"], r["tf2_p"])]:
            key = (r["catlas_cluster"], gene)
            if key in gene_p and abs(gene_p[key] - float(p)) > 1e-10:
                raise RuntimeError("Inconsistent gene-level p-value in one cluster")
            gene_p[key] = float(p)
    keys = sorted(gene_p)
    gene_q = dict(zip(keys, bh([gene_p[k] for k in keys])))
    pair_q = bh([float(r["joint_p"]) for r in tested])
    for r, q in zip(tested, pair_q):
        r["tf1_q"] = gene_q[(r["catlas_cluster"], r["tf1_ensembl"])]
        r["tf2_q"] = gene_q[(r["catlas_cluster"], r["tf2_ensembl"])]
        r["joint_q"] = q
    counts = Counter()
    for r in results:
        if r["analysis_status"] != "tested":
            status = r["analysis_status"]
        else:
            pass_screen = (float(r["tf1_q"]) <= args.q and float(r["tf2_q"]) <= args.q and
                           float(r["joint_q"]) <= args.q and
                           float(r["tf1_odds_ratio"]) > 1 and float(r["tf2_odds_ratio"]) > 1 and
                           float(r["joint_odds_ratio"]) > 1 and int(r["joint_fg_hits"]) >= args.min_joint and
                           float(r["joint_lift_fg"]) > float(r["joint_lift_bg"]))
            if r["mapping_tier"] == "exact_reviewed":
                status = "catlas_coenriched" if pass_screen else "not_coenriched"
            elif r["mapping_tier"].startswith("proxy_"):
                status = "proxy_coenriched" if pass_screen else "proxy_not_coenriched"
            else:
                status = "provisional_coenriched" if pass_screen else "provisional_not_coenriched"
        r["catlas_result"] = status
        counts[status] += 1
    agg = out / "aggregate"
    if (agg / "catlas_pair_results.tsv.gz").exists():
        raise RuntimeError("Aggregate already exists; use a new plan directory or move old aggregate aside")
    write_tsv(agg / "catlas_pair_results.tsv.gz", results, FINAL_FIELDS)
    write_json(agg / "summary.json", {"n_input_candidates": len({r["candidate_id"] for r in plan}),
                "n_plan_rows": len(plan), "n_result_rows": len(results), "n_tested_rows": len(tested),
                "n_gene_tests": len(gene_q), "n_joint_tests": len(pair_q),
                "n_analysis_shards": n_shards, "result_counts": dict(counts),
                "q_threshold": args.q, "min_joint_peaks": args.min_joint,
                "interpretation": "Motif-sequence co-enrichment in accessible cCREs, not direct TF binding"})
    print(json.dumps(dict(counts), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--out", required=True)
    f.set_defaults(func=command_fetch)
    p = sub.add_parser("plan")
    p.add_argument("--candidate-summary", required=True)
    p.add_argument("--catlas-h5ad", required=True)
    p.add_argument("--ontology", required=True)
    p.add_argument("--meme", required=True)
    p.add_argument("--cluster-map", help="Optional reviewed CATlas-cluster to RNA-cell-type crosswalk")
    p.add_argument("--include-status", action="append", default=["replicated_donor_supported"])
    p.add_argument("--out", required=True)
    p.set_defaults(func=command_plan)
    s = sub.add_parser("scan")
    s.add_argument("--out", required=True)
    s.add_argument("--catlas-h5ad", required=True)
    s.add_argument("--genome-fasta", required=True)
    s.add_argument("--fimo", default="fimo")
    s.add_argument("--bedtools", default="bedtools")
    s.add_argument("--fimo-p", type=float, default=1e-4)
    s.add_argument("--max-stored-scores", type=int, default=5000000)
    s.add_argument("--num-shards", type=int, default=1)
    s.add_argument("--shard-index", type=int, default=0)
    s.add_argument("--tmpdir")
    s.set_defaults(func=command_scan)
    c = sub.add_parser("compile")
    c.add_argument("--out", required=True)
    c.set_defaults(func=command_compile)
    a = sub.add_parser("analyze")
    a.add_argument("--out", required=True)
    a.add_argument("--catlas-h5ad", required=True)
    a.add_argument("--num-shards", type=int, default=1)
    a.add_argument("--shard-index", type=int, default=0)
    a.add_argument("--min-matched", type=int, default=1000)
    a.add_argument("--min-motif-hits", type=int, default=20)
    a.set_defaults(func=command_analyze)
    g = sub.add_parser("aggregate")
    g.add_argument("--out", required=True)
    g.add_argument("--q", type=float, default=0.05)
    g.add_argument("--min-joint", type=int, default=20)
    g.set_defaults(func=command_aggregate)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
