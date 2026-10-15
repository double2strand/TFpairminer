#!/usr/bin/env python3
"""Remote, resumable CELLxGENE Census TF-expression discovery.

Commands: plan (network), run (network, sharded), aggregate (local).
Outputs are tab-delimited gzip files so they can be inspected without a
specialized file format or a database server.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import sqlite3
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy import stats as sstats


OBS_COLUMNS = [
    "soma_joinid", "dataset_id", "cell_type_ontology_term_id", "cell_type",
    "tissue_general_ontology_term_id", "tissue_general",
    "disease_ontology_term_id", "disease", "assay_ontology_term_id", "assay",
    "donor_id", "raw_sum", "is_primary_data",
]
KEY_COLUMNS = ["cell_type_id", "cell_type", "tissue_id", "tissue", "disease_id", "disease", "assay_id", "assay"]
PROFILE_COLUMNS = ["dataset_id", "collection_id", *KEY_COLUMNS, "ensembl_id", "symbol", "hgnc_id", "category", "n_cells", "n_donors", "n_detected", "detection_fraction", "mean_log1p_cptt", "mean_raw_count"]
PAIR_COLUMNS = ["dataset_id", "collection_id", *KEY_COLUMNS, "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol", "n_cells", "n_donors", "tf1_detection_fraction", "tf2_detection_fraction", "co_detected_cells", "co_detection_lift", "cell_pearson_r", "donor_pearson_r",
               "resid_cell_r", "resid_cell_df", "metacell_r", "metacell_df", "n_metacells",
               "pseudobulk_r", "pseudobulk_df", "primary_stat", "primary_r", "primary_df", "p_value", "q_value_bh",
               "n_tested_pairs", "passed"]
SCHEMA_VERSION = 2
POSITIVE_CONTROLS = [("MYC", "MAX"), ("JUN", "FOS"), ("JUNB", "FOS"), ("JUN", "ATF3"), ("MAX", "MXD4"),
                     ("HIF1A", "ARNT"), ("CLOCK", "BMAL1"), ("STAT1", "STAT2"), ("SMAD3", "SMAD4"),
                     ("RELA", "NFKB1"), ("TCF3", "TCF12"), ("HNF1A", "HNF4A"), ("FOXA2", "HNF4A"),
                     ("SPI1", "CEBPA"), ("GATA1", "TAL1"), ("PBX1", "MEIS1"), ("SOX2", "POU5F1")]
AVAIL_COLUMNS = ["dataset_id", "ensembl_id", "symbol", "census_feature_present"]
ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def census_module():
    try:
        import cellxgene_census
    except ImportError as exc:
        raise SystemExit("Install requirements.txt before plan/run: python -m pip install -r requirements.txt") from exc
    return cellxgene_census


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


@contextmanager
def tsv_writer(path, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    os.close(fd)
    tmp = Path(name)
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(tmp, mode="wt", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            yield writer
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    else:
        os.replace(tmp, path)


def write_tsv(path, rows, fields):
    with tsv_writer(path, fields) as writer:
        writer.writerows(rows)


def open_tsv(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", newline="")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json_atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    tmp = Path(name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
            fh.write("\n")
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def load_register(path, tf_set):
    fields = {"strict": "Strict_TF_set", "exploratory": "Exploratory_TF_set", "expanded": "Expanded_regulator_set"}
    genes = {}
    no_ensembl = []
    for row in read_tsv(path):
        gid = row.get("Ensembl_gene_ID", "").split(".")[0]
        if not gid.startswith("ENSG"):
            no_ensembl.append({"ensembl_id": "", "symbol": row["Symbol"], "hgnc_id": row["HGNC_ID"],
                               "category": row["Primary_category"], "is_tf": row[fields[tf_set]] == "Yes",
                               "reason": "No Ensembl gene ID in compendium"})
            continue
        if gid in genes:
            raise ValueError(f"Duplicate Ensembl ID in register: {gid}")
        genes[gid] = {"ensembl_id": gid, "symbol": row["Symbol"], "hgnc_id": row["HGNC_ID"],
                      "category": row["Primary_category"], "is_tf": row[fields[tf_set]] == "Yes"}
    if not genes or not any(g["is_tf"] for g in genes.values()):
        raise ValueError("No usable TFs found in register")
    return genes, no_ensembl


def clean(value):
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def dataset_path(outdir, dataset_id):
    if not ID_RE.fullmatch(dataset_id):
        raise ValueError(f"Unexpected dataset ID: {dataset_id!r}")
    return Path(outdir) / "datasets" / dataset_id


def command_plan(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "manifest.json").exists():
        raise SystemExit("Plan already exists. Choose a new --out directory; this preserves release provenance.")
    genes, no_ensembl = load_register(args.register, args.tf_set)
    census_api = census_module()
    with census_api.open_soma(census_version=args.census_version) as census:
        summary = census["census_info"]["summary"].read().concat().to_pandas()
        matches = summary.loc[summary["label"] == "census_build_date", "value"]
        if matches.empty:
            raise RuntimeError("Census build date not found; refusing an unpinned run")
        release = str(matches.iloc[0])[:10]
        var = census["census_data"]["homo_sapiens"].ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()
        var = var[var["feature_id"].isin(genes)].copy()
        if var["feature_id"].duplicated().any():
            raise RuntimeError("Census contains duplicate feature_id values")
        mapped = set(var["feature_id"])
        if not mapped:
            raise RuntimeError("No compendium genes mapped to this Census release")
        feature_rows = [{**genes[r.feature_id], "soma_joinid": int(r.soma_joinid), "census_name": clean(r.feature_name)} for r in var.itertuples(index=False)]
        data = census["census_info"]["datasets"].read(column_names=["soma_joinid", "dataset_id", "collection_id", "collection_name", "dataset_title", "dataset_total_cell_count"]).concat().to_pandas()
        if args.dataset_id:
            data = data[data["dataset_id"].isin(args.dataset_id)]
            missing = sorted(set(args.dataset_id) - set(data["dataset_id"]))
            if missing:
                raise ValueError(f"Dataset IDs absent from Census: {missing}")
        if args.exclude_collection_id:
            data = data[~data["collection_id"].isin(args.exclude_collection_id)]
        if args.exclude_collection_name:
            data = data[~data["collection_name"].isin(args.exclude_collection_name)]
        data = data.sort_values("dataset_id")
        if args.max_datasets:
            data = data.head(args.max_datasets)
        dataset_rows = [{k: clean(v) for k, v in r.items()} for r in data.to_dict("records")]
    write_tsv(out / "features.tsv", sorted(feature_rows, key=lambda r: r["ensembl_id"]), ["ensembl_id", "symbol", "hgnc_id", "category", "is_tf", "soma_joinid", "census_name"])
    unmapped = [{**genes[g], "reason": "Not in pinned Census human RNA features"} for g in sorted(genes.keys() - mapped)]
    write_tsv(out / "unmapped.tsv", [*unmapped, *no_ensembl], ["ensembl_id", "symbol", "hgnc_id", "category", "is_tf", "reason"])
    write_tsv(out / "datasets.tsv", dataset_rows, ["soma_joinid", "dataset_id", "collection_id", "collection_name", "dataset_title", "dataset_total_cell_count"])
    manifest = {"requested_census_version": args.census_version, "resolved_census_version": release,
                "register_sha256": sha256(args.register), "register_path": str(Path(args.register).resolve()),
                "tf_set": args.tf_set, "n_register_genes": len(genes) + len(no_ensembl),
                "n_with_ensembl_id": len(genes), "n_without_ensembl_id": len(no_ensembl),
                "n_mapped_genes": len(mapped),
                "n_tf_mapped": sum(r["is_tf"] for r in feature_rows), "n_datasets": len(dataset_rows),
                "dataset_selection": "explicit_ids" if args.dataset_id else ("first_n" if args.max_datasets else "all"),
                "normalization": "log1p(10000 * raw_count / obs.raw_sum); raw_sum > 0",
                "excluded_collection_ids": sorted(set(args.exclude_collection_id)),
                "excluded_collection_names": sorted(set(args.exclude_collection_name)),
                "filters": "Homo sapiens; is_primary_data == True; cell type/tissue/disease/assay stratified"}
    write_json_atomic(out / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


def read_plan(out):
    out = Path(out)
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    features = [{**r, "soma_joinid": int(r["soma_joinid"]), "is_tf": r["is_tf"] == "True"} for r in read_tsv(out / "features.tsv")]
    datasets = list(read_tsv(out / "datasets.tsv"))
    return manifest, features, datasets


def get_present_feature_ids(census, dataset_joinid, features):
    matrix = census["census_data"]["homo_sapiens"].ms["RNA"]["feature_dataset_presence_matrix"]
    coords = (np.asarray([dataset_joinid], dtype=np.int64), np.asarray([f["soma_joinid"] for f in features], dtype=np.int64))
    present_joinids = set()
    for table in matrix.read(coords=coords).tables():
        data = table.to_pydict()
        for jid, present in zip(data["soma_dim_1"], data["soma_data"]):
            if present:
                present_joinids.add(int(jid))
    return present_joinids


def load_obs_to_sqlite(census, dataset_id, conn):
    conn.execute("CREATE TABLE cells (jid INTEGER PRIMARY KEY, cell_type_id TEXT, cell_type TEXT, tissue_id TEXT, tissue TEXT, disease_id TEXT, disease TEXT, assay_id TEXT, assay TEXT, donor_id TEXT, raw_sum REAL)")
    obs = census["census_data"]["homo_sapiens"].obs
    if not ID_RE.fullmatch(dataset_id):
        raise ValueError("Invalid dataset ID")
    filt = f"is_primary_data == True and dataset_id == '{dataset_id}'"
    n = 0
    insert = "INSERT INTO cells VALUES (?,?,?,?,?,?,?,?,?,?,?)"
    for table in obs.read(value_filter=filt, column_names=OBS_COLUMNS):
        for arrow_batch in table.to_batches(max_chunksize=10000):
            batch = []
            for r in arrow_batch.to_pylist():
                raw_sum = float(r["raw_sum"] or 0)
                if raw_sum <= 0 or not math.isfinite(raw_sum):
                    continue
                batch.append((int(r["soma_joinid"]), clean(r["cell_type_ontology_term_id"]), clean(r["cell_type"]),
                              clean(r["tissue_general_ontology_term_id"]), clean(r["tissue_general"]),
                              clean(r["disease_ontology_term_id"]), clean(r["disease"]),
                              clean(r["assay_ontology_term_id"]), clean(r["assay"]), clean(r["donor_id"]), raw_sum))
            conn.executemany(insert, batch)
            n += len(batch)
    conn.commit()
    conn.execute("CREATE INDEX group_idx ON cells (cell_type_id, tissue_id, disease_id, assay_id)")
    return n


def strata(conn, min_cells, assay_ids):
    sql = "SELECT cell_type_id,cell_type,tissue_id,tissue,disease_id,disease,assay_id,assay,COUNT(*) AS n FROM cells GROUP BY cell_type_id,tissue_id,disease_id,assay_id HAVING COUNT(*) >= ? ORDER BY n DESC"
    for row in conn.execute(sql, (min_cells,)):
        if not row[0] or row[0] == "unknown" or (assay_ids and row[6] not in assay_ids):
            continue
        yield dict(zip([*KEY_COLUMNS, "n"], row))


def cell_batches(conn, key, batch_size):
    sql = "SELECT jid,raw_sum,donor_id FROM cells WHERE cell_type_id=? AND tissue_id=? AND disease_id=? AND assay_id=? ORDER BY jid"
    cur = conn.execute(sql, (key["cell_type_id"], key["tissue_id"], key["disease_id"], key["assay_id"]))
    while batch := cur.fetchmany(batch_size):
        yield batch


def bh_qvalues(p):
    """Benjamini-Hochberg q-values for a 1-D array (NaN-safe: NaNs stay NaN, excluded from m)."""
    p = np.asarray(p, dtype=np.float64)
    q = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    m = len(ok)
    if not m:
        return q
    order = ok[np.argsort(p[ok], kind="stable")]
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1.0)
    return q


def corr_pvalues(r, df):
    """Two-sided p for Pearson/partial correlation r with df residual degrees of freedom (t test)."""
    r = np.clip(np.asarray(r, dtype=np.float64), -0.999999999, 0.999999999)
    if df is None or df < 1:
        return np.full(r.shape, np.nan)
    t = r * np.sqrt(df / (1.0 - r * r))
    return 2.0 * sstats.t.sf(np.abs(t), df)


def scatter_to_corr(w):
    d = np.sqrt(np.maximum(np.diag(w), 0))
    den = np.outer(d, d)
    return np.clip(np.divide(w, den, out=np.zeros_like(w), where=den > 0), -1, 1)


def grouped_partial_corr(y, lib, groups):
    """Partial correlation of columns of y (units x t) after removing group means (fixed effects)
    and a linear library-size covariate. Returns (corr matrix, residual df) or (None, None)."""
    y = np.asarray(y, dtype=np.float64)
    lib = np.asarray(lib, dtype=np.float64)
    groups = np.asarray(groups)
    n = y.shape[0]
    uniq, inv = np.unique(groups, return_inverse=True)
    g = len(uniq)
    yc = y.copy(); lc = lib.copy()
    for k in range(g):
        idx = inv == k
        yc[idx] -= yc[idx].mean(axis=0)
        lc[idx] -= lc[idx].mean()
    k_cov = 0
    if np.dot(lc, lc) > 1e-12:
        beta = (lc @ yc) / np.dot(lc, lc)
        yc -= np.outer(lc, beta)
        k_cov = 1
    df = n - g - k_cov - 1
    if df < 2:
        return None, None
    return scatter_to_corr(yc.T @ yc), df


class Accumulator:
    """Streaming sufficient statistics for one stratum.

    cell_raw   : cell-level Pearson r on log1p CPTT (original v1 statistic)
    resid_cell : cell-level partial r, within-donor (donor fixed effects) and log-library-size adjusted
    metacell   : partial r over within-donor metacells of `metacell_size` cells (donor FE + library size)
    pseudobulk : partial r over donor means (donors with >= min_donor_cells cells; library size adjusted)
    """

    def __init__(self, features, metacell_size=20, min_donor_cells=10, track_donor_covariance=False):
        self.features = features
        self.track_donor_covariance = track_donor_covariance
        self.donor_sum_sq = {}
        self.donor_gram = {}
        self.metacell_size = metacell_size
        self.min_donor_cells = min_donor_cells
        self.tf_indices = np.array([i for i, f in enumerate(features) if f["is_tf"]], dtype=np.int64)
        self.n = 0
        self.n_detected = np.zeros(len(features), dtype=np.int64)
        self.sum_log = np.zeros(len(features), dtype=np.float64)
        self.sum_raw = np.zeros(len(features), dtype=np.float64)
        t = len(self.tf_indices)
        self.sum_log_tf = np.zeros(t, dtype=np.float64)
        self.sum_sq_tf = np.zeros(t, dtype=np.float64)
        self.gram = np.zeros((t, t), dtype=np.float64)
        self.codetect = np.zeros((t, t), dtype=np.int64)
        self.sum_l = 0.0
        self.sum_ll = 0.0
        self.sum_yl = np.zeros(t, dtype=np.float64)
        # per group (donor; unknown donors pooled as one group): [n, sum_y, sum_l]
        self.groups = defaultdict(lambda: [0, np.zeros(t, dtype=np.float64), 0.0])
        # metacells: key (group, index) -> [n, sum_y, sum_l]
        self.mc = defaultdict(lambda: [0, np.zeros(t, dtype=np.float64), 0.0])

    @property
    def donors(self):
        return {d: v for d, v in self.groups.items() if d != "__unknown__"}

    @staticmethod
    def group_of(donor):
        return "__unknown__" if str(donor).lower() in ("", "unknown", "na", "nan", "none") else str(donor)

    def add(self, x, raw_sums, donor_ids):
        x = sparse.csr_matrix(x, dtype=np.float64)
        if x.shape[1] != len(self.features):
            raise ValueError("Expression columns do not match requested features")
        raw_sums = np.asarray(raw_sums, dtype=np.float64)
        if len(raw_sums) != x.shape[0] or np.any(raw_sums <= 0):
            raise ValueError("Invalid raw_sum alignment")
        y = x.multiply((10000.0 / raw_sums)[:, None]).tocsr()
        y.data = np.log1p(y.data)
        lib = np.log(raw_sums)
        self.n += x.shape[0]
        self.n_detected += np.asarray((x > 0).sum(axis=0)).ravel().astype(np.int64)
        self.sum_log += np.asarray(y.sum(axis=0)).ravel()
        self.sum_raw += np.asarray(x.sum(axis=0)).ravel()
        yt = y[:, self.tf_indices].tocsr()
        self.sum_log_tf += np.asarray(yt.sum(axis=0)).ravel()
        self.sum_sq_tf += np.asarray(yt.power(2).sum(axis=0)).ravel()
        self.gram += (yt.T @ yt).toarray()
        binary = (x[:, self.tf_indices] > 0).astype(np.int64)
        self.codetect += (binary.T @ binary).toarray()
        self.sum_l += float(lib.sum())
        self.sum_ll += float((lib * lib).sum())
        self.sum_yl += np.asarray(yt.T @ lib).ravel()
        # assign each cell to a group and a within-group metacell (deterministic, in jid order)
        gkeys = [self.group_of(d) for d in donor_ids]
        mkeys = []
        for gk in gkeys:
            cnt = self.groups[gk][0]
            self.groups[gk][0] = cnt + 1
            mkeys.append((gk, cnt // self.metacell_size))
        if self.track_donor_covariance:   # per-donor second moments (Tabula Sapiens validation)
            for donor in set(gkeys):
                if donor == "__unknown__":
                    continue
                rows = np.array([i for i, g in enumerate(gkeys) if g == donor], dtype=int)
                dy = yt[rows]
                t = len(self.tf_indices)
                if donor not in self.donor_sum_sq:
                    self.donor_sum_sq[donor] = np.zeros(t); self.donor_gram[donor] = np.zeros((t, t))
                self.donor_sum_sq[donor] += np.asarray(dy.power(2).sum(axis=0)).ravel()
                self.donor_gram[donor] += (dy.T @ dy).toarray()
        for keys, store in ((gkeys, self.groups), (mkeys, self.mc)):
            uniq = {}
            rows = np.array([uniq.setdefault(k, len(uniq)) for k in keys], dtype=np.int64)
            ind = sparse.csr_matrix((np.ones(len(rows)), (rows, np.arange(len(rows)))), shape=(len(uniq), len(rows)))
            sy = (ind @ yt).toarray()
            sl = ind @ lib
            cnt = np.asarray(ind.sum(axis=1)).ravel()
            for k, i in uniq.items():
                cur = store[k]
                if store is self.mc:
                    cur[0] += int(cnt[i])
                cur[1] += sy[i]
                cur[2] += float(sl[i])

    def profiles(self, metadata):
        n_donors = len([d for d, v in self.donors.items() if v[0] > 0])
        for i, gene in enumerate(self.features):
            yield {**metadata, "ensembl_id": gene["ensembl_id"], "symbol": gene["symbol"],
                   "hgnc_id": gene["hgnc_id"], "category": gene["category"], "n_cells": self.n,
                   "n_donors": n_donors, "n_detected": int(self.n_detected[i]),
                   "detection_fraction": self.n_detected[i] / self.n,
                   "mean_log1p_cptt": self.sum_log[i] / self.n,
                   "mean_raw_count": self.sum_raw[i] / self.n}

    # ---- statistics -------------------------------------------------------
    def cell_raw_corr(self):
        n, sums = self.n, self.sum_log_tf
        ss = self.sum_sq_tf - sums * sums / n
        denom = np.sqrt(np.maximum(ss[:, None], 0) * np.maximum(ss[None, :], 0))
        cov = self.gram - np.outer(sums, sums) / n
        return np.clip(np.divide(cov, denom, out=np.zeros_like(cov), where=denom > 0), -1, 1)

    def resid_cell_corr(self):
        """Within-group scatter with library-size partialled out, from sufficient statistics."""
        groups = [v for v in self.groups.values() if v[0] > 0]
        w = self.gram.copy()
        wyl = self.sum_yl.copy()
        wll = self.sum_ll
        for cnt, sy, sl in groups:
            w -= np.outer(sy, sy) / cnt
            wyl -= sy * sl / cnt
            wll -= sl * sl / cnt
        k_cov = 0
        if wll > 1e-12:
            w -= np.outer(wyl, wyl) / wll
            k_cov = 1
        df = self.n - len(groups) - k_cov - 1
        if df < 2:
            return None, None
        return scatter_to_corr(w), df

    def metacell_corr(self):
        keep = [(k, v) for k, v in self.mc.items() if v[0] >= max(2, self.metacell_size // 2)]
        if len(keep) < 5:
            return None, None, len(keep)
        y = np.vstack([v[1] / v[0] for _, v in keep])
        lib = np.array([v[2] / v[0] for _, v in keep])
        grp = np.array([k[0] for k, _ in keep])
        r, df = grouped_partial_corr(y, lib, grp)
        return r, df, len(keep)

    def pseudobulk_corr(self):
        keep = [v for d, v in self.donors.items() if v[0] >= self.min_donor_cells]
        if len(keep) < 4:
            return None, None, len(keep)
        y = np.vstack([v[1] / v[0] for v in keep])
        lib = np.array([v[2] / v[0] for v in keep])
        r, df = grouped_partial_corr(y, lib, np.zeros(len(keep)))
        return r, df, len(keep)

    def pairs(self, metadata, min_detection, min_rho, min_lift, min_codetected, top_k,
              primary="metacell", max_q=0.05, store_p=0.0, store_min_r=0.1):
        """Yield pair rows. Rows with passed=True are screened hits; with store_p > 0, additional
        tested pairs with primary p <= store_p and primary r >= store_min_r are kept (passed=False)
        so that the aggregate meta-analysis is less biased by per-dataset selection."""
        if self.n < 3 or len(self.tf_indices) < 2:
            return
        n = self.n
        raw = self.cell_raw_corr()
        detected = self.n_detected[self.tf_indices]
        frac = detected / n
        lift_denom = np.outer(detected, detected)
        lift = np.divide(self.codetect * n, lift_denom, out=np.zeros_like(raw), where=lift_denom > 0)
        det_ok = (frac[:, None] >= min_detection) & (frac[None, :] >= min_detection)
        det_ok = np.triu(det_ok, 1)
        res_r, res_df = self.resid_cell_corr()
        mc_r, mc_df, n_mc = self.metacell_corr()
        pb_r, pb_df, n_pb = self.pseudobulk_corr()
        avail = {"cell_raw": (raw, n - 2), "resid_cell": (res_r, res_df), "metacell": (mc_r, mc_df), "pseudobulk": (pb_r, pb_df)}
        chain = {"cell_raw": ["cell_raw"], "resid_cell": ["resid_cell", "cell_raw"],
                 "metacell": ["metacell", "resid_cell", "cell_raw"], "pseudobulk": ["pseudobulk", "metacell", "resid_cell", "cell_raw"]}[primary]
        used = next(s for s in chain if avail[s][0] is not None)
        pr, pdf = avail[used]
        ii, jj = np.nonzero(det_ok)
        pvals = corr_pvalues(pr[ii, jj], pdf)
        qvals = bh_qvalues(pvals)
        n_tested = len(ii)
        if primary == "cell_raw":
            passed = (pr[ii, jj] >= min_rho) & (lift[ii, jj] >= min_lift) & (self.codetect[ii, jj] >= min_codetected)
        else:
            passed = ((pr[ii, jj] >= min_rho) & (qvals <= max_q) & (lift[ii, jj] >= min_lift)
                      & (self.codetect[ii, jj] >= min_codetected))
        if top_k:
            keep = np.zeros(len(ii), dtype=bool)
            for t in range(len(self.tf_indices)):
                sel = np.flatnonzero(passed & ((ii == t) | (jj == t)))
                keep[sel[np.argsort(-pr[ii[sel], jj[sel]], kind="stable")[:top_k]]] = True
            passed &= keep
        stored = passed.copy()
        if store_p > 0:
            stored |= (pvals <= store_p) & (pr[ii, jj] >= store_min_r)
        g = lambda m, i, j: "" if m is None else float(m[i, j])
        for k in np.flatnonzero(stored):
            i, j = int(ii[k]), int(jj[k])
            a = self.features[self.tf_indices[i]]
            b = self.features[self.tf_indices[j]]
            yield {**metadata, "tf1_ensembl": a["ensembl_id"], "tf1_symbol": a["symbol"],
                   "tf2_ensembl": b["ensembl_id"], "tf2_symbol": b["symbol"],
                   "n_cells": n, "n_donors": n_pb, "tf1_detection_fraction": frac[i],
                   "tf2_detection_fraction": frac[j], "co_detected_cells": int(self.codetect[i, j]),
                   "co_detection_lift": lift[i, j], "cell_pearson_r": float(raw[i, j]),
                   "donor_pearson_r": g(pb_r, i, j),
                   "resid_cell_r": g(res_r, i, j), "resid_cell_df": res_df or "",
                   "metacell_r": g(mc_r, i, j), "metacell_df": mc_df or "", "n_metacells": n_mc,
                   "pseudobulk_r": g(pb_r, i, j), "pseudobulk_df": pb_df or "",
                   "primary_stat": used, "primary_r": float(pr[i, j]), "primary_df": pdf,
                   "p_value": float(pvals[k]), "q_value_bh": float(qvals[k]),
                   "n_tested_pairs": n_tested, "passed": bool(passed[k])}


def process_dataset(census, dataset, features, args):
    dataset_id = dataset["dataset_id"]
    target = dataset_path(args.out, dataset_id)
    if (target / "_SUCCESS").exists():
        prior = json.loads((target / "_SUCCESS").read_text(encoding="utf-8"))
        if prior["run_config"] != args.run_config:
            raise RuntimeError(f"Completed dataset {dataset_id} used different run parameters; use a new plan/output directory")
        print(f"SKIP completed {dataset_id}", flush=True)
        return
    target.mkdir(parents=True, exist_ok=True)
    present_joinids = get_present_feature_ids(census, int(dataset["soma_joinid"]), features)
    present = [f for f in features if f["soma_joinid"] in present_joinids]
    metadata_base = {"dataset_id": dataset_id, "collection_id": dataset["collection_id"]}
    write_tsv(target / "feature_availability.tsv.gz", ({"dataset_id": dataset_id, "ensembl_id": f["ensembl_id"], "symbol": f["symbol"], "census_feature_present": f["soma_joinid"] in present_joinids} for f in features), AVAIL_COLUMNS)
    if not present:
        write_tsv(target / "profiles.tsv.gz", [], PROFILE_COLUMNS)
        write_tsv(target / "pairs.tsv.gz", [], PAIR_COLUMNS)
        write_json_atomic(target / "_SUCCESS", {"run_config": args.run_config, "n_primary_cells": 0, "n_strata": 0, "n_profile_rows": 0, "n_pair_rows": 0})
        return
    with tempfile.TemporaryDirectory(prefix="tf_census_", dir=args.tmpdir) as tempdir:
        conn = sqlite3.connect(str(Path(tempdir) / "obs.sqlite"))
        n_meta = load_obs_to_sqlite(census, dataset_id, conn)
        if not n_meta:
            write_tsv(target / "profiles.tsv.gz", [], PROFILE_COLUMNS)
            write_tsv(target / "pairs.tsv.gz", [], PAIR_COLUMNS)
            write_json_atomic(target / "_SUCCESS", {"run_config": args.run_config, "n_primary_cells": 0, "n_strata": 0, "n_profile_rows": 0, "n_pair_rows": 0})
            return
        feature_ids = [f["soma_joinid"] for f in present]
        group_count = profile_count = pair_count = pass_count = 0
        with tsv_writer(target / "profiles.tsv.gz", PROFILE_COLUMNS) as profile_writer, tsv_writer(target / "pairs.tsv.gz", PAIR_COLUMNS) as pair_writer:
            for key in strata(conn, args.min_cells, set(args.assay_id)):
                if args.max_strata and group_count >= args.max_strata:
                    break
                acc = Accumulator(present, metacell_size=args.metacell_size, min_donor_cells=args.min_donor_cells)
                for batch in cell_batches(conn, key, args.batch_size):
                    ids = np.asarray([r[0] for r in batch], dtype=np.int64)
                    cells = {r[0]: r for r in batch}
                    adata = census_module().get_anndata(
                        census, organism="Homo sapiens", measurement_name="RNA", X_name="raw",
                        obs_coords=ids, var_coords=np.asarray(feature_ids, dtype=np.int64),
                        obs_column_names=["soma_joinid"], var_column_names=["soma_joinid", "feature_id"])
                    obs_ids = [int(v) for v in adata.obs["soma_joinid"]]
                    var_ids = [int(v) for v in adata.var["soma_joinid"]]
                    if len(obs_ids) != len(batch) or set(obs_ids) != set(ids) or len(var_ids) != len(feature_ids) or set(var_ids) != set(feature_ids):
                        raise RuntimeError(f"Census coordinate mismatch for {dataset_id}; stopping instead of misaligning counts")
                    var_pos = {jid: i for i, jid in enumerate(var_ids)}
                    x = adata.X[:, [var_pos[jid] for jid in feature_ids]]
                    acc.add(x, [cells[j][1] for j in obs_ids], [cells[j][2] for j in obs_ids])
                if acc.n != key["n"]:
                    raise RuntimeError(f"Incomplete stratum {dataset_id}: {acc.n} versus {key['n']}")
                group_metadata = {**metadata_base, **{k: key[k] for k in KEY_COLUMNS}}
                for row in acc.profiles(group_metadata):
                    profile_writer.writerow(row)
                    profile_count += 1
                for row in acc.pairs(group_metadata, args.min_detection, args.min_rho, args.min_lift, args.min_codetected, args.top_k,
                                          primary=args.primary_stat, max_q=args.max_q, store_p=args.store_p, store_min_r=args.store_min_r):
                    pair_writer.writerow(row)
                    pair_count += 1
                    pass_count += bool(row.get("passed", True))
                group_count += 1
                print(f"{dataset_id}: {group_count} strata; {acc.n} cells; {pass_count} passing pairs ({pair_count} stored rows)", flush=True)
        write_json_atomic(target / "_SUCCESS", {"run_config": args.run_config, "n_primary_cells": n_meta, "n_strata": group_count, "n_profile_rows": profile_count, "n_pair_rows": pair_count, "n_passing_pairs": pass_count})
        conn.close()


def command_run(args):
    manifest, features, datasets = read_plan(args.out)
    if args.num_shards < 1 or args.batch_size < 1 or args.min_cells < 1 or args.top_k < 0:
        raise ValueError("Shard count, batch size, and minimum cells must be positive; top-k must be nonnegative")
    if not 0 <= args.min_detection <= 1 or not -1 <= args.min_rho <= 1 or args.min_lift < 0 or args.min_codetected < 0 or args.max_strata < 0:
        raise ValueError("Invalid screening threshold")
    if not 0 <= args.shard_index < args.num_shards:
        raise ValueError("--shard-index must be in [0, --num-shards)")
    selected = [d for d in datasets if int(hashlib.sha256(d["dataset_id"].encode()).hexdigest(), 16) % args.num_shards == args.shard_index]
    args.run_config = {"min_cells": args.min_cells, "min_detection": args.min_detection,
                       "min_rho": args.min_rho, "min_lift": args.min_lift,
                       "min_codetected": args.min_codetected, "top_k": args.top_k,
                       "assay_ids": sorted(set(args.assay_id)), "max_strata": args.max_strata,
                       "schema_version": SCHEMA_VERSION, "primary_stat": args.primary_stat, "max_q": args.max_q,
                       "metacell_size": args.metacell_size, "min_donor_cells": args.min_donor_cells,
                       "store_p": args.store_p, "store_min_r": args.store_min_r}
    print(f"Pinned Census {manifest['resolved_census_version']}; shard {args.shard_index}/{args.num_shards}; {len(selected)} datasets", flush=True)
    api = census_module()
    with api.open_soma(census_version=manifest["resolved_census_version"]) as census:
        for d in selected:
            process_dataset(census, d, features, args)


def meta_analyse(rows):
    """Fisher-z meta-analysis over (r, df) per dataset. Weight = df - 1 (~ n_eff - 3).
    Returns [k, r_fixed, p_fixed, r_random, p_random, q placeholder, tau2, I2]."""
    zs, ws = [], []
    for r, df in rows:
        if df is None or df < 2:
            continue
        zs.append(math.atanh(max(-0.999999, min(0.999999, r))))
        ws.append(df - 1.0)
    k = len(zs)
    if not k:
        return [0, "", "", "", "", "", "", ""]
    z, w = np.array(zs), np.array(ws)
    zf = float(np.sum(w * z) / np.sum(w)); sef = 1 / math.sqrt(np.sum(w))
    pf = float(2 * sstats.norm.sf(abs(zf / sef)))
    q = float(np.sum(w * (z - zf) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    wr = 1 / (1 / w + tau2)
    zr = float(np.sum(wr * z) / np.sum(wr)); ser = 1 / math.sqrt(np.sum(wr))
    pr = float(2 * sstats.norm.sf(abs(zr / ser)))
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    return [k, math.tanh(zf), pf, math.tanh(zr), pr, "", tau2, i2]


def command_aggregate(args):
    root = Path(args.out)
    completed = sorted(p.parent for p in (root / "datasets").glob("*/_SUCCESS"))
    if not completed:
        raise SystemExit("No completed dataset outputs found")
    config_set = {json.dumps(json.loads((folder / "_SUCCESS").read_text(encoding="utf-8"))["run_config"], sort_keys=True) for folder in completed}
    if len(config_set) != 1:
        raise RuntimeError("Completed datasets used different run parameters; refusing to mix them")
    run_config = json.loads(next(iter(config_set)))
    planned = len(list(read_tsv(root / "datasets.tsv")))
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    ensembl_hgnc = {r["ensembl_id"]: r["hgnc_id"] for r in read_tsv(root / "features.tsv")}
    prior_lookup = defaultdict(lambda: {"classes": set(), "evidence": set(), "pmids": set()})
    if getattr(args, "interaction_priors", None):
        for r in read_tsv(args.interaction_priors):
            if not r["A_HGNC_ID"] or not r["B_HGNC_ID"]:
                continue
            key = tuple(sorted((r["A_HGNC_ID"], r["B_HGNC_ID"])))
            prior_lookup[key]["classes"].add(r["Source_pair_class"])
            prior_lookup[key]["evidence"].add(r["Best_source_evidence"])
            prior_lookup[key]["pmids"].update(v.strip() for v in re.split(r"[,;|]", r["PMIDs"]) if v.strip())
    target = root / "aggregate"
    target.mkdir(exist_ok=True)
    db = target / "aggregate.sqlite"
    if db.exists():
        raise SystemExit("Aggregate database already exists; use a fresh output directory or move the old aggregate aside")
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("CREATE TABLE profiles (dataset_id TEXT, collection_id TEXT, cell_type_id TEXT, cell_type TEXT, tissue_id TEXT, tissue TEXT, disease_id TEXT, disease TEXT, assay_id TEXT, assay TEXT, ensembl_id TEXT, symbol TEXT, hgnc_id TEXT, category TEXT, n_cells INTEGER, n_donors INTEGER, n_detected INTEGER, sum_log REAL, sum_raw REAL)")
    pair_cols = ["dataset_id", "collection_id", *KEY_COLUMNS, "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
                 "n_cells", "n_donors", "co_detected_cells", "co_detection_lift", "cell_pearson_r", "donor_pearson_r",
                 "resid_cell_r", "metacell_r", "pseudobulk_r", "primary_stat", "primary_r", "primary_df", "p_value", "q_value_bh", "passed"]
    conn.execute(f"CREATE TABLE pairs ({', '.join(pair_cols)})")
    for folder in completed:
        with open_tsv(folder / "profiles.tsv.gz") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                conn.execute("INSERT INTO profiles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                             tuple(r[k] for k in ["dataset_id", "collection_id", *KEY_COLUMNS, "ensembl_id", "symbol", "hgnc_id", "category", "n_cells", "n_donors", "n_detected"])
                             + (int(r["n_cells"]) * float(r["mean_log1p_cptt"]), int(r["n_cells"]) * float(r["mean_raw_count"])))
        with open_tsv(folder / "pairs.tsv.gz") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                vals = []
                for k in pair_cols:
                    v = r.get(k, "")
                    if k == "passed":
                        v = 1 if v in ("True", "") else 0   # v1 outputs only contained hits
                    elif k == "primary_stat" and not v:
                        v = "cell_raw"
                    elif k == "primary_r" and v == "":
                        v = r["cell_pearson_r"]
                    vals.append(None if v == "" else v)
                conn.execute(f"INSERT INTO pairs VALUES ({','.join('?' * len(pair_cols))})", vals)
        conn.commit()
    profile_sql = """SELECT cell_type_id,cell_type,tissue_id,tissue,disease_id,disease,assay_id,assay,ensembl_id,symbol,hgnc_id,category,
      COUNT(DISTINCT dataset_id),COUNT(DISTINCT collection_id),SUM(n_cells),SUM(n_detected),SUM(n_detected)*1.0/SUM(n_cells),SUM(sum_log)/SUM(n_cells),SUM(sum_raw)/SUM(n_cells)
      FROM profiles GROUP BY cell_type_id,tissue_id,disease_id,assay_id,ensembl_id ORDER BY cell_type_id,tissue_id,disease_id,assay_id,ensembl_id"""
    with open(target / "cell_type_profiles.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([*KEY_COLUMNS, "ensembl_id", "symbol", "hgnc_id", "category", "n_datasets", "n_collections", "n_cells", "n_detected", "detection_fraction", "mean_log1p_cptt", "mean_raw_count"])
        writer.writerows(conn.execute(profile_sql))
    cand_header = [*KEY_COLUMNS, "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol", "n_datasets_with_hit",
                   "n_collections_with_hit", "n_cells_in_hit_datasets", "co_detected_cells", "cell_r_weighted_by_cells",
                   "mean_co_detection_lift", "mean_donor_r_when_available", "primary_stat", "best_primary_r", "min_q_within_stratum",
                   "n_datasets_in_meta", "meta_r_fixed", "meta_p_fixed", "meta_r_random", "meta_p_random", "meta_q_random",
                   "meta_tau2", "meta_i2", "physical_prior_recorded", "prior_best_evidence", "prior_source_classes", "prior_pmids"]
    key_sql = "cell_type_id,tissue_id,disease_id,assay_id,tf1_ensembl,tf2_ensembl"
    cur = conn.execute(f"SELECT {', '.join(pair_cols)} FROM pairs ORDER BY {key_sql}")
    ix = {c: i for i, c in enumerate(pair_cols)}
    out_rows = []
    import itertools
    keyf = lambda r: tuple(r[ix[c]] for c in key_sql.split(","))
    for key, grp in itertools.groupby(cur, key=keyf):
        grp = list(grp)
        hits = [r for r in grp if int(r[ix["passed"]]) == 1]
        meta = meta_analyse([(float(r[ix["primary_r"]]), float(r[ix["primary_df"]]) if r[ix["primary_df"]] not in (None, "") else None) for r in grp])
        if not hits:
            continue
        f = grp[0]
        ncell = sum(int(r[ix["n_cells"]]) for r in hits)
        donors = [float(r[ix["donor_pearson_r"]]) for r in hits if r[ix["donor_pearson_r"]] not in (None, "")]
        qs = [float(r[ix["q_value_bh"]]) for r in hits if r[ix["q_value_bh"]] not in (None, "")]
        row = [*(f[ix[c]] for c in KEY_COLUMNS), f[ix["tf1_ensembl"]], f[ix["tf1_symbol"]], f[ix["tf2_ensembl"]], f[ix["tf2_symbol"]],
               len({r[ix["dataset_id"]] for r in hits}), len({r[ix["collection_id"]] for r in hits}), ncell,
               sum(int(r[ix["co_detected_cells"]]) for r in hits),
               sum(float(r[ix["cell_pearson_r"]]) * int(r[ix["n_cells"]]) for r in hits) / ncell,
               sum(float(r[ix["co_detection_lift"]]) for r in hits) / len(hits),
               (sum(donors) / len(donors)) if donors else "",
               "|".join(sorted({r[ix["primary_stat"]] for r in hits})),
               max(float(r[ix["primary_r"]]) for r in hits), min(qs) if qs else "", *meta]
        out_rows.append(row)
    # BH over random-effects meta p within stratum
    by_stratum = defaultdict(list)
    for i, row in enumerate(out_rows):
        by_stratum[tuple(row[0:8:2])].append(i)
    p_col = cand_header.index("meta_p_random")
    for idxs in by_stratum.values():
        ps = [out_rows[i][p_col] if out_rows[i][p_col] != "" else np.nan for i in idxs]
        for i, q in zip(idxs, bh_qvalues(ps)):
            out_rows[i][p_col + 1] = "" if not np.isfinite(q) else float(q)
    nc, rc = cand_header.index("n_collections_with_hit"), cand_header.index("cell_r_weighted_by_cells")
    out_rows.sort(key=lambda r: (-r[nc], -r[rc]))
    with open(target / "cell_type_pair_candidates.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(cand_header)
        for row in out_rows:
            hgnc_a = ensembl_hgnc.get(row[8], "")
            hgnc_b = ensembl_hgnc.get(row[10], "")
            prior = prior_lookup.get(tuple(sorted((hgnc_a, hgnc_b)))) if hgnc_a and hgnc_b else None
            if prior:
                best = "Direct-interaction annotation" if "Direct-interaction annotation" in prior["evidence"] else "Association only"
                annotation = ["Yes", best, " | ".join(sorted(prior["classes"])), " | ".join(sorted(prior["pmids"]))]
            else:
                annotation = ["No", "", "", ""]
            writer.writerow([*row, *annotation])
    # positive-control check
    sym_rows = defaultdict(list)
    for row in out_rows:
        sym_rows[tuple(sorted((row[9], row[11])))].append(row)
    controls = []
    for a_, b_ in POSITIVE_CONTROLS:
        rs = sym_rows.get(tuple(sorted((a_, b_))), [])
        controls.append({"pair": f"{a_}-{b_}", "n_contexts_with_hit": len(rs),
                         "max_collections": max((r[nc] for r in rs), default=0),
                         "best_primary_r": max((r[cand_header.index("best_primary_r")] for r in rs), default=None),
                         "min_q": min((r[cand_header.index("min_q_within_stratum")] for r in rs if r[cand_header.index("min_q_within_stratum")] != ""), default=None)})
    write_tsv(target / "positive_controls.tsv", controls, ["pair", "n_contexts_with_hit", "max_collections", "best_primary_r", "min_q"])
    conn.close()
    plan_complete = len(completed) == planned
    broad_scope_complete = plan_complete and manifest.get("dataset_selection") == "all" and run_config.get("max_strata", 0) == 0
    held_out = bool(manifest.get("excluded_collection_ids") or manifest.get("excluded_collection_names"))
    write_json_atomic(target / "summary.json", {"n_planned_datasets": planned, "n_completed_datasets": len(completed),
                                                "complete": plan_complete,
                                                "atlas_wide": broad_scope_complete and not held_out,
                                                "atlas_wide_except_heldout_collections": broad_scope_complete and held_out,
                                                "excluded_collection_ids": manifest.get("excluded_collection_ids", []),
                                                "excluded_collection_names": manifest.get("excluded_collection_names", []),
                                                "dataset_selection": manifest.get("dataset_selection", "unknown"),
                                                "run_config": run_config, "schema_version": SCHEMA_VERSION,
                                                "n_candidate_rows": len(out_rows),
                                                "positive_controls_recovered": sum(1 for c in controls if c["n_contexts_with_hit"]),
                                                "positive_controls_tested": len(controls)})
    print(f"Aggregated {len(completed)} datasets into {target}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    plan = sub.add_parser("plan", help="Resolve a Census release, genes, and datasets")
    plan.add_argument("--register", required=True, help="gene_register.tsv from TF compendium v0.2")
    plan.add_argument("--out", required=True)
    plan.add_argument("--census-version", default="stable", help="stable or an explicit Census release date")
    plan.add_argument("--tf-set", choices=["strict", "exploratory", "expanded"], default="exploratory")
    plan.add_argument("--dataset-id", action="append", default=[], help="Repeat to select pilot datasets")
    plan.add_argument("--exclude-collection-id", action="append", default=[], help="Exclude a source collection from discovery; repeatable")
    plan.add_argument("--exclude-collection-name", action="append", default=[], help="Exclude a source collection by its exact release-specific name")
    plan.add_argument("--max-datasets", type=int, default=0, help="Pilot only; first N datasets by ID")
    plan.set_defaults(func=command_plan)
    run = sub.add_parser("run", help="Mine dataset shards from a pinned Census")
    run.add_argument("--out", required=True)
    run.add_argument("--shard-index", type=int, default=0)
    run.add_argument("--num-shards", type=int, default=1)
    run.add_argument("--batch-size", type=int, default=1000)
    run.add_argument("--tmpdir", default=None, help="Node-local scratch directory for temporary cell metadata")
    run.add_argument("--min-cells", type=int, default=100)
    run.add_argument("--min-detection", type=float, default=0.05)
    run.add_argument("--min-rho", type=float, default=0.25)
    run.add_argument("--min-lift", type=float, default=1.1)
    run.add_argument("--min-codetected", type=int, default=25)
    run.add_argument("--primary-stat", choices=["metacell", "pseudobulk", "resid_cell", "cell_raw"], default="metacell",
                     help="Statistic driving p/q and selection; falls back down the chain if unavailable. cell_raw = v1 behaviour (no FDR filter)")
    run.add_argument("--max-q", type=float, default=0.05, help="BH q-value cutoff within stratum (ignored for cell_raw)")
    run.add_argument("--metacell-size", type=int, default=20, help="Cells per within-donor metacell")
    run.add_argument("--min-donor-cells", type=int, default=10, help="Min cells for a donor to enter pseudobulk")
    run.add_argument("--store-p", type=float, default=0.01, help="Also store non-passing tested pairs with primary p <= this (for meta-analysis); 0 disables")
    run.add_argument("--store-min-r", type=float, default=0.1, help="...and primary r >= this")
    run.add_argument("--top-k", type=int, default=0, help="Optional per-TF hit cap; 0 retains all passing pairs")
    run.add_argument("--assay-id", action="append", default=[], help="Optional assay ontology ID filter; repeatable")
    run.add_argument("--max-strata", type=int, default=0, help="Pilot only; cap strata per dataset")
    run.set_defaults(func=command_run)
    agg = sub.add_parser("aggregate", help="Combine all completed dataset outputs")
    agg.add_argument("--out", required=True)
    agg.add_argument("--interaction-priors", help="Optional interaction_priors.tsv from TF compendium v0.2")
    agg.set_defaults(func=command_aggregate)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
