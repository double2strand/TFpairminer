#!/usr/bin/env python3
"""Targeted donor-aware RNA co-expression test of ChromBPNet candidate pairs.

Box-scale complement to run_remote.sh: instead of an atlas-wide uncapped screen, it
tests only the TF pairs listed in pairs_cbp.tsv, only in a declared set of RNA cell
types (the context matched to the ATAC cell line), with the same donor-aware statistic
and thresholds as tf_census.py (within-donor metacell partial r; BH q<=0.05, r>=0.25,
co-detection lift>=1.1, >=25 co-detected cells, detection>=5%, >=100 cells/stratum).

BH is applied per stratum over the candidate pairs that pass detection there (a much
smaller family than an atlas screen, so q is less conservative than tf_census; reported).

Screens:  cellxgene = all Census datasets except the Tabula Sapiens collection
          tabula    = Tabula Sapiens "All Cells" only (non-overlapping replication)
Checkpoints: one JSON per stratum under OUT/strata/; reruns skip finished strata.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, sys, time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tf_census as tc  # noqa: E402

CONTEXTS = {
    # label, mapping level, CL ids. Cell lines are cancer lines: matches are PROXY, never exact.
    "hepatocyte": ("proxy: HepG2 (hepatoblastoma line) -> normal hepatocytes", [
        "CL:0000182", "CL:0019026", "CL:0019028", "CL:0019029"]),
    "erythroid_myeloid_progenitor": ("proxy: K562 (CML erythroleukemia line) -> normal erythroid/MEP/CMP progenitors", [
        "CL:0000038", "CL:0000049", "CL:0000050", "CL:0000547", "CL:0000549", "CL:0000550",
        "CL:0000552", "CL:0000764", "CL:0000765", "CL:0001066"]),
}


def candidate_pairs(paths, sym2id):
    """Return {(ensA, ensB)} from k=2 rows and all sub-pairs of k=3 rows; plus label rejects."""
    pairs, rejects = set(), []
    for p in paths:
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                names = row["tfs"].split("+")
                if any("::" in n for n in names):
                    rejects.append((row["tfs"], "composite_motif_label")); continue
                ids = [sym2id.get(n.strip().upper()) for n in names]
                if any(i is None for i in ids):
                    rejects.append((row["tfs"], "symbol_not_in_register")); continue
                for i in range(len(ids)):
                    for j in range(i + 1, len(ids)):
                        if ids[i] != ids[j]:
                            pairs.add(tuple(sorted((ids[i], ids[j]))))
    return pairs, rejects


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs-cbp", nargs="+", required=True)
    ap.add_argument("--register", required=True)
    ap.add_argument("--context", choices=sorted(CONTEXTS), required=True)
    ap.add_argument("--screen", choices=["cellxgene", "tabula"], required=True)
    ap.add_argument("--census-version", default="2025-11-08")
    ap.add_argument("--tabula-dataset-id", default="53d208b0-2cfd-4366-9866-c3c6114081bc")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cells", type=int, default=100)
    ap.add_argument("--min-detection", type=float, default=0.05)
    ap.add_argument("--min-rho", type=float, default=0.25)
    ap.add_argument("--min-lift", type=float, default=1.1)
    ap.add_argument("--min-codetected", type=int, default=25)
    ap.add_argument("--max-q", type=float, default=0.05)
    ap.add_argument("--metacell-size", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=20000)
    ap.add_argument("--normal-only", action="store_true", default=True)
    a = ap.parse_args()
    out = Path(a.out); (out / "strata").mkdir(parents=True, exist_ok=True)

    reg = {}
    sym2id = {}
    for row in tc.read_tsv(a.register):
        gid = row.get("Ensembl_gene_ID", "").split(".")[0]
        if gid.startswith("ENSG"):
            reg[gid] = row["Symbol"]; sym2id.setdefault(row["Symbol"].upper(), gid)
    pairs, rejects = candidate_pairs(a.pairs_cbp, sym2id)
    genes = sorted({g for p in pairs for g in p})
    label, cl_ids = CONTEXTS[a.context]
    cc = tc.census_module()
    with cc.open_soma(census_version=a.census_version) as census:
        hs = census["census_data"]["homo_sapiens"]
        ds = census["census_info"]["datasets"].read().concat().to_pandas()
        ts_ids = set(ds[ds.collection_name.str.contains("Tabula Sapiens", case=False)].dataset_id)
        if a.screen == "tabula":
            if a.tabula_dataset_id not in set(ds.dataset_id):
                raise SystemExit("Tabula dataset id not in this Census release")
        var = hs.ms["RNA"].var.read(column_names=["soma_joinid", "feature_id"]).concat().to_pandas()
        var["gid"] = var.feature_id.str.split(".").str[0]
        var = var[var.gid.isin(genes)]
        feats = [{"ensembl_id": r.gid, "symbol": reg[r.gid], "hgnc_id": "", "category": "", "is_tf": True,
                  "soma_joinid": int(r.soma_joinid)} for r in var.itertuples()]
        fidx = {f["ensembl_id"]: i for i, f in enumerate(feats)}
        cts = ", ".join(f"'{c}'" for c in cl_ids)
        filt = f"is_primary_data == True and cell_type_ontology_term_id in [{cts}]"
        if a.normal_only:
            filt += " and disease_ontology_term_id == 'PATO:0000461'"
        if a.screen == "tabula":
            filt += f" and dataset_id == '{a.tabula_dataset_id}'"
        cols = ["soma_joinid", "dataset_id", "cell_type_ontology_term_id", "cell_type", "tissue_general_ontology_term_id",
                "tissue_general", "assay_ontology_term_id", "assay", "donor_id", "raw_sum"]
        obs = hs.obs.read(value_filter=filt, column_names=cols).concat().to_pandas()
        if a.screen == "cellxgene":
            obs = obs[~obs.dataset_id.isin(ts_ids)]
        obs = obs[obs.raw_sum > 0]
        meta = ds.set_index("dataset_id")
        manifest = {"screen": a.screen, "context": a.context, "mapping": label, "cl_ids": cl_ids,
                    "census_version": a.census_version, "filter": filt,
                    "excluded_tabula_datasets": len(ts_ids) if a.screen == "cellxgene" else 0,
                    "n_cells": int(len(obs)), "n_datasets": int(obs.dataset_id.nunique()),
                    "n_candidate_pairs": len(pairs), "n_genes_in_census": len(feats), "rejected_labels": rejects,
                    "datasets": sorted(obs.dataset_id.unique()),
                    "thresholds": {k: getattr(a, k) for k in ("min_cells", "min_detection", "min_rho", "min_lift",
                                                               "min_codetected", "max_q", "metacell_size")}}
        tc.write_json_atomic(out / "manifest.json", manifest)
        print(f"{a.screen}/{a.context}: {len(obs)} cells, {obs.dataset_id.nunique()} datasets, {len(pairs)} pairs", flush=True)
        keys = ["dataset_id", "cell_type_ontology_term_id", "tissue_general_ontology_term_id", "assay_ontology_term_id"]
        groups = [(k, g) for k, g in obs.groupby(keys, observed=True) if len(g) >= a.min_cells]
        groups.sort(key=lambda kg: -len(kg[1]))
        fids = np.asarray([f["soma_joinid"] for f in feats], dtype=np.int64)
        for gi, (k, g) in enumerate(groups):
            sid = hashlib.sha1("|".join(k).encode()).hexdigest()[:16]
            done = out / "strata" / f"{sid}.json"
            if done.exists():
                continue
            t0 = time.time()
            acc = tc.Accumulator(feats, metacell_size=a.metacell_size)
            g = g.sort_values("soma_joinid")
            for s in range(0, len(g), a.batch_size):
                b = g.iloc[s:s + a.batch_size]
                ad = cc.get_anndata(census, organism="Homo sapiens", measurement_name="RNA", X_name="raw",
                                    obs_coords=b.soma_joinid.to_numpy(np.int64), var_coords=fids,
                                    obs_column_names=["soma_joinid"], var_column_names=["soma_joinid"])
                opos = {int(j): i for i, j in enumerate(ad.obs.soma_joinid)}
                vpos = {int(j): i for i, j in enumerate(ad.var.soma_joinid)}
                x = ad.X[[opos[int(j)] for j in b.soma_joinid]][:, [vpos[int(j)] for j in fids]]
                acc.add(x, b.raw_sum.to_numpy(), b.donor_id.astype(str).to_numpy())
            n = acc.n
            frac = acc.n_detected / n
            mc_r, mc_df, n_mc = acc.metacell_corr()
            stat = "metacell"
            r_mat, df = mc_r, mc_df
            if r_mat is None:
                r_mat, df, _ = acc.pseudobulk_corr(); stat = "pseudobulk"
            rows = []
            for (ga, gb) in sorted(pairs):
                if ga not in fidx or gb not in fidx:
                    continue
                i, j = fidx[ga], fidx[gb]
                if min(frac[i], frac[j]) < a.min_detection:
                    continue
                cod = int(acc.codetect[i, j])
                lift = cod * n / (acc.n_detected[i] * acc.n_detected[j])
                rows.append({"tf1": ga, "tf2": gb, "co_detected": cod, "lift": float(lift),
                             "det1": float(frac[i]), "det2": float(frac[j]),
                             "r": None if r_mat is None else float(r_mat[i, j])})
            if r_mat is not None and rows:
                p = tc.corr_pvalues(np.array([x["r"] for x in rows]), df)
                q = tc.bh_qvalues(p)
                for x, pv, qv in zip(rows, p, q):
                    x["p"], x["q"] = float(pv), float(qv)
                    x["hit"] = bool(x["r"] >= a.min_rho and qv <= a.max_q and x["lift"] >= a.min_lift
                                    and x["co_detected"] >= a.min_codetected)
            dmeta = meta.loc[k[0]] if k[0] in meta.index else {}
            res = {"dataset_id": k[0], "collection_name": str(dmeta.get("collection_name", "")),
                   "dataset_title": str(dmeta.get("dataset_title", "")), "cell_type_id": k[1],
                   "cell_type": str(g.cell_type.iloc[0]), "tissue_id": k[2], "tissue": str(g.tissue_general.iloc[0]),
                   "assay_id": k[3], "assay": str(g.assay.iloc[0]), "n_cells": int(n), "n_donors": int(g.donor_id.nunique()),
                   "n_metacells": int(n_mc), "stat": stat if r_mat is not None else "none", "df": df,
                   "n_tested": len(rows), "pairs": rows}
            tc.write_json_atomic(done, res)
            print(f"[{gi+1}/{len(groups)}] {k[0][:8]} {res['cell_type']}/{res['tissue']}/{res['assay']}: {n} cells,"
                  f" {res['n_donors']} donors, {len(rows)} tested, {sum(x.get('hit', False) for x in rows)} hits ({time.time()-t0:.0f}s)", flush=True)
    tc.write_json_atomic(out / "_SUCCESS", {"n_strata": len(groups)})
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
