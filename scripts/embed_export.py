#!/usr/bin/env python3
"""Optional per-cell layer for the site: UMAP + per-cell TF expression for datasets in a published run.

For each dataset in <run>/datasets.tsv that has a completed output (or --dataset-id), it:
  1. reads primary-data cells from the pinned Census (obs), downsampled to --max-cells, stratified by cell type;
  2. reads raw counts for all mapped compendium genes (features.tsv), log1p(1e4*count/raw_sum);
  3. computes PCA(30) -> kNN(15) -> UMAP with scanpy on those genes (Census ships no embeddings
     for this release through the 1.18 API, so the embedding is computed here and labelled as such);
  4. writes <site>/runs/<run_id>/embed/<dataset_id>/cells.json and genes/<SYMBOL>.json (sparse per-cell
     log1p CPTT for TFs in candidate pairs + positive controls, capped by --max-genes),
  and records the layer in run.json. Real Census data only; no imputation.
"""
import argparse, csv, json, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "mining"))
import tf_census as m  # noqa: E402

safe = lambda s: re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True)
    p.add_argument("--site-data", default="docs/data")
    p.add_argument("--run-id", required=True, help="published run id (docs/data/runs/<run_id>)")
    p.add_argument("--dataset-id", action="append", default=[])
    p.add_argument("--max-datasets", type=int, default=20)
    p.add_argument("--max-cells", type=int, default=8000, help="per-dataset cell cap (size control)")
    p.add_argument("--max-genes", type=int, default=300)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    import scanpy as sc
    import anndata as ad
    from scipy import sparse
    run = Path(a.run); out = Path(a.site_data) / "runs" / a.run_id
    info_path = out / "run.json"
    info = json.loads(info_path.read_text())
    manifest = json.loads((run / "manifest.json").read_text())
    feats = list(m.read_tsv(run / "features.tsv"))
    sym2feat = {f["symbol"]: f for f in feats}
    genes = []
    cand = run / "aggregate" / "cell_type_pair_candidates.tsv"
    if cand.exists():
        for r in m.read_tsv(cand):
            for s in (r["tf1_symbol"], r["tf2_symbol"]):
                if s not in genes:
                    genes.append(s)
    for a_, b_ in m.POSITIVE_CONTROLS:
        for s in (a_, b_):
            if s not in genes and s in sym2feat:
                genes.append(s)
    genes = [g for g in genes if g in sym2feat][:a.max_genes]
    strata = json.loads((out / "strata.json").read_text())
    skey = {tuple(r[:8]): i for i, r in enumerate(strata["rows"])}
    datasets = list(m.read_tsv(run / "datasets.tsv"))
    if a.dataset_id:
        datasets = [d for d in datasets if d["dataset_id"] in a.dataset_id]
    else:
        datasets = [d for d in datasets if (run / "datasets" / d["dataset_id"] / "_SUCCESS").exists()][:a.max_datasets]
    rng = np.random.default_rng(a.seed)
    layer = []
    api = m.census_module()
    with api.open_soma(census_version=manifest["resolved_census_version"]) as census:
        for d in datasets:
            did = d["dataset_id"]
            obs = census["census_data"]["homo_sapiens"].obs.read(
                value_filter=f"is_primary_data == True and dataset_id == '{did}'",
                column_names=["soma_joinid", "cell_type_ontology_term_id", "cell_type", "tissue_general_ontology_term_id", "tissue_general",
                              "disease_ontology_term_id", "disease", "assay_ontology_term_id", "assay", "raw_sum"]).concat().to_pandas()
            obs = obs[obs.raw_sum > 0]
            n_total = len(obs)
            if n_total < 30:
                continue
            if n_total > a.max_cells:   # stratified by cell type, >= min(50, n) per type
                keep = []
                frac = a.max_cells / n_total
                for _, g in obs.groupby("cell_type"):
                    k = min(len(g), max(50, int(round(len(g) * frac))))
                    keep.append(g.sample(n=k, random_state=a.seed))
                obs = __import__("pandas").concat(keep)
                if len(obs) > a.max_cells:
                    obs = obs.sample(n=a.max_cells, random_state=a.seed)
            obs = obs.sort_values("soma_joinid")
            var_ids = np.array(sorted(int(f["soma_joinid"]) for f in feats), dtype=np.int64)
            adata = api.get_anndata(census, organism="Homo sapiens", measurement_name="RNA", X_name="raw",
                                    obs_coords=obs.soma_joinid.to_numpy(), var_coords=var_ids,
                                    obs_column_names=["soma_joinid"], var_column_names=["soma_joinid", "feature_id"])
            adata.obs.index = adata.obs.soma_joinid.astype(str)
            o = obs.set_index(obs.soma_joinid.astype(str)).loc[adata.obs.index]
            X = sparse.csr_matrix(adata.X, dtype=np.float32)
            X = sparse.csr_matrix(X.multiply((1e4 / o.raw_sum.to_numpy())[:, None]))
            X.data = np.log1p(X.data)
            emb = ad.AnnData(X=X)
            emb.var_names = adata.var["feature_id"].astype(str).to_numpy()
            # keep only genes measured in this dataset (presence matrix) to avoid structural zeros
            present = m.get_present_feature_ids(census, int(d["soma_joinid"]), [{"soma_joinid": int(v)} for v in var_ids])
            pres_mask = np.isin(adata.var["soma_joinid"].astype(int).to_numpy(), list(present))
            e2 = emb[:, pres_mask].copy()
            sc.pp.highly_variable_genes(e2, n_top_genes=min(2000, e2.n_vars), flavor="seurat")
            e2 = e2[:, e2.var.highly_variable].copy()
            sc.pp.scale(e2, max_value=10)
            sc.tl.pca(e2, n_comps=min(30, e2.n_vars - 1, e2.n_obs - 1), random_state=a.seed)
            sc.pp.neighbors(e2, n_neighbors=15, random_state=a.seed)
            sc.tl.umap(e2, random_state=a.seed)
            xy = e2.obsm["X_umap"]
            cts = sorted(o.cell_type.unique())
            cti = {c: i for i, c in enumerate(cts)}
            sidx = [skey.get(tuple(r), -1) for r in o[["cell_type_ontology_term_id", "cell_type", "tissue_general_ontology_term_id", "tissue_general",
                                                        "disease_ontology_term_id", "disease", "assay_ontology_term_id", "assay"]].itertuples(index=False)]
            dout = out / "embed" / did
            (dout / "genes").mkdir(parents=True, exist_ok=True)
            (dout / "cells.json").write_text(json.dumps({
                "dataset_id": did, "dataset_title": d.get("dataset_title", ""), "collection_name": d.get("collection_name", ""),
                "n_primary_cells": n_total, "n_cells": len(o), "method": "scanpy: log1p CPTT on Census-present compendium genes, HVG(<=2000) -> scale -> PCA(30) -> kNN(15) -> UMAP",
                "cell_types": cts, "cols": ["x", "y", "ct", "s"],
                "rows": [[round(float(x), 3), round(float(y), 3), cti[c], s] for (x, y), c, s in zip(xy, o.cell_type, sidx)]}, separators=(",", ":")))
            col = {fid: i for i, fid in enumerate(adata.var["feature_id"].astype(str))}
            Xc = X.tocsc()
            written = []
            for g in genes:
                j = col.get(sym2feat[g]["ensembl_id"])
                if j is None:
                    continue
                colv = Xc[:, j]
                nz = colv.nonzero()[0]
                (dout / "genes" / f"{safe(g)}.json").write_text(json.dumps({"gene": g, "measured": int(sym2feat[g]["soma_joinid"]) in present,
                    "i": nz.tolist(), "v": [round(float(v), 2) for v in colv[nz].toarray().ravel()]}, separators=(",", ":")))
                written.append(g)
            layer.append({"dataset_id": did, "title": d.get("dataset_title", ""), "n_cells": len(o), "n_primary_cells": n_total, "genes": written})
            print(f"{did}: {len(o)}/{n_total} cells, {len(written)} genes", flush=True)
    info["embedding"] = {"status": "done" if layer else "not_run", "datasets": layer, "max_cells": a.max_cells}
    info_path.write_text(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
