#!/usr/bin/env python3
"""Generate a tiny synthetic dataset that exercises every pipeline stage offline (minutes on a laptop).

Writes into OUT (default: $TFPM_V2P_SMOKE or /tmp/tfpm_v2p_smoke/data):
  genome.fa.gz              synthetic chr21/chr22 (600 kb each) with planted motif pairs in 'peaks'
  blacklist.bed.gz          empty blacklist
  ccla/SMOKE_rep{1,2}.fragments.tsv.gz      CCLA-style line, 2 replicates (HNF4A+FOXA2 peaks)
  corpus/SMK1.h5ad + SMK1_fragments.tar.gz  Corpus-style sample: Hepatocyte (HNF4A+FOXA2) and Macrophage (SPI1+CEBPB)
  smoke.chrombpnet_nobias.h5                tiny Keras model with the ChromBPNet I/O (2114x4 -> 1000 profile, 1 count)
  rna_census_like.h5ad                      CELLxGENE-schema cells (2 datasets incl. Tabula Sapiens) with correlated TF pairs
  samples.tsv, context_map.tsv              manifest + context map in the production format
Requires numpy, pandas, anndata, scipy (+ tensorflow for the model; skipped with --no-model).
"""
import argparse, gzip, hashlib, io, json, os, re, sys, tarfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MEME = ROOT / "mining/motif_cooccurrence/resources/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
PAIRS = {"liver": ("MA0114.5", "MA0047.4"), "myeloid": ("MA0080.7", "MA0466.4")}   # HNF4A+FOXA2, Spi1+CEBPB
ENSG = {"HNF4A": "ENSG00000101076", "FOXA2": "ENSG00000125798", "SPI1": "ENSG00000066336", "CEBPB": "ENSG00000172216",
        "GATA1": "ENSG00000102145", "KLF1": "ENSG00000105610"}
TABULA = "53d208b0-2cfd-4366-9866-c3c6114081bc"
CHR_LEN = 600_000


def consensus(mid):
    lines = MEME.read_text().splitlines(); i = next(k for k, l in enumerate(lines) if l.startswith(f"MOTIF {mid}"))
    rows = []
    for l in lines[i + 2:]:
        p = l.split()
        if not p or l.startswith("URL") or l.startswith("MOTIF"):
            if rows: break
            continue
        try: rows.append([float(x) for x in p])
        except ValueError: continue
    return "".join("ACGT"[int(np.argmax(r))] for r in rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.environ.get("TFPM_V2P_SMOKE", "/tmp/tfpm_v2p_smoke/data"))
    ap.add_argument("--no-model", action="store_true"); a = ap.parse_args()
    out = Path(a.out); (out / "ccla").mkdir(parents=True, exist_ok=True); (out / "corpus").mkdir(exist_ok=True)
    rng = np.random.default_rng(1)
    cons = {k: (consensus(m1), consensus(m2)) for k, (m1, m2) in PAIRS.items()}
    print("planted:", cons)
    # ---- genome + peaks
    seqs, peaks = {}, {"liver": [], "myeloid": [], "shared": []}
    for c in ("chr21", "chr22"):
        s = np.array(list("ACGT"))[rng.choice(4, CHR_LEN, p=[.29, .21, .21, .29])]
        for k, centre in enumerate(range(5000, CHR_LEN - 5000, 2000)):
            kind = ("liver", "myeloid", "shared")[k % 3]
            if kind != "shared" and rng.random() < 0.85:
                m1, m2 = cons[kind]; gap = int(rng.integers(8, 40))
                p1 = centre - (len(m1) + gap + len(m2)) // 2; p2 = p1 + len(m1) + gap
                s[p1:p1 + len(m1)] = list(m1); s[p2:p2 + len(m2)] = list(m2)
            peaks[kind].append((c, centre))
        seqs[c] = "".join(s)
    with gzip.open(out / "genome.fa.gz", "wt") as f:
        for c, s in seqs.items():
            f.write(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    with gzip.open(out / "blacklist.bed.gz", "wt") as f:
        f.write("")
    # ---- CCLA-style fragments (liver-like line)
    def frags(peak_list, n_per_peak, n_bg, barcodes):
        rows = []
        for c, centre in peak_list:
            for _ in range(rng.poisson(n_per_peak)):
                st = int(centre + rng.normal(0, 60)); ln = int(rng.integers(60, 250))
                rows.append((c, st, st + ln, barcodes[rng.integers(len(barcodes))]))
        for _ in range(n_bg):
            c = ("chr21", "chr22")[rng.integers(2)]; st = int(rng.integers(1000, CHR_LEN - 1000))
            rows.append((c, st, st + int(rng.integers(60, 250)), barcodes[rng.integers(len(barcodes))]))
        rows.sort(key=lambda r: (r[0], r[1])); return rows
    for rep in (1, 2):
        bcs = [f"L{rep}_{i:04d}" for i in range(500)]
        rows = frags(peaks["liver"] + peaks["shared"], 30, 20000, bcs)
        with gzip.open(out / "ccla" / f"SMOKE_rep{rep}.fragments.tsv.gz", "wt") as f:
            f.writelines(f"{c}\t{s}\t{e}\t{b}\t1\n" for c, s, e, b in rows)
    # ---- Corpus-style sample
    import anndata as ad, pandas as pd, scipy.sparse as sp
    allp = sorted(peaks["liver"] + peaks["myeloid"] + peaks["shared"])
    var = [f"{c}:{x - 250}-{x + 250}" for c, x in allp]
    kinds = {p: k for k in peaks for p in peaks[k]}
    cts = ["Hepatocyte"] * 180 + ["Macrophage"] * 170
    obs_names = [f"SMK1#{'ACGT'[i % 4]}{i:06d}-1" for i in range(len(cts))]
    X = np.zeros((len(cts), len(allp)), dtype=np.float32); frag_rows = []
    for i, ct in enumerate(cts):
        own = "liver" if ct == "Hepatocyte" else "myeloid"
        pr = np.array([0.35 if kinds[p] == own else 0.25 if kinds[p] == "shared" else 0.02 for p in allp])
        hit = rng.random(len(allp)) < pr; X[i, hit] = 1
        bc = obs_names[i].split("#")[1]
        for j in np.flatnonzero(hit):
            c, x = allp[j]; st = int(x + rng.normal(0, 60)); frag_rows.append((c, st, st + int(rng.integers(60, 250)), bc))
    A = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"Cell_type (HSC)": cts}, index=obs_names),
                   var=pd.DataFrame(index=var))
    A.write_h5ad(out / "corpus" / "SMK1.h5ad")
    frag_rows.sort(key=lambda r: (r[0], r[1])); buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as g:
        g.write("".join(f"{c}\t{s}\t{e}\t{b}\t1\n" for c, s, e, b in frag_rows).encode())
    with tarfile.open(out / "corpus" / "SMK1_fragments.tar.gz", "w:gz") as t:
        ti = tarfile.TarInfo("SMK1/SMK1_fragments.tsv.gz"); data = buf.getvalue(); ti.size = len(data)
        t.addfile(ti, io.BytesIO(data))
    # ---- RNA (CELLxGENE schema subset)
    genes = list(ENSG.values()) + [f"ENSG9999900{k:04d}" for k in range(40)]
    rows, mats = [], []
    for ds, coll, title in (("smoke-ds-1", "Smoke Liver+Blood Atlas", "smoke dataset 1"),
                            (TABULA, "Tabula Sapiens", "Tabula Sapiens - All Cells")):
        for cl, lab, (g1, g2) in (("CL:0000182", "hepatocyte", ("HNF4A", "FOXA2")), ("CL:0000235", "macrophage", ("SPI1", "CEBPB"))):
            for d in range(3):
                n = 120; z = np.repeat(rng.gamma(2.0, 1.0, n // 20), 20)   # program shared by metacell-sized blocks
                M = rng.poisson(0.3, (n, len(genes))).astype(np.float32)
                for gname in (g1, g2):
                    M[:, genes.index(ENSG[gname])] = rng.poisson(0.4 + 1.2 * z)
                mats.append(M)
                rows += [{"dataset_id": ds, "collection_name": coll, "dataset_title": title, "donor_id": f"{ds[:6]}_d{d}",
                          "cell_type_ontology_term_id": cl, "cell_type": lab, "tissue_ontology_term_id": "UBERON:0002107",
                          "tissue": "liver", "tissue_general": "liver", "disease_ontology_term_id": "PATO:0000461",
                          "disease": "normal", "assay": "10x 3' v3", "assay_ontology_term_id": "EFO:0009922",
                          "is_primary_data": True, "suspension_type": "cell", "organism": "Homo sapiens"} for _ in range(n)]
    R = ad.AnnData(X=sp.csr_matrix(np.vstack(mats)), obs=pd.DataFrame(rows, index=[f"c{i}" for i in range(len(rows))]),
                   var=pd.DataFrame({"feature_name": [next((k for k, v in ENSG.items() if v == g), g) for g in genes]}, index=genes))
    R.var.index.name = "feature_id"; R.write_h5ad(out / "rna_census_like.h5ad")
    # ---- tiny ChromBPNet-shaped model
    if not a.no_model:
        import tensorflow as tf
        inp = tf.keras.Input((2114, 4), name="sequence")
        x = tf.keras.layers.Conv1D(8, 21, padding="valid", activation="relu")(inp)
        x = tf.keras.layers.Cropping1D((x.shape[1] - 1000) // 2)(x); x = tf.keras.layers.Cropping1D((0, x.shape[1] - 1000))(x)
        prof = tf.keras.layers.Flatten(name="logits_profile_predictions")(tf.keras.layers.Conv1D(1, 1)(x))
        cnt = tf.keras.layers.Dense(1, name="logcount_predictions")(tf.keras.layers.GlobalAveragePooling1D()(x))
        tf.keras.Model(inp, [prof, cnt]).save(out / "smoke.chrombpnet_nobias.h5")
    # ---- manifest + context map
    def sz(p): return str(Path(p).stat().st_size)
    def md5(p): return hashlib.md5(Path(p).read_bytes()).hexdigest()
    H = ("profile_id atlas source kind context tissue inputs sizes md5 checksum_note fragments fragment_sizes celltype_col "
         "n_cells n_cell_types accession pmid cbp_pretrained cbp_pretrained_sha256 rna_context mapping_level").split()
    c1, c2 = out / "ccla/SMOKE_rep1.fragments.tsv.gz", out / "ccla/SMOKE_rep2.fragments.tsv.gz"
    h5, tg = out / "corpus/SMK1.h5ad", out / "corpus/SMK1_fragments.tar.gz"
    rows = [{"profile_id": "CCLA_SMOKE", "atlas": "CCLA (scatlas HCL cell lines)", "source": "synthetic", "kind": "fragments",
             "context": "SMOKE liver-like line", "tissue": "liver", "inputs": f"{c1},{c2}", "sizes": f"{sz(c1)},{sz(c2)}",
             "md5": f"{md5(c1)},{md5(c2)}", "checksum_note": "synthetic", "fragments": "", "fragment_sizes": "",
             "celltype_col": "", "n_cells": "1000", "n_cell_types": "1", "accession": "smoke", "pmid": "",
             "cbp_pretrained": "", "cbp_pretrained_sha256": "", "rna_context": "hepatocyte", "mapping_level": "proxy"},
            {"profile_id": "SMK1", "atlas": "Human scATAC Corpus", "source": "synthetic", "kind": "corpus",
             "context": "SMOKE liver+macrophage tissue", "tissue": "liver", "inputs": str(h5), "sizes": sz(h5),
             "md5": md5(h5), "checksum_note": "synthetic", "fragments": str(tg), "fragment_sizes": sz(tg),
             "celltype_col": "Cell_type (HSC)", "n_cells": str(len(cts)), "n_cell_types": "2", "accession": "smoke",
             "pmid": "", "cbp_pretrained": "", "cbp_pretrained_sha256": "", "rna_context": "auto", "mapping_level": "auto"}]
    with open(out / "samples.tsv", "w") as f:
        f.write("\t".join(H) + "\n"); f.writelines("\t".join(r[h] for h in H) + "\n" for r in rows)
    CM = [("profile", "CCLA_SMOKE", "", "hepatocyte", "CL:0000182", "UBERON:0002107", "PATO:0000461", "proxy", "smoke line"),
          ("celltype", "Hepatocytes?", "", "hepatocyte", "CL:0000182", "", "PATO:0000461", "exact", "smoke"),
          ("celltype", "Macrophages?", "", "macrophage", "CL:0000235", "", "PATO:0000461", "exact", "smoke"),
          ("celltype", ".*", "", "", "", "", "", "unmapped", "default")]
    with open(out / "context_map.tsv", "w") as f:
        f.write("match_on\tpattern\tsample_pattern\trna_context\tcl_ids\ttissue_id\tdisease_id\tmapping_level\tnote\n")
        f.writelines("\t".join(r) + "\n" for r in CM)
    print("smoke data in", out)


if __name__ == "__main__":
    main()
