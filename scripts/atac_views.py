#!/usr/bin/env python3
"""Export CATlas (scATAC level B) views for the site: per-pair cross-cluster co-enrichment over all CATlas clusters,
motif hit fractions with CIs (mapped cluster), JASPAR PFMs + motif similarity, and optional within-cCRE motif spacing.

Reuses mining/catlas_pipeline.py (matched_indices, motif_presence, odds_p, bh, parse_meme) on the compiled motif hits of a
completed `tfpm atac` run, so it is the same matched accessible-vs-inaccessible distal cCRE test, run across all clusters.
Output is sharded: <site_run>/atac/index.json, atac/motifs.json, atac/pairs/<A>__<B>.json. Per-cluster matching is cached
in <atac_dir>/views_cache/ so reruns and full runs only pay the matching once per cluster (shard with --cluster-shards/-index).
Same-cCRE motif co-enrichment only: not co-binding and not a TF-TF interaction test."""
import argparse, hashlib, json, math, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "mining"))
import catlas_pipeline as cp  # noqa: E402

CAVEAT = ("Same-cCRE JASPAR motif co-enrichment in accessible vs GC/width-matched inaccessible distal CATlas cCREs. "
          "Not co-binding, not a TF-TF interaction test, no target-gene assignment. Similar motifs (e.g. AP-1 family) inflate joint hits.")


def wilson(k, n, z=1.96):
    if not n:
        return [None, None, None]
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(p, 5), round(max(0, c - h), 5), round(min(1, c + h), 5)]


def pfm(block):
    rows = []
    for line in block.splitlines():
        parts = line.split()
        if len(parts) == 4:
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                pass
    return np.array(rows)


def similarity(a, b, min_overlap=5):
    """Max over offsets and strands of the mean per-column Pearson r (overlap >= min_overlap); 1 = identical motif."""
    best = -1.0
    for bb in (b, b[::-1, ::-1]):  # reverse complement: reverse columns order, ACGT -> TGCA
        for off in range(-(len(bb) - min_overlap), len(a) - min_overlap + 1):
            ia = range(max(0, off), min(len(a), off + len(bb)))
            if len(ia) < min_overlap:
                continue
            rs = []
            for i in ia:
                x, y = a[i] - a[i].mean(), bb[i - off] - bb[i - off].mean()
                den = math.sqrt((x * x).sum() * (y * y).sum())
                rs.append((x * y).sum() / den if den else 0.0)
            best = max(best, float(np.mean(rs)))
    return round(best, 3)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True); p.add_argument("--site-run-dir", required=True)
    p.add_argument("--catlas-h5ad", required=True)
    p.add_argument("--atac-dir", help="atac dir with compiled/ hits (default: atac_catlas if scanned, else atac_catlas_proxy_fetal)")
    p.add_argument("--genome-fasta"); p.add_argument("--fimo", default="fimo"); p.add_argument("--bedtools", default="bedtools")
    p.add_argument("--min-matched", type=int, default=1000); p.add_argument("--spacing-max-peaks", type=int, default=3000)
    a = p.parse_args()
    run = Path(a.run); site = Path(a.site_run_dir) / "atac"
    cands = [Path(a.atac_dir)] if a.atac_dir else [run / "atac_catlas", run / "atac_catlas_proxy_fetal"]
    src = next((d for d in cands if (d / "compiled" / "hits.npy").exists()), None)
    if src is None:
        print("No compiled CATlas motif scan in this run; nothing to export (coverage-gap runs have no scan)."); return
    manifest, plan = cp.load_plan(src)
    tier = "proxy_fetal" if src.name.endswith("proxy_fetal") else "exact"
    # mapped clusters per RNA cell type (from plans of both atac dirs)
    mapping = defaultdict(list)
    for d, t in ((run / "atac_catlas", "exact"), (run / "atac_catlas_proxy_fetal", "proxy_fetal")):
        if (d / "plan.tsv").exists():
            for r in cp.read_tsv(d / "plan.tsv"):
                if r["catlas_cluster"]:
                    m = {"cluster": r["catlas_cluster"], "tier": r["mapping_tier"], "stage": r["catlas_stage"], "run": t}
                    if m not in mapping[r["cell_type_id"]]:
                        mapping[r["cell_type_id"]].append(m)
    # unique TF pairs with direct motifs
    pairs = {}
    for r in plan:
        if r["motif_status"] != "both_exact_motifs":
            continue
        k = tuple(sorted([(r["tf1_symbol"], r["tf1_motifs"]), (r["tf2_symbol"], r["tf2_motifs"])]))
        pairs.setdefault(k, set()).add(r["cell_type_id"])
    import anndata as ad
    h5 = ad.read_h5ad(a.catlas_h5ad, backed="r")
    names = list(h5.obs["cell type"].astype(str)); stages = ["fetal" if n.lower().startswith("fetal ") else "adult" for n in names]
    width = (h5.var["end"].to_numpy() - h5.var["start"].to_numpy()).astype(np.int32)
    distal = h5.var["Class"].astype(str).str.lower().eq("distal").to_numpy()
    packed = np.load(src / "compiled" / "hits.npy", mmap_mode="r"); gc = np.load(src / "compiled" / "gc.npy", mmap_mode="r")
    midx = {r["motif_id"]: int(r["row_index"]) for r in cp.read_tsv(src / "compiled" / "motif_ids.tsv")}
    eligible = distal & np.isfinite(gc) & (width > 0); n = manifest["n_peaks"]
    pres = {}
    for k in pairs:
        for sym, mot in k:
            pres.setdefault(sym, cp.motif_presence(mot, packed, midx, n))
    cache = src / "views_cache"; cache.mkdir(exist_ok=True)
    X = h5.X
    stats = {k: [] for k in pairs}
    for ci, cl in enumerate(names):
        cf = cache / f"{hashlib.sha256(cl.encode()).hexdigest()[:16]}.npz"
        if cf.exists():
            z = np.load(cf); fg, bg, n_acc = z["fg"], z["bg"], int(z["n_acc"])
        else:
            row = X[ci:ci + 1]; acc = np.zeros(n, bool); acc[row.indices if hasattr(row, "indices") else np.flatnonzero(np.asarray(row.todense()).ravel())] = True
            n_acc = int((acc & eligible).sum())
            seed = int(hashlib.sha256((cl + "|CATlas").encode()).hexdigest()[:16], 16)   # same seed as the pipeline
            fg, bg = cp.matched_indices(acc, eligible, gc, width, seed)
            np.savez(cf, fg=fg, bg=bg, n_acc=n_acc)
        for k in pairs:
            A, B = pres[k[0][0]], pres[k[1][0]]
            af, ab, bf, bb = int(A[fg].sum()), int(A[bg].sum()), int(B[fg].sum()), int(B[bg].sum())
            J = A & B; jf, jb = int(J[fg].sum()), int(J[bg].sum()); nf = len(fg)
            s = {"c": ci, "n_acc": n_acc, "n": nf}
            if nf >= a.min_matched and min(af, bf) >= 20:
                s.update(or_a=cp.odds_p(af, nf, ab, nf)[0], or_b=cp.odds_p(bf, nf, bb, nf)[0])
                s["or_j"], s["p_j"] = cp.odds_p(jf, nf, jb, nf)
                s.update(lift_fg=jf * nf / (af * bf) if af and bf else 0.0, lift_bg=jb * nf / (ab * bb) if ab and bb else 0.0,
                         fa=wilson(af, nf), fb=wilson(bf, nf), fj=wilson(jf, nf), ba=wilson(ab, nf), bb_=wilson(bb, nf), bj=wilson(jb, nf),
                         jf=jf, jb=jb, acc_w=round(jf / max(1, int((J & eligible).sum())), 4))
            stats[k].append(s)
        if ci % 20 == 0:
            print(f"cluster {ci + 1}/{len(names)}", flush=True)
    h5.file.close()
    # motifs + similarity
    _, blocks = cp.parse_meme(src / "selected_motifs.meme")
    motifs = {mid: {"name": nm, "pfm": [[round(x, 3) for x in r] for r in pfm(bl).tolist()]} for mid, (nm, bl) in blocks.items()}
    mats = {mid: np.array(m["pfm"]) for mid, m in motifs.items()}
    sim = {f"{x}|{y}": similarity(mats[x], mats[y]) for x in mats for y in mats if x < y}
    site.mkdir(parents=True, exist_ok=True); (site / "pairs").mkdir(exist_ok=True)
    cp.write_json(site / "motifs.json", {"motifs": motifs, "similarity": sim, "source": "JASPAR 2026 CORE vertebrates (non-redundant)",
                                         "similarity_method": "max mean column Pearson r over offsets/strands (overlap>=5)"})
    # optional spacing in the mapped cluster's joint accessible cCREs
    def spacing(k, cl_index):
        if not a.genome_fasta or cl_index is None:
            return None
        cf = cache / f"{hashlib.sha256(names[cl_index].encode()).hexdigest()[:16]}.npz"; fg = np.load(cf)["fg"]
        A, B = pres[k[0][0]], pres[k[1][0]]; peaks = fg[(A & B)[fg]][: a.spacing_max_peaks]
        if not len(peaks):
            return None
        import anndata as ad2
        v = ad2.read_h5ad(a.catlas_h5ad, backed="r").var
        mids = {k[0][0]: set(k[0][1].split(",")), k[1][0]: set(k[1][1].split(","))}
        with tempfile.TemporaryDirectory() as t:
            bed, fa, meme = Path(t) / "p.bed", Path(t) / "p.fa", Path(t) / "m.meme"
            with open(bed, "w") as fh:
                for pi in peaks:
                    fh.write(f"{v['chrom'].iloc[pi]}\t{int(v['start'].iloc[pi])}\t{int(v['end'].iloc[pi])}\tpeak_{pi}\n")
            with open(fa, "w") as fh:
                subprocess.run([a.bedtools, "getfasta", "-fi", a.genome_fasta, "-bed", str(bed), "-nameOnly"], stdout=fh, check=True)
            hdr, _ = cp.parse_meme(src / "selected_motifs.meme")
            meme.write_text(hdr + "".join(blocks[m][1] for m in sorted(mids[k[0][0]] | mids[k[1][0]])))
            out = subprocess.run([a.fimo, "--text", "--no-qvalue", "--thresh", "1e-4", str(meme), str(fa)], capture_output=True, text=True, check=True).stdout
        sites = defaultdict(lambda: defaultdict(list))
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) < 6 or f[0] == "motif_id" or line.startswith("#"):
                continue
            mid, pk, st, en = f[0], f[2], int(f[3]), int(f[4])
            for sym, ms in mids.items():
                if mid in ms:
                    sites[pk][sym].append((st + en) / 2)
        dists, posA, posB, overl = [], [], [], 0
        sa, sb = k[0][0], k[1][0]
        for pk, d in sites.items():
            if d[sa] and d[sb]:
                dd = min(abs(x - y) for x in d[sa] for y in d[sb]); dists.append(dd); overl += dd < 4
                posA += [x - 200 for x in d[sa]]; posB += [x - 200 for x in d[sb]]
        h = lambda v, lo, hi, w: np.histogram(v, bins=np.arange(lo, hi + w, w))[0].tolist()
        return {"cluster": names[cl_index], "n_peaks": int(len(peaks)), "n_with_both": len(dists),
                "frac_overlapping_sites": round(overl / max(1, len(dists)), 3),
                "dist_bins": list(range(0, 201, 5)), "dist_hist": h(dists, 0, 200, 5),
                "pos_bins": list(range(-200, 201, 10)), "pos_a": h(posA, -200, 200, 10), "pos_b": h(posB, -200, 200, 10)}
    idx_pairs = []
    for k, cts in pairs.items():
        S = stats[k]; tested = [s for s in S if "p_j" in s]
        for s, q in zip(tested, cp.bh([s["p_j"] for s in tested])):
            s["q_j"] = q
        for s in S:
            for f in ("or_a", "or_b", "or_j", "p_j", "q_j", "lift_fg", "lift_bg"):
                if f in s and s[f] is not None and math.isfinite(s[f]):
                    s[f] = float(f"{s[f]:.4g}")
                elif f in s:
                    s[f] = None
        mapped = [m for ct in cts for m in mapping.get(ct, [])]
        mi = next((names.index(m["cluster"]) for m in mapped if m["cluster"] in names), None)
        (sa, ma), (sb, mb) = k
        sims = [sim.get("|".join(sorted([x, y])), 1.0 if x == y else None) for x in ma.split(",") for y in mb.split(",")]
        fn = f"{sa}__{sb}.json"
        cp.write_json(site / "pairs" / fn, {"tf_a": sa, "tf_b": sb, "motifs_a": ma.split(","), "motifs_b": mb.split(","),
                      "motif_similarity_max": max([x for x in sims if x is not None], default=None),
                      "rna_cell_types": sorted(cts), "mapped": mapped, "clusters": S, "spacing": spacing(k, mi),
                      "scan_source": src.name, "scan_tier": tier, "caveat": CAVEAT})
        idx_pairs.append({"a": sa, "b": sb, "file": fn})
    cp.write_json(site / "index.json", {"clusters": names, "stages": stages, "mapping": mapping, "pairs": idx_pairs,
                  "scan_source": src.name, "scan_tier": tier, "caveat": CAVEAT, "min_matched": a.min_matched,
                  "note": "Cross-cluster stats re-run the pipeline's matched test for every CATlas cluster. Only the RNA-mapped cluster "
                          "(exact or labelled proxy) bears on the RNA context; other clusters are context, not validation."})
    print(f"Exported {len(idx_pairs)} pairs x {len(names)} clusters -> {site}")


if __name__ == "__main__":
    main()
