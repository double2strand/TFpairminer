#!/usr/bin/env python3
"""Join ChromBPNet pairs_cbp.tsv with two targeted RNA screens (targeted_rna.py output).

Same four evidence views as integrate.py, but contexts are cell-line -> normal-tissue PROXY
matches (labelled in every row), and RNA is targeted rather than atlas-wide. ATAC atlases stay
separate; anything not run is 'unassessed', never negative.
"""
from __future__ import annotations

import argparse, csv, json
from collections import defaultdict
from pathlib import Path

from .integrate import module_rows

VIEWS = {"atac_only": "01_atac_only.tsv", "atac_rna1": "02_atac_plus_rna1.tsv",
         "atac_rna2": "03_atac_plus_rna2.tsv", "rna_only": "04_rna_replicated_cbp_no_support.tsv"}
FIELDS = ["atac_dataset", "atac_context", "rna_context", "mapping_level", "mapping_note", "k", "tfs", "motifs",
          "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol", "atac_status", "cbp_median_site_pct",
          "cbp_synergy", "cbp_synergy_q", "cbp_cooperative",
          "rna1_status", "rna1_best_r", "rna1_min_q", "rna1_n_hit_strata", "rna1_n_tested_strata", "rna1_best_stratum",
          "rna2_status", "rna2_best_r", "rna2_min_q", "rna2_n_hit_strata", "rna2_n_tested_strata", "rna2_best_stratum",
          "evidence_view", "note", "cell_type_id", "tissue_id", "disease_id"]


def load_screen(root: Path):
    man = json.loads((root / "manifest.json").read_text())
    if not (root / "_SUCCESS").exists():
        return man, None
    per = defaultdict(list)
    for f in sorted((root / "strata").glob("*.json")):
        s = json.loads(f.read_text())
        lab = f"{s['cell_type']} | {s['tissue']} | {s['assay']} | {s['dataset_id'][:8]} ({s['n_cells']} cells, {s['n_donors']} donors)"
        for p in s["pairs"]:
            if "hit" not in p:
                continue
            per[(p["tf1"], p["tf2"])].append((p["hit"], p["r"], p["q"], lab))
    return man, per


def rna_status(per, pairs):
    if per is None:
        return {"status": "unassessed"}
    results = [per.get(p, []) for p in pairs]
    if any(not r for r in results):
        return {"status": "unassessed"}   # not detected / not in any >=100-cell stratum
    hit_each = [[x for x in r if x[0]] for r in results]
    best = [max(r, key=lambda x: x[1]) for r in (h or r0 for h, r0 in zip(hit_each, results))]
    out = {"best_r": round(min(b[1] for b in best), 4), "min_q": max(min(x[2] for x in (h or r0)) for h, r0 in zip(hit_each, results)),
           "n_hit_strata": min(len(h) for h in hit_each), "n_tested_strata": min(len(r) for r in results),
           "best_stratum": best[0][3]}
    out["min_q"] = float(f"{out['min_q']:.3g}")
    out["status"] = "hit" if all(hit_each) else "screened_no_hit"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--line", action="append", required=True,
                    help="ATAC_DATASET:CONTEXT:pairs_cbp.tsv:RNA1_DIR:RNA2_DIR")
    ap.add_argument("--register", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sym = {}
    for row in csv.DictReader(open(a.register), delimiter="\t"):
        g = row.get("Ensembl_gene_ID", "").split(".")[0]
        if g.startswith("ENSG"):
            sym.setdefault(row["Symbol"].upper(), g)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rows, provenance = [], []
    for spec in a.line:
        ds, ctx, cbp, r1dir, r2dir = spec.split(":")
        m1, p1 = load_screen(Path(r1dir)); m2, p2 = load_screen(Path(r2dir))
        provenance.append({"atac_dataset": ds, "atac_context": ctx, "pairs_cbp": cbp, "rna1": m1 | {"datasets": len(m1["datasets"])},
                           "rna2": m2 | {"datasets": len(m2["datasets"])}})
        for r in csv.DictReader(open(cbp), delimiter="\t"):
            names = r["tfs"].split("+")
            ids = [sym.get(n.upper()) if "::" not in n else None for n in names]
            sub = [tuple(sorted((ids[i], ids[j]))) for i in range(len(ids)) for j in range(i + 1, len(ids))] if all(ids) else []
            a_stat = ("supported" if r["contribution_supported"] == "yes" else "tested_no_support") if r["status"] == "tested" else "unassessed"
            e1 = rna_status(p1, sub) if sub else {"status": "unassessed"}
            e2 = rna_status(p2, sub) if sub else {"status": "unassessed"}
            s1, s2 = e1["status"], e2["status"]
            views = []
            if a_stat == "supported":
                if s1 == "hit": views.append("atac_rna1")
                if s2 == "hit": views.append("atac_rna2")
                if s1 == s2 == "screened_no_hit": views.append("atac_only")
            elif a_stat == "tested_no_support" and s1 == s2 == "hit":
                views.append("rna_only")
            note = ("composite/unknown TF label: RNA unassessed" if not sub else
                    "both RNA screens hit" if s1 == s2 == "hit" else
                    "an RNA screen is unassessed (gene undetected or no >=100-cell stratum)" if "unassessed" in (s1, s2) else "")
            pct = [r.get(f"median_site_pct_{i}") for i in (1, 2, 3) if r.get(f"median_site_pct_{i}")]
            row = {"atac_dataset": ds, "atac_context": ctx, "rna_context": m1["context"], "mapping_level": "proxy",
                   "mapping_note": m1["mapping"], "k": r["k"], "tfs": r["tfs"], "motifs": r["motifs"],
                   "tf1_ensembl": ids[0] or "", "tf1_symbol": names[0], "tf2_ensembl": ids[1] or "", "tf2_symbol": names[1],
                   "atac_status": a_stat, "cbp_median_site_pct": ";".join(pct), "cbp_synergy": r["synergy"],
                   "cbp_synergy_q": r["synergy_q"], "cbp_cooperative": r["cooperative"],
                   "evidence_view": "|".join(views) or "unassessed", "note": note,
                   "cell_type_id": m1["context"], "tissue_id": "", "disease_id": "PATO:0000461"}
            for tag, e in (("rna1", e1), ("rna2", e2)):
                for k in ("status", "best_r", "min_q", "n_hit_strata", "n_tested_strata", "best_stratum"):
                    row[f"{tag}_{k}"] = e.get(k, "")
            rows.append(row)

    def w(name, rr, f=FIELDS):
        with open(out / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, f, delimiter="\t", extrasaction="ignore"); wr.writeheader(); wr.writerows(rr)
    w("master.tsv", rows)
    counts = {}
    for v, fn in VIEWS.items():
        sel = [r for r in rows if v in r["evidence_view"].split("|")]
        w(fn, sel); counts[v] = {d: sum(r["atac_dataset"] == d for r in sel) for d in sorted({r["atac_dataset"] for r in rows})}
    w("05_unassessed.tsv", [r for r in rows if r["evidence_view"] == "unassessed"])
    pair_rows = [r for r in rows if r["k"] == "2" and r["tf1_ensembl"] and r["tf2_ensembl"]]
    mods = module_rows(pair_rows)
    w("modules.tsv", mods, ["module_id", "view", "atac_dataset", "atac_context", "cell_type_id", "tissue_id", "disease_id",
                            "n_tfs", "n_pairs", "tf_ensembl", "tf_symbols"])
    summ = {"counts": counts, "master_rows": len(rows), "n_modules": len(mods), "lines": provenance}
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(counts, indent=1))


if __name__ == "__main__":
    main()
