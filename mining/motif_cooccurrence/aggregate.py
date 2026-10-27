"""Cross-sample aggregation -> static JSON for docs/v2 (TFpairminer v2 site)."""
import json, math, re, time
from collections import defaultdict
from pathlib import Path

import numpy as np

from .engine import read_tsv

SIG = ("priority_distinct", "family_level")
ST_CODE = {"priority_distinct": "D", "family_level": "F", "exploratory": "E"}
NOTE = ("Motif-sequence co-occurrence hypotheses from scATAC accessible windows. Not ChromBPNet "
        "contribution scores, not TF occupancy, not evidence of physical interaction.")


def _f(x, nd=4):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    return float(f"{v:.{nd}g}")


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def rating(c):
    d, f, units = c["n_distinct_rep"], c["n_sig"], c["n_parent_sig"]
    if d >= 3 and units >= 2:
        return "strong", "Strong hypothesis"
    if d >= 2 or (d >= 1 and (c["best_perm_p"] or 1) <= 0.01):
        return "moderate", "Moderate hypothesis"
    if d >= 1:
        return "preliminary", "Preliminary hypothesis"
    if c["n_family_sig"] >= 1 and c["n_distinct_any"] == 0:
        return "family", "Family-level only"
    return "weak", "Weak / single-replicate"


def summary_text(c, units):
    a, b = c["tf_a"], c["tf_b"]
    ctx = sorted({units[u]["context"] for u in c["sig_units"]})
    s = (f"{a} and {b} motifs co-occur in the same accessible ~300 bp windows more often than expected "
         f"after matching GC, accessibility and repeat content, in {c['n_sig']} of {c['n_units']} profiles analysed")
    if ctx:
        s += " (" + "; ".join(ctx[:4]) + ("; ..." if len(ctx) > 4 else "") + ")"
    s += f". Median lift {c['median_lift_sig']}x"
    if c["median_gap"] is not None:
        s += f", typical spacing ~{int(c['median_gap'])} bp between sites"
    s += ". "
    if c["n_distinct_rep"]:
        s += f"In {c['n_distinct_rep']} profile(s) the enrichment is seen in both replicate halves and the motifs are distinguishable. "
    if "broad" in c["family_flags"]:
        s += "It is enriched in most profiles, so it may reflect general sequence composition of accessible DNA rather than context-specific cooperation. "
    ff = [f for f in c["family_flags"] if f != "broad"]
    if ff:
        s += "Caution: " + ", ".join(sorted(ff)).replace("_", " ") + " — the signal may reflect one family/shared motif rather than two distinct factors. "
    s += "This is a sequence-level hypothesis to test (e.g. ChIP/CUT&RUN co-binding, motif mutagenesis), not a confirmed interaction."
    return s


def cbp_text(c, ev):
    tested = {u: e for u, e in ev.items() if e["st"] == "tested"}
    if not tested:
        return (" Evidence: motif co-occurrence only (ChromBPNet not yet run for this pair)." if not ev else
                " Evidence: motif co-occurrence only (ChromBPNet run, but too few windows with both sites to test).")
    sup = [u for u, e in tested.items() if e["sup"] == "yes"]
    syn = [u for u, e in tested.items() if e.get("coop") == "synergistic"]
    t = (f" Evidence: ChromBPNet contribution scores support both motif sites in {len(sup)} of {len(tested)} profiles tested"
         if sup else f" Evidence: motif co-occurrence only; ChromBPNet contribution scores did not highlight both sites ({len(tested)} profile(s) tested)")
    t += (f"; joint in-silico ablation is synergistic (non-additive) in {len(syn)}." if syn else "; no significant ablation synergy.")
    return t


def aggregate(out, site, meme=None, max_pairs=12000, max_examples=6):
    out, site = Path(out), Path(site)
    for d in ("tf", "samples"):
        (site / d).mkdir(parents=True, exist_ok=True)
    units = {}
    for sj in sorted(out.glob("samples/**/summary.json")):
        s = json.loads(sj.read_text()); d = sj.parent
        if not (d / "pairs.tsv").exists():
            continue
        u = s["sample"]
        units[u] = {"id": u, "dataset": s["dataset"], "context": s["context"], "parent": s.get("parent_sample", u),
                    "kind": s.get("kind", ""), "dir": d, "summary": {k: v for k, v in s.items() if k not in ("params",)},
                    "params": s.get("params", {})}
        wm = d / "windows_meta.json"
        if not wm.exists() and (d.parent.parent / "windows_meta.json").exists():
            wm = d.parent.parent / "windows_meta.json"
        units[u]["windows_meta"] = json.loads(wm.read_text()) if wm.exists() else {}
    U = sorted(units)
    pairs = defaultdict(dict); examples = defaultdict(list); sets = []; mods = []
    for u in U:
        d = units[u]["dir"]
        for r in read_tsv(d / "pairs.tsv"):
            a, b = sorted([r["tf_a"], r["tf_b"]])
            pairs[(a, b)][u] = {"l": _f(r["lift"]), "q": _f(r["q"], 3), "s": ST_CODE.get(r["status"], "E"),
                                "r": r["replicate_support"][:1], "n": int(r["n_both"]), "e": _f(r["expected_both"], 3),
                                "pp": _f(r.get("perm_p"), 3), "g": _f(r.get("median_gap_bp"), 3),
                                "amb": r.get("ambiguity", ""), "r1": _f(r.get("rep1_lift"), 3), "r2": _f(r.get("rep2_lift"), 3),
                                "fa": r["families_a"] if a == r["tf_a"] else r["families_b"],
                                "fb": r["families_b"] if a == r["tf_a"] else r["families_a"],
                                "ma": r["motif_a"] if a == r["tf_a"] else r["motif_b"],
                                "mb": r["motif_b"] if a == r["tf_a"] else r["motif_a"]}
        if (d / "pair_site_examples.tsv").exists():
            for r in read_tsv(d / "pair_site_examples.tsv"):
                a, b = sorted([r["tf_a"], r["tf_b"]])
                if len([e for e in examples[(a, b)] if e["u"] == u]) < 3:
                    sw = r["tf_a"] != a
                    examples[(a, b)].append({"u": u, "chrom": r["chrom"],
                        "a": [int(r["b_start" if sw else "a_start"]), int(r["b_end" if sw else "a_end"])],
                        "b": [int(r["a_start" if sw else "b_start"]), int(r["a_end" if sw else "b_end"])],
                        "w": r["window_id"], "gap": int(float(r["gap_bp"])), "cuts": int(r["n_cuts"]), "rc": r["rep_cuts"]})
        for r in read_tsv(d / "motif_sets.tsv") if (d / "motif_sets.tsv").exists() else []:
            if r["status"] == "synergistic_set" or float(r["q_synergy"]) <= 0.05:
                sets.append({"u": u, "tfs": r["tfs"].split("+"), "n": int(r["n_all3"]), "li": _f(r["lift_vs_independent"], 3),
                             "sl": _f(r["synergy_lift"], 3), "q": _f(r["q_synergy"], 3), "s": r["status"]})
        for r in read_tsv(d / "modules.tsv") if (d / "modules.tsv").exists() else []:
            mods.append({"u": u, "tfs": r["module"].split("+"), "n": int(r["n_windows_all"])})
    # ---- ChromBPNet evidence layer (optional; samples/<u>/chrombpnet/)
    cbp, cbp_tracks, cbp_sets, cbp_units = defaultdict(dict), defaultdict(list), defaultdict(dict), {}
    for u in U:
        cd = units[u]["dir"] / "chrombpnet"
        if not (cd / "pairs_cbp.tsv").exists():
            continue
        src = json.loads((cd / "model" / "SOURCE.json").read_text()) if (cd / "model" / "SOURCE.json").exists() else {}
        cbp_units[u] = {**(json.loads((cd / "cbp_summary.json").read_text()) if (cd / "cbp_summary.json").exists() else {}),
                        "model": src.get("source") or src.get("pretrained", ""), "method": src.get("note", ""),
                        "modisco": [{k: r[k] for k in ("pattern", "n_seqlets", "consensus", "best_jaspar", "best_jaspar_id")}
                                    for r in (read_tsv(cd / "modisco_patterns.tsv") if (cd / "modisco_patterns.tsv").exists() else [])][:25]}
        for r in read_tsv(cd / "pairs_cbp.tsv"):
            tfs = tuple(sorted(r["tfs"].split("+")))
            e = {"st": r["status"], "n": int(float(r["n_windows"])) if r.get("n_windows") else 0}
            if r["status"] == "tested":
                e |= {"sup": r["contribution_supported"], "syn": _f(r["synergy"], 4), "p": _f(r["synergy_p"], 3),
                      "q": _f(r.get("synergy_q"), 3), "coop": r.get("cooperative", ""),
                      "pct": [_f(r[k], 3) for k in sorted(r) if k.startswith("median_site_pct_") and r[k] != ""],
                      "dj": _f(r["delta_joint"], 4), "ds": [_f(r[k], 4) for k in sorted(r) if k.startswith("delta_") and k[6:].isdigit() and r[k] != ""]}
            (cbp if len(tfs) == 2 else cbp_sets)[tfs][u] = e
        if (cd / "tracks.json").exists():
            for t in json.loads((cd / "tracks.json").read_text()):
                tfs = tuple(sorted(t["tfs"].split("+")))
                if len(cbp_tracks[tfs]) < 4:
                    cbp_tracks[tfs].append({**t, "u": u})

    # ---- ATAC+RNA evidence layer (optional; out/atac_rna/ from tf_pair_evidence.targeted_integrate)
    rna = defaultdict(dict)
    ar = out / "atac_rna"
    if (ar / "master.tsv").exists():
        for r in read_tsv(ar / "master.tsv"):
            if r["k"] == "2":
                rna[tuple(sorted(r["tfs"].split("+")))][r["atac_dataset"]] = {
                    k: r[k] for k in ("evidence_view", "atac_status", "rna1_status", "rna1_best_r", "rna1_min_q",
                                      "rna2_status", "rna2_best_r", "rna2_min_q", "mapping_note", "rna_context")}
        rows = list(read_tsv(ar / "master.tsv"))
        (site / "atac_rna.json").write_text(json.dumps({
            "summary": json.loads((ar / "summary.json").read_text()), "rows": rows,
            "modules": list(read_tsv(ar / "modules.tsv"))}, separators=(",", ":"), default=str))

    def rna_label(tfs):
        views = {v for e in rna.get(tfs, {}).values() for v in e["evidence_view"].split("|")}
        for v in ("atac_rna1|atac_rna2", "atac_rna1", "atac_rna2", "rna_only", "atac_only"):
            if set(v.split("|")) <= views:
                return {"atac_rna1|atac_rna2": "atac_rna_both"}.get(v, v)
        return "unassessed" if tfs in rna else "not_run"

    def rna_text(tfs):
        ev = rna.get(tfs)
        if not ev:
            return ""
        bits = []
        for u, e in sorted(ev.items()):
            st = lambda k: {"hit": "co-expressed", "screened_no_hit": "no significant co-expression",
                            "unassessed": "unassessed"}.get(e[k], e[k])
            bits.append(f"{u.replace('CCLA_', '')} ({e['rna_context'].replace('_', ' ')}, proxy match): "
                        f"CELLxGENE {st('rna1_status')}, Tabula Sapiens {st('rna2_status')}")
        return " RNA layer (donor-aware within-donor metacell r): " + "; ".join(bits) + "."

    def cbp_label(tfs):
        ev = cbp.get(tfs) or cbp_sets.get(tfs) or {}
        tested = [e for e in ev.values() if e["st"] == "tested"]
        if not tested:
            return "not_run" if not ev else "untestable"
        return "supported" if any(e["sup"] == "yes" for e in tested) else "not_supported"
    cons = []
    for (a, b), per in pairs.items():
        sigu = [u for u, v in per.items() if v["s"] in "DF" and v["s"] != "E"]
        if not sigu:
            continue
        c = {"tf_a": a, "tf_b": b, "n_units": len(U), "sig_units": sigu, "n_sig": len(sigu),
             "n_distinct_any": sum(per[u]["s"] == "D" for u in sigu),
             "n_distinct_rep": sum(per[u]["s"] == "D" and per[u]["r"] == "y" for u in sigu),
             "n_family_sig": sum(per[u]["s"] == "F" for u in sigu),
             "n_parent_sig": len({units[u]["parent"] for u in sigu}),
             "datasets": sorted({units[u]["dataset"] for u in sigu}),
             "median_lift_sig": _f(np.median([per[u]["l"] for u in sigu]), 3),
             "best_q": min(per[u]["q"] for u in sigu),
             "best_perm_p": min([per[u]["pp"] for u in sigu if per[u]["pp"] is not None], default=None),
             "median_gap": _f(np.median([per[u]["g"] for u in sigu if per[u]["g"] is not None]), 3)
             if any(per[u]["g"] is not None for u in sigu) else None,
             "family_flags": sorted({f.split(":")[0] for u in sigu for f in per[u]["amb"].split(";") if f}),
             "fam_a": next(iter(per.values()))["fa"], "fam_b": next(iter(per.values()))["fb"]}
        if c["n_units"] >= 4 and c["n_sig"] >= 0.75 * c["n_units"]:
            c["family_flags"] = c["family_flags"] + ["broad"]
        c["rating"], c["rating_label"] = rating(c)
        c["cbp"] = cbp_label((a, b))
        c["rna"] = rna_label((a, b))
        c["cbp_coop"] = sorted({e.get("coop") for e in cbp.get((a, b), {}).values() if e.get("coop")})
        c["score"] = round(c["n_distinct_rep"] * 3 + c["n_distinct_any"] + 0.5 * c["n_sig"]
                           + min(10, -math.log10(max(c["best_q"], 1e-300))) / 10 + (c["median_lift_sig"] or 1) / 10, 3)
        cons.append(c)
    order = {"strong": 0, "moderate": 1, "preliminary": 2, "family": 3, "weak": 4}
    cons.sort(key=lambda c: (order[c["rating"]], -c["score"]))
    cons = cons[:max_pairs]
    for c in cons:
        c["summary"] = summary_text(c, units) + cbp_text(c, cbp.get((c["tf_a"], c["tf_b"]), {})) + rna_text((c["tf_a"], c["tf_b"]))
    for s_ in sets:
        e = cbp_sets.get(tuple(sorted(s_["tfs"])), {})
        s_["cbp"] = cbp_label(tuple(sorted(s_["tfs"])))
        if s_["u"] in e and e[s_["u"]]["st"] == "tested":
            s_["cbp_e"] = e[s_["u"]]
    # per-TF detail files
    by_tf = defaultdict(list)
    for c in cons:
        by_tf[c["tf_a"]].append(c); by_tf[c["tf_b"]].append(c)
    set_by_tf = defaultdict(list)
    for s_ in sets:
        for t in s_["tfs"]:
            set_by_tf[t].append(s_)
    mod_by_tf = defaultdict(list)
    for m in mods:
        for t in m["tfs"]:
            mod_by_tf[t].append(m)
    for tf, cl in by_tf.items():
        det = {"tf": tf, "pairs": [], "sets": sorted(set_by_tf.get(tf, []), key=lambda s: s["q"] or 1)[:200],
               "set_tracks": {"+".join(k): v for k, v in cbp_tracks.items() if len(k) > 2 and tf in k},
               "modules": mod_by_tf.get(tf, [])[:100]}
        for c in cl:
            key = (c["tf_a"], c["tf_b"])
            per = {u: {k: v for k, v in e.items() if k not in ("fa", "fb", "ma", "mb")} for u, e in pairs[key].items()}
            for e in per.values():
                if not e["amb"]: del e["amb"]
            first = next(iter(pairs[key].values()))
            det["pairs"].append({**{k: v for k, v in c.items() if k != "sig_units"}, "ma": first["ma"], "mb": first["mb"],
                                 "per": per, "examples": examples.get(key, [])[:max_examples],
                                 "cbp_per": cbp.get(key, {}), "cbp_tracks": cbp_tracks.get(key, []), "rna_per": rna.get(key, {})})
        (site / "tf" / f"{slug(tf)}.json").write_text(json.dumps(det, separators=(",", ":")))
    # per-unit sample files
    for u in U:
        top = sorted([(k, v[u]) for k, v in pairs.items() if u in v and v[u]["s"] in "DF"],
                     key=lambda kv: (kv[1]["s"] != "D", kv[1]["r"] != "y", kv[1]["q"]))[:600]
        su = {"unit": {k: v for k, v in units[u].items() if k not in ("dir",)},
              "pairs": [{"tf_a": k[0], "tf_b": k[1], **v, "cbp": cbp_label(k), **({"cbp_e": cbp[k][u]} if u in cbp.get(k, {}) else {})} for k, v in top],
              "cbp": cbp_units.get(u),
              "sets": sorted([s_ for s_ in sets if s_["u"] == u], key=lambda s: s["q"] or 1)[:300],
              "modules": [m for m in mods if m["u"] == u][:100]}
        (site / "samples" / f"{slug(u)}.json").write_text(json.dumps(su, separators=(",", ":"), default=str))
    index = {"generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()), "note": NOTE,
             "units": [{"id": u, "slug": slug(u), "dataset": units[u]["dataset"], "context": units[u]["context"],
                        "parent": units[u]["parent"], "kind": units[u]["kind"],
                        **{k: units[u]["summary"].get(k) for k in ("n_windows", "n_motifs_tested", "n_pairs_tested",
                           "n_priority_distinct", "n_family_level", "n_replicate_supported_sig",
                           "n_sets_tested", "n_synergistic_sets", "n_modules")}} for u in U],
             "counts": {"pairs_listed": len(cons), "tfs": len(by_tf),
                        "by_rating": {r: sum(c["rating"] == r for c in cons) for r in order}},
             "pairs": [[c["tf_a"], c["tf_b"], c["rating"], c["n_sig"], c["n_distinct_rep"], c["median_lift_sig"],
                        c["best_q"], ",".join(c["datasets"]), "|".join(c["family_flags"]), c["median_gap"], c["cbp"],
                        ",".join(c["cbp_coop"]), c["rna"]] for c in cons],
             "pair_fields": ["tf_a", "tf_b", "rating", "n_sig", "n_distinct_rep", "median_lift", "best_q", "datasets",
                             "flags", "median_gap", "cbp", "cbp_coop", "rna"],
             "rna_counts": {k: sum(c["rna"] == k for c in cons) for k in ("atac_rna_both", "atac_rna1", "atac_rna2", "rna_only", "atac_only", "unassessed", "not_run")},
             "cbp_units": cbp_units,
             "cbp_counts": {k: sum(c["cbp"] == k for c in cons) for k in ("supported", "not_supported", "untestable", "not_run")},
             "top_sets": sorted([s_ for s_ in sets if s_["s"] == "synergistic_set"], key=lambda s: s["q"] or 1)[:300],
             "tf_slugs": {tf: slug(tf) for tf in by_tf}}
    (site / "index.json").write_text(json.dumps(index, separators=(",", ":"), default=str))
    print(json.dumps({"units": len(U), "pairs_listed": len(cons), "tfs": len(by_tf), "sets": len(sets)}))
    return index
