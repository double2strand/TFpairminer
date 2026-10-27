"""TFpairminer v2 - ChromBPNet evidence layer.

Stages (each runnable alone, each shardable where it matters):
  train    fragments -> MACS3 peaks -> GC-matched non-peaks -> bias model -> ChromBPNet  (official `chrombpnet` CLI, GPU)
           or --pretrained <url|path>: fetch an ENCODE/Kundaje-lab ChromBPNet model instead of training.
  contrib  contribution scores of the bias-corrected *counts* head over accessible windows
           (--method ig: integrated gradients from dinucleotide-shuffled references, TF2/CPU-friendly;
            --method deepshap: official `chrombpnet contribs_bw` DeepSHAP, needs the chrombpnet package/GPU).
  modisco  TF-MoDISco (modisco-lite) on the hypothetical contributions; patterns mapped to JASPAR by PWM similarity.
  pairs    per candidate motif pair/set: contribution support of each motif site vs. matched null positions,
           and in-silico marginalisation (shuffle site A, B, A+B[, C ...]) -> cooperative (non-additive) effect.

Outputs are evidence about what a *sequence model trained on accessibility* has learned; they are not TF occupancy
and not proof of physical interaction.
"""
import argparse, json, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np

from . import engine as E

INLEN, OUTLEN = 2114, 1000
CENTRE0 = (INLEN - OUTLEN) // 2           # first input bp of the central 1000 bp


# ------------------------------------------------------------------ helpers
def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def onehot(codes):
    x = np.zeros((len(codes), 4), np.float32); ok = codes < 4
    x[np.arange(len(codes))[ok], codes[ok]] = 1
    return x


def dinuc_shuffle(x, rng):
    """Dinucleotide-preserving shuffle of a one-hot sequence (Altschul-Erickson via random Eulerian walk)."""
    seq = x.argmax(1); seq[x.sum(1) == 0] = 4
    n = len(seq)
    if n < 3:
        return x.copy()
    nxt = {b: [] for b in range(5)}
    for i in range(n - 1):
        nxt[seq[i]].append(seq[i + 1])
    last = seq[-1]
    # Altschul-Erickson: random spanning arborescence into `last` (random walk / Wilson), its edges used last
    present = {b for b in range(5) if nxt[b]} | {last}
    in_tree = {last}; succ = {}
    for b in present:
        u = b
        while u not in in_tree:
            succ[u] = nxt[u][rng.integers(len(nxt[u]))]; u = succ[u]
        u = b
        while u not in in_tree:
            in_tree.add(u); u = succ[u]
    for b in present:
        if not nxt[b]:
            continue
        e = list(nxt[b])
        if b != last:
            e.remove(succ[b]); rng.shuffle(e); e.append(succ[b])
        else:
            rng.shuffle(e)
        nxt[b] = e
    out = [seq[0]]; ptr = {b: 0 for b in nxt}
    for _ in range(n - 1):
        c = out[-1]; out.append(nxt[c][ptr[c]]); ptr[c] += 1
    return onehot(np.array(out[:n]))


def find_model(path):
    p = Path(path)
    if p.is_file():
        return p
    for cand in sorted(p.rglob("*.h5")):
        if "nobias" in cand.name:
            return cand
    raise FileNotFoundError(f"no *nobias*.h5 ChromBPNet model under {path}")


def load_model(path):
    import tensorflow as tf
    km = tf.keras.models.load_model(str(find_model(path)), compile=False)

    @tf.function(input_signature=[tf.TensorSpec([None, INLEN, 4], tf.float32)])
    def grad(x):
        with tf.GradientTape() as g:
            g.watch(x); y = tf.reduce_sum(km(x, training=False)[1])
        return g.gradient(y, x)

    @tf.function(input_signature=[tf.TensorSpec([None, INLEN, 4], tf.float32)])
    def counts(x):
        return km(x, training=False)[1][:, 0]

    return km, grad, counts


def batched(f, X, bs=32):
    import tensorflow as tf
    return np.concatenate([f(tf.constant(X[i:i + bs])).numpy() for i in range(0, len(X), bs)], 0) if len(X) else np.zeros(0)


def window_inputs(genome, rows):
    out = np.zeros((len(rows), INLEN, 4), np.float32)
    for i, r in enumerate(rows):
        c = (int(r["start"]) + int(r["end"])) // 2; s = c - INLEN // 2
        codes = np.asarray(genome.get(r["chrom"])[max(0, s):s + INLEN]).astype(np.int64)
        if len(codes) < INLEN:
            codes = np.r_[codes, np.full(INLEN - len(codes), 4)]
        out[i] = onehot(codes)
    return out


def win_offset(r):
    """input-coordinate of window bp 0"""
    c = (int(r["start"]) + int(r["end"])) // 2
    return INLEN // 2 - (c - int(r["start"]))


def shard_slice(n, shard, nshards):
    per = math.ceil(n / nshards); return slice(shard * per, min(n, (shard + 1) * per))


# ------------------------------------------------------------------ candidates
WATCH = [("FOXA2", "HNF1A"), ("FOXA1", "HNF4A"), ("CEBPB", "HNF1A"), ("GATA1", "TAL1::TCF3"), ("FOS::JUN", "TEAD1"),
         ("HNF4A", "HNF1A"), ("FOXA2", "HNF4A"), ("GATA1", "KLF1"), ("CEBPA", "FOXA2"), ("NFE2", "GATA1")]


def plan_candidates(sample_dir, site_index=None, n_sample=60, n_site=60, max_pairs=150, max_sets=30):
    rows = E.read_tsv(Path(sample_dir) / "pairs.tsv")
    byname = {frozenset((r["tf_a"], r["tf_b"])): r for r in rows}
    chosen, why = [], {}

    def add(r, reason):
        k = frozenset((r["tf_a"], r["tf_b"]))
        if k not in why:
            chosen.append(r); why[k] = reason
    good = [r for r in rows if r["status"] == "priority_distinct" and r["replicate_support"] == "yes" and int(r["n_both"]) >= 30]
    for r in sorted(good, key=lambda r: -float(r["lift"]))[:n_sample]:
        add(r, "sample_top")
    if site_index and Path(site_index).exists():
        d = json.loads(Path(site_index).read_text()); F = d["pair_fields"]; k = 0
        for p in d["pairs"]:
            p = dict(zip(F, p))
            if "broad" in (p.get("flags") or "") or p["rating"] not in ("strong", "moderate"):
                continue
            r = byname.get(frozenset((p["tf_a"], p["tf_b"])))
            if r and r["q"] and float(r["q"]) <= 0.05:
                add(r, "cross_profile"); k += 1
            if k >= n_site:
                break
        for p in d["pairs"][:20]:                       # a few broad pairs as positive-control-for-composition
            p = dict(zip(F, p)); r = byname.get(frozenset((p["tf_a"], p["tf_b"])))
            if r and "broad" in (p.get("flags") or ""):
                add(r, "broad_control")
    for a, b in WATCH:
        r = byname.get(frozenset((a, b)))
        if r:
            add(r, "watchlist")
    chosen = chosen[:max_pairs]
    sets = [s for s in E.read_tsv(Path(sample_dir) / "motif_sets.tsv") if s["status"] == "synergistic_set"]
    sets = sorted(sets, key=lambda s: float(s["q_synergy"]))[:max_sets]
    return [{"tfs": [r["tf_a"], r["tf_b"]], "motifs": [r["motif_a"], r["motif_b"]], "reason": why[frozenset((r["tf_a"], r["tf_b"]))],
             "motif_lift": r["lift"], "motif_q": r["q"], "motif_status": r["status"]} for r in chosen] + \
           [{"tfs": s["tfs"].split("+"), "motifs": s["motifs"].split("+"), "reason": "synergistic_set",
             "motif_lift": s["synergy_lift"], "motif_q": s["q_synergy"], "motif_status": s["status"]} for s in sets]


def best_sites(hits, motif_idx, mlen):
    """dict window -> {motif: (pos, len, score)} keeping the best-scoring hit per motif/window"""
    want = set(motif_idx)
    sel = np.isin(hits["motif"], list(want))
    W, M, P, S = hits["window"][sel], hits["motif"][sel], hits["pos"][sel], hits["score"][sel]
    out = {}
    for w, m, p, s in zip(W, M, P, S):
        d = out.setdefault(int(w), {})
        if m not in d or s > d[m][2]:
            d[int(m)] = (int(p), int(mlen[m]), float(s))
    return out


def joint_windows(sites, mids, need_nonoverlap=True):
    res = []
    for w, d in sites.items():
        if all(m in d for m in mids):
            iv = sorted((d[m][0], d[m][0] + d[m][1]) for m in mids)
            if need_nonoverlap and any(iv[i + 1][0] < iv[i][1] for i in range(len(iv) - 1)):
                continue
            res.append(w)
    return res


# ------------------------------------------------------------------ stages
def stage_train(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.pretrained:
        src = a.pretrained; dst = out / "model" / "chrombpnet_nobias.h5"; dst.parent.mkdir(parents=True, exist_ok=True)
        if src.startswith("http"):
            subprocess.run(["curl", "-fsSL", "-o", str(dst), src], check=True)
        else:
            import shutil; shutil.copy(find_model(src), dst)
        (out / "model" / "SOURCE.json").write_text(json.dumps({"pretrained": src, "sha256": E.sha256(dst), "note":
            "pretrained ChromBPNet (bias-corrected) - trained on bulk ATAC of the matching cell line, not on these fragments"}))
        log("pretrained model ->", dst); return
    frags = " ".join(a.fragments)
    bias_block = "BIAS=" + a.bias_model if a.bias_model else f'''chrombpnet bias pipeline -ifrag $OUT/fragments.tsv.gz -d ATAC -g $GENOME -c $SIZES \\
  -p $OUT/peaks/peaks.no_bl.narrowPeak -n $OUT/peaks/bg_negatives.bed -fl $FOLD -b 0.5 -o $OUT/bias -fp {a.name}
BIAS=$OUT/bias/models/{a.name}_bias.h5'''
    sh = f"""#!/bin/bash
# Generated by tfpm v2 chrombpnet train - official ChromBPNet workflow for one profile. Needs a GPU.
set -euo pipefail
OUT={out}; mkdir -p $OUT/{{peaks,bias,chrombpnet}}
GENOME={a.fasta}; SIZES={a.chrom_sizes}; BL={a.blacklist}; FOLD={a.fold}
# 1. merge replicate fragments (5-col hg38 fragments)
zcat {frags} | sort -k1,1 -k2,2n -S 4G --parallel=8 | (command -v bgzip >/dev/null && bgzip -@8 || gzip -c) > $OUT/fragments.tsv.gz
# 2. MACS3 peaks on Tn5 cut sites (ENCODE ATAC style), blacklist-filtered, top {a.max_peaks}
zcat $OUT/fragments.tsv.gz | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$2+1,"N",".","+"; print $1,$3-1,$3,"N",".","-"}}' > $OUT/peaks/cuts.bed
MACS=$(command -v macs3 || command -v macs2)   # kundajelab/chrombpnet image ships macs2
$MACS callpeak -t $OUT/peaks/cuts.bed -f BED -n {a.name} -g hs -p 0.01 --shift -75 --extsize 150 --nomodel \\
  --keep-dup all --call-summits --outdir $OUT/peaks
bedtools slop -i $BL -g $SIZES -b 1057 > $OUT/peaks/bl_slop.bed
bedtools intersect -v -a $OUT/peaks/{a.name}_peaks.narrowPeak -b $OUT/peaks/bl_slop.bed | sort -k8,8gr | head -n {a.max_peaks} \\
  > $OUT/peaks/peaks.no_bl.narrowPeak
rm -f $OUT/peaks/cuts.bed
# 3. GC-matched background regions
chrombpnet prep nonpeaks -g $GENOME -p $OUT/peaks/peaks.no_bl.narrowPeak -c $SIZES -fl $FOLD -br $BL -o $OUT/peaks/bg
# 4. bias model (or pass --bias-model to reuse a published ATAC bias model)
{bias_block}
# 5. ChromBPNet (writes bigwigs of observed/predicted signal + the bias-corrected *_nobias.h5)
chrombpnet pipeline -ifrag $OUT/fragments.tsv.gz -d ATAC -g $GENOME -c $SIZES -p $OUT/peaks/peaks.no_bl.narrowPeak \\
  -n $OUT/peaks/bg_negatives.bed -fl $FOLD -b $BIAS -o $OUT/chrombpnet
mkdir -p $OUT/model && cp $OUT/chrombpnet/models/chrombpnet_nobias.h5 $OUT/model/chrombpnet_nobias.h5
"""
    (out / "train.sh").write_text(sh); os.chmod(out / "train.sh", 0o755)
    log("wrote", out / "train.sh")
    if a.execute:
        subprocess.run(["bash", str(out / "train.sh")], check=True)


def stage_contrib(a):
    sd = Path(a.sample_dir); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cand_f = out / "candidates.json"
    if not cand_f.exists():
        cand_f.write_text(json.dumps(plan_candidates(sd, a.site_index), indent=0))
    cands = json.loads(cand_f.read_text())
    rows = E.read_tsv(sd / "windows.tsv"); hits = dict(np.load(sd / "hits.npz"))
    names = list(hits["motif_ids"]); idx = {m: i for i, m in enumerate(names)}
    # window universe: most-accessible pooled windows + joint windows of every candidate
    pooled = [i for i, r in enumerate(rows) if "pooled" in r["sets"]]
    top = sorted(pooled, key=lambda i: -int(rows[i]["n_cuts"]))[:a.n_top]
    allm = sorted({idx[m] for c in cands for m in c["motifs"] if m in idx})
    sites = best_sites(hits, allm, hits["motif_len"]); rng = np.random.default_rng(0); chosen = set(top)
    for c in cands:
        mids = [idx[m] for m in c["motifs"] if m in idx]
        jw = [w for w in joint_windows(sites, mids) if "pooled" in rows[w]["sets"]]
        rng.shuffle(jw); chosen.update(jw[:a.per_candidate])
    wins = sorted(chosen)
    (out / "contrib_windows.json").write_text(json.dumps(wins))
    sl = shard_slice(len(wins), a.shard, a.nshards); my = wins[sl]
    log(f"contrib: {len(wins)} windows total, shard {a.shard}/{a.nshards}: {len(my)}; method={a.method}")
    if a.method == "deepshap":
        bed = out / f"contrib.shard{a.shard}.bed"
        with open(bed, "w") as f:
            for w in my:
                r = rows[w]; c = (int(r["start"]) + int(r["end"])) // 2
                f.write(f"{r['chrom']}\t{c - 500}\t{c + 500}\t{r['window_id']}\t0\t.\t0\t0\t0\t500\n")
        subprocess.run(["chrombpnet", "contribs_bw", "-m", str(find_model(a.model)), "-r", str(bed), "-g", a.fasta,
                        "-c", a.chrom_sizes, "-op", str(out / f"deepshap.shard{a.shard}"), "-pc", "counts"], check=True)
        log("DeepSHAP h5 written; run `contrib --method import-deepshap` equivalents on the cluster (see README)"); return
    km, grad, counts = load_model(a.model); g = E.Genome(a.genome)
    X = window_inputs(g, [rows[w] for w in my])
    t0 = time.time(); attr = np.zeros_like(X); hyp = np.zeros_like(X)
    for r_i in range(a.n_refs):
        R = np.stack([dinuc_shuffle(x, np.random.default_rng(1000 * r_i + i)) for i, x in enumerate(X)])
        G = np.zeros_like(X)
        for s in range(a.steps):
            al = (s + 0.5) / a.steps
            G += batched(grad, R + al * (X - R), a.batch)
        G /= a.steps; hyp += G / a.n_refs; attr += (G * (X - R)) / a.n_refs
        log(f"  ref {r_i + 1}/{a.n_refs} done ({time.time() - t0:.0f}s)")
    pred = batched(counts, X, a.batch)
    c0, c1 = CENTRE0, CENTRE0 + OUTLEN
    proj = attr.sum(-1)[:, c0:c1]                                    # projected contributions (bp) on centre 1 kb
    hypc = (hyp - hyp.mean(-1, keepdims=True))[:, c0:c1]             # mean-normalised hypothetical contributions
    np.savez_compressed(out / f"contrib.shard{a.shard}.npz", windows=np.array(my, np.int32),
                        onehot=X[:, c0:c1].astype(np.int8), proj=proj.astype(np.float16), hyp=hypc.astype(np.float16),
                        pred_logcount=pred.astype(np.float32), method=f"IG(dinuc-shuffled refs={a.n_refs}, steps={a.steps}) counts head")
    log(f"contrib done: {len(my)} windows in {time.time() - t0:.0f}s")


def load_contrib(out):
    parts = sorted(Path(out).glob("contrib.shard*.npz"))
    if not parts:
        raise FileNotFoundError(f"no contrib shards in {out}")
    zs = [np.load(p) for p in parts]
    return {k: np.concatenate([z[k] for z in zs]) for k in ("windows", "onehot", "proj", "hyp", "pred_logcount")} | \
           {"method": str(zs[0]["method"])}


def stage_modisco(a):
    out = Path(a.out); C = load_contrib(out)
    import modiscolite
    oh = C["onehot"].astype(np.float32); hyp = C["hyp"].astype(np.float32)
    log(f"modisco on {len(oh)} x {oh.shape[1]} bp")
    pos, neg = modiscolite.tfmodisco.TFMoDISco(hypothetical_contribs=hyp, one_hot=oh, max_seqlets_per_metacluster=a.max_seqlets,
                                               sliding_window_size=20, flank_size=5, target_seqlet_fdr=0.05, n_leiden_runs=2)
    modiscolite.io.save_hdf5(str(out / "modisco.h5"), pos, neg, 5)
    motifs = E.read_meme(a.meme); rows = []
    for sign, pats in (("pos", pos), ("neg", neg)):
        for k, p in enumerate(pats or []):
            ppm = np.asarray(p.sequence); cwm = np.asarray(p.contrib_scores)
            imp = np.abs(cwm).sum(1); keep = np.where(imp >= 0.3 * imp.max())[0]
            ppm_t = ppm[keep.min():keep.max() + 1] if len(keep) else ppm
            sims = sorted(((E.pwm_similarity(ppm_t, m["pwm"]), m["name"], m["id"]) for m in motifs), reverse=True)[:3]
            rows.append({"pattern": f"{sign}_{k}", "n_seqlets": len(p.seqlets), "width": len(ppm_t),
                         "consensus": "".join("ACGT"[i] for i in ppm_t.argmax(1)),
                         "best_jaspar": sims[0][1], "best_jaspar_id": sims[0][2], "similarity": round(sims[0][0], 3),
                         "matches_top3": ";".join(f"{n}:{s:.2f}" for s, n, _ in sims)})
    E.write_tsv(out / "modisco_patterns.tsv", rows, list(rows[0].keys()) if rows else ["pattern"])
    log(f"modisco: {len(rows)} patterns")


def _site_contrib(proj, w_off, p, L):
    s = w_off + p - CENTRE0
    return float(proj[s:s + L].sum()) if 0 <= s and s + L <= proj.shape[0] else np.nan


def stage_pairs(a):
    sd = Path(a.sample_dir); out = Path(a.out)
    rows = E.read_tsv(sd / "windows.tsv"); hits = dict(np.load(sd / "hits.npz"))
    names = list(hits["motif_ids"]); idx = {m: i for i, m in enumerate(names)}; mlen = hits["motif_len"]
    cands = json.loads((out / "candidates.json").read_text())
    sl = shard_slice(len(cands), a.shard, a.nshards); my = list(enumerate(cands))[sl]
    C = load_contrib(out); crow = {int(w): i for i, w in enumerate(C["windows"])}
    proj = C["proj"].astype(np.float32)
    allm = sorted({idx[m] for c in cands for m in c["motifs"] if m in idx})
    sites = best_sites(hits, allm, mlen)
    km, grad, counts = load_model(a.model); g = E.Genome(a.genome); rng = np.random.default_rng(7)
    # null for site contribution: random positions inside the 300-bp window, per length
    nullcache = {}

    def null_dist(L):
        if L not in nullcache:
            vals = []
            for i in rng.choice(len(C["windows"]), size=min(600, len(C["windows"])), replace=False):
                r = rows[int(C["windows"][i])]; wl = int(r["end"]) - int(r["start"]); off = win_offset(r)
                for p in rng.integers(0, wl - L, 3):
                    vals.append(_site_contrib(proj[i], off, int(p), L))
            nullcache[L] = np.sort(np.array([v for v in vals if np.isfinite(v)]))
        return nullcache[L]
    res, tracks = [], []
    for ci, c in my:
        mids = [idx.get(m) for m in c["motifs"]]
        rec = {"cand": ci, "tfs": "+".join(c["tfs"]), "motifs": "+".join(c["motifs"]), "k": len(c["tfs"]), "reason": c["reason"],
               "motif_lift": c["motif_lift"], "motif_q": c["motif_q"]}
        if None in mids:
            res.append(rec | {"status": "motif_missing"}); continue
        jw = [w for w in joint_windows(sites, mids) if w in crow]
        rec["n_windows"] = len(jw)
        if len(jw) < a.min_windows:
            res.append(rec | {"status": "too_few_windows"}); continue
        # ---- contribution support
        pct = np.zeros((len(jw), len(mids)))
        for j, w in enumerate(jw):
            r = rows[w]; off = win_offset(r)
            for k, m in enumerate(mids):
                p, L, _ = sites[w][m]; nd = null_dist(L)
                pct[j, k] = np.searchsorted(nd, _site_contrib(proj[crow[w]], off, p, L)) / max(1, len(nd))
        med = np.median(pct, 0)
        both = float(np.mean((pct >= a.site_pct).all(1)))
        rec |= {f"median_site_pct_{k + 1}": round(float(v), 3) for k, v in enumerate(med)}
        rec["frac_windows_all_sites_high"] = round(both, 3)
        supported = bool((med >= a.site_pct).all())
        # ---- in-silico marginalisation (shuffle site bases; average over seeds)
        jw_s = list(jw); rng.shuffle(jw_s); jw_s = jw_s[:a.ablate_windows]
        X = window_inputs(g, [rows[w] for w in jw_s])
        subsets = [()] + [(k,) for k in range(len(mids))] + ([tuple(range(len(mids)))] if len(mids) > 1 else [])
        if len(mids) == 3:
            subsets += [(0, 1), (0, 2), (1, 2)]
        pred = {s: np.zeros(len(X)) for s in subsets}
        base = batched(counts, X, a.batch)
        for seed in range(a.seeds):
            r2 = np.random.default_rng(100 + seed)
            for s in subsets[1:]:
                Xs = X.copy()
                for j, w in enumerate(jw_s):
                    off = win_offset(rows[w])
                    for k in s:
                        p, L, _ = sites[w][mids[k]]; a0 = off + p
                        Xs[j, a0:a0 + L] = Xs[j, a0:a0 + L][r2.permutation(L)]
                pred[s] += batched(counts, Xs, a.batch) / a.seeds
        delta = {s: base - pred[s] for s in subsets[1:]}                  # loss of log-counts when sites destroyed
        full = tuple(range(len(mids)))
        syn = delta[full] - sum(delta[(k,)] for k in range(len(mids)))     # >0: joint effect larger than additive
        from scipy.stats import wilcoxon
        try:
            p_syn = float(wilcoxon(syn).pvalue) if np.any(syn != 0) else 1.0
        except ValueError:
            p_syn = 1.0
        rec |= {"n_ablated": len(X), "pred_logcount_wt": round(float(base.mean()), 4),
                **{f"delta_{k + 1}": round(float(np.mean(delta[(k,)])), 4) for k in range(len(mids))},
                "delta_joint": round(float(np.mean(delta[full])), 4), "synergy": round(float(np.mean(syn)), 4),
                "synergy_median": round(float(np.median(syn)), 4), "synergy_p": p_syn,
                "frac_windows_synergy_pos": round(float(np.mean(syn > 0)), 3),
                "contribution_supported": "yes" if supported else "no", "status": "tested"}
        if len(mids) == 3:
            rec["synergy_vs_best_pair"] = round(float(np.mean(delta[full] - max((delta[(0, 1)], delta[(0, 2)], delta[(1, 2)]), key=np.mean)
                                                              - 0)), 4)
        res.append(rec)
        # example tracks: 2 windows with the highest joint contribution
        sc = [(proj[crow[w]].sum(), w) for w in jw]
        for _, w in sorted(sc, reverse=True)[:a.n_tracks]:
            r = rows[w]; off = win_offset(r); i = crow[w]; wl = int(r["end"]) - int(r["start"])
            s0 = off - CENTRE0; seg = proj[i, s0:s0 + wl]; ohs = C["onehot"][i, s0:s0 + wl]
            tracks.append({"cand": ci, "tfs": rec["tfs"], "window": r["window_id"], "chrom": r["chrom"], "start": int(r["start"]),
                           "seq": "".join("ACGTN"[int(np.argmax(b)) if b.sum() else 4] for b in ohs),
                           "contrib": [round(float(v), 4) for v in seg],
                           "sites": [{"tf": c["tfs"][k], "pos": sites[w][m][0], "len": sites[w][m][1]} for k, m in enumerate(mids)]})
        log(f"[{ci}] {rec['tfs']}: n={len(jw)} supported={rec['contribution_supported']} synergy={rec['synergy']:+.4f} p={p_syn:.2g}")
    (out / f"pairs_cbp.shard{a.shard}.json").write_text(json.dumps({"results": res, "tracks": tracks}))
    if a.nshards == 1 or a.merge:
        merge(out)


def merge(out):
    out = Path(out); res, tracks = [], []
    for p in sorted(out.glob("pairs_cbp.shard*.json")):
        d = json.loads(p.read_text()); res += d["results"]; tracks += d["tracks"]
    tested = [r for r in res if r.get("status") == "tested"]
    q = E.bh(np.array([r["synergy_p"] for r in tested])) if tested else []
    for r, qq in zip(tested, q):
        r["synergy_q"] = float(qq)
        r["cooperative"] = ("synergistic" if qq <= 0.05 and r["synergy"] > 0 else
                            "sub-additive" if qq <= 0.05 and r["synergy"] < 0 else "additive/ns")
    fields = sorted({k for r in res for k in r}, key=lambda k: (k not in ("tfs", "motifs", "k", "status"), k))
    E.write_tsv(out / "pairs_cbp.tsv", res, fields)
    (out / "tracks.json").write_text(json.dumps(tracks))
    s = {"n_candidates": len(res), "n_tested": len(tested),
         "n_contribution_supported": sum(r["contribution_supported"] == "yes" for r in tested),
         "n_synergistic": sum(r.get("cooperative") == "synergistic" for r in tested),
         "n_subadditive": sum(r.get("cooperative") == "sub-additive" for r in tested)}
    (out / "cbp_summary.json").write_text(json.dumps(s, indent=1)); log("merged", s)


# ------------------------------------------------------------------ CLI
def main(argv=None):
    ap = argparse.ArgumentParser(prog="tfpm v2 chrombpnet")
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train"); t.add_argument("--out", required=True); t.add_argument("--name", default="profile")
    t.add_argument("--pretrained"); t.add_argument("--fragments", nargs="*", default=[]); t.add_argument("--fasta")
    t.add_argument("--chrom-sizes"); t.add_argument("--blacklist"); t.add_argument("--fold"); t.add_argument("--bias-model")
    t.add_argument("--max-peaks", type=int, default=300000); t.add_argument("--execute", action="store_true")
    for n in ("contrib", "modisco", "pairs", "merge"):
        s = sub.add_parser(n); s.add_argument("--out", required=True)
        if n in ("contrib", "pairs"):
            s.add_argument("--sample-dir", required=True); s.add_argument("--model", required=True)
            s.add_argument("--genome", required=True); s.add_argument("--shard", type=int, default=0)
            s.add_argument("--nshards", type=int, default=1); s.add_argument("--batch", type=int, default=32)
        if n == "contrib":
            s.add_argument("--site-index"); s.add_argument("--n-top", type=int, default=800)
            s.add_argument("--per-candidate", type=int, default=25); s.add_argument("--method", default="ig", choices=["ig", "deepshap"])
            s.add_argument("--n-refs", type=int, default=1); s.add_argument("--steps", type=int, default=6)
            s.add_argument("--fasta"); s.add_argument("--chrom-sizes")
        if n == "modisco":
            s.add_argument("--meme", required=True); s.add_argument("--max-seqlets", type=int, default=20000)
        if n == "pairs":
            s.add_argument("--min-windows", type=int, default=8); s.add_argument("--ablate-windows", type=int, default=40)
            s.add_argument("--seeds", type=int, default=2); s.add_argument("--site-pct", type=float, default=0.75)
            s.add_argument("--n-tracks", type=int, default=2); s.add_argument("--merge", action="store_true")
    a = ap.parse_args(argv)
    {"train": stage_train, "contrib": stage_contrib, "modisco": stage_modisco, "pairs": stage_pairs,
     "merge": lambda a: merge(a.out)}[a.cmd](a)


if __name__ == "__main__":
    main()
