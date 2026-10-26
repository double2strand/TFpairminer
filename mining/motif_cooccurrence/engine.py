"""TFpairminer v2 — scATAC fragment/cCRE motif co-occurrence engine.

Generalises the HepG2 chr22 pilot (chrombpnet_pair_discovery/pilot.py) to
genome-wide, multi-sample, full-JASPAR runs with pairs AND higher-order motif
sets.  Everything produced here is a MOTIF-SEQUENCE CO-OCCURRENCE HYPOTHESIS:
not ChromBPNet attribution, not occupancy, not a physical TF interaction.

Stages (each idempotent, one sample per job so wr/local can parallelise):
  windows-frag   replicate-aware Tn5 insertion hotspots from *.fragments.tsv.gz
  windows-corpus per-cell-type pseudobulk top cCREs from a Corpus cell x cCRE h5ad
  scan           numba PWM scan of all windows (full JASPAR or a panel)
  pairs          GC x insertion-stratified co-occurrence test (+ permutations
                 for the top pairs, replicate-split support, spacing, ambiguity)
  sets           triplets beyond their best pair + modules (cliques) of pairs
"""
import gzip, hashlib, json, math, os, re, sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

BIN_BP = 50
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i; LUT[_b + 32] = _i  # soft-masked lower case too

# Curated motif classes whose PWMs cannot separate paralogous proteins.
FAMILY_PREFIXES = [
    ("AP-1 (bZIP FOS/JUN/ATF)", r"^(FOS|FOSL\d|JUN|JUNB|JUND|ATF\d|BATF\d?|BACH\d|NFE2.*|MAF.*)$"),
    ("CEBP (bZIP)", r"^(CEBP[A-G]|DBP|HLF|TEF|NFIL3)$"),
    ("E-box bHLH", r"^(MYC|MYCN|MAX|MNT|MLX.*|USF\d|TFE.*|MITF|ARNT.*|HEY\d|HES\d|BHLHE\d+|CLOCK|NPAS\d|TCF[34]|TCF12|ASCL\d|NEUROD\d|ATOH\d|TWIST\d|HAND\d|MYOD1|MYOG|MYF\d|SREBF\d)$"),
    ("GC-box SP/KLF", r"^(SP\d|KLF\d+|EGR\d|ZNF148|ZNF281|PATZ1|WT1|MAZ|VEZF1)$"),
    ("ETS", r"^(ETS\d|ETV\d|ELK\d|ELF\d|ERG|FLI1|GABPA|SPI[1BC]|EHF|ELF5|FEV|ERF)$"),
    ("Forkhead", r"^(FOX[A-S]\d*)$"),
    ("SOX/HMG", r"^(SOX\d+|SRY|LEF1|TCF7.*|HMG.*)$"),
    ("Nuclear receptor", r"^(NR\d[A-Z]\d|HNF4[AG]|RXR[ABG]|RAR[ABG]|ESRR[ABG]|PPAR[ADG]|THR[AB]|VDR|ESR\d|AR|PGR|NR3C\d|RORA|RORB|RORC|RARA)$"),
    ("GATA", r"^(GATA\d|TRPS1)$"),
    ("RUNX", r"^RUNX\d$"), ("TEAD", r"^TEAD\d$"), ("STAT", r"^STAT\d.*$"),
    ("IRF", r"^IRF\d$"), ("NFKB/REL", r"^(NFKB\d|REL[AB]?)$"), ("NFI", r"^NFI[ABCX]$"),
    ("Homeobox", r"^(HOX[A-D]\d+|CDX\d|PDX1|NKX\d-\d|LHX\d|ISL\d|PAX\d|OTX\d|CRX|PITX\d|DLX\d|MSX\d|MEIS\d|PBX\d|ONECUT\d|HNF1[AB]|POU\d?F\d|VSX\d|EMX\d|GSC|ALX\d|PHOX2[AB]|SIX\d|BARX\d|VAX\d|NANOG)$"),
    ("CTCF-like", r"^(CTCF|CTCFL)$"), ("RFX", r"^RFX\d$"), ("E2F", r"^E2F\d$"),
    ("SMAD", r"^SMAD\d$"), ("ZBTB", r"^ZBTB\d+.*$"), ("CREB/ATF", r"^(CREB\d|CREM|CREB3.*|XBP1)$"),
    ("NRF/NFY", r"^(NRF1|NFY[ABC])$"), ("TP53 family", r"^(TP53|TP63|TP73)$"),
    ("MEF2", r"^MEF2[A-D]$"), ("HSF", r"^HSF\d$"), ("GRHL/TFCP2", r"^(GRHL\d|TFCP2.*|UBP1)$"),
]
_FAM_RE = [(n, re.compile(p)) for n, p in FAMILY_PREFIXES]


def tf_family(name):
    parts = name.upper().replace("(VAR.2)", "").replace("(VAR.3)", "").split("::")
    fams = {f for p in parts for f, r in _FAM_RE if r.match(p.strip())}
    return sorted(fams)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path, rows, fields):
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(f) is None else (f"{r[f]:.6g}" if isinstance(r.get(f), float) else str(r[f])) for f in fields) + "\n")
    os.replace(tmp, path)


def read_tsv(path):
    with open_text(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


# ---------------------------------------------------------------- genome
def prep_genome(fasta, outdir, chroms=MAIN_CHROMS):
    """FASTA(.gz) -> one uint8 code array per chromosome (A=0,C=1,G=2,T=3,N=4)."""
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    want = set(chroms); name, buf = None, []

    def flush():
        if name in want and buf:
            raw = np.frombuffer(b"".join(buf), dtype=np.uint8)
            np.save(outdir / f"{name}.npy", LUT[raw])
            np.save(outdir / f"{name}.rep.npy", np.packbits(raw >= 97))  # soft-masked = repeat
    with gzip.open(fasta, "rb") if str(fasta).endswith(".gz") else open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                flush(); name = line[1:].split()[0].decode(); buf = []
            elif name in want:
                buf.append(line.rstrip())
        flush()
    (outdir / "SOURCE.json").write_text(json.dumps({"fasta": str(fasta), "chroms": sorted(want)}))


class Genome:
    def __init__(self, d):
        self.d = Path(d); self.cache = {}

    def chroms(self):
        return [c for c in MAIN_CHROMS if (self.d / f"{c}.npy").exists()]

    def repeat_frac(self, chrom, s, e):
        f = self.d / f"{chrom}.rep.npy"
        if not f.exists():
            return float("nan")
        k = ("rep", chrom)
        if k not in self.cache:
            self.cache[k] = np.load(f, mmap_mode="r")
        b = np.unpackbits(np.asarray(self.cache[k][s // 8:(e + 7) // 8]))[s % 8:s % 8 + (e - s)]
        return float(b.mean()) if len(b) else float("nan")

    def get(self, chrom):
        if chrom not in self.cache:
            self.cache[chrom] = np.load(self.d / f"{chrom}.npy", mmap_mode="r")
        return self.cache[chrom]


def read_blacklist(path):
    out = defaultdict(list)
    if path:
        with open_text(path) as fh:
            for l in fh:
                f = l.split("\t")
                if len(f) >= 3 and not l.startswith("#"):
                    out[f[0]].append((int(f[1]), int(f[2])))
    return out


# ---------------------------------------------------------------- windows from fragments
def insertion_bins(path, genome, chroms):
    """Stream a 5-column fragments file once -> {chrom: int32 bins of Tn5 insertions}."""
    import pandas as pd
    sizes = {c: len(genome.get(c)) for c in chroms}
    counts = {c: np.zeros(math.ceil(sizes[c] / BIN_BP), dtype=np.int32) for c in chroms}
    n_frag, barcodes = 0, set()
    for chunk in pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3, 4],
                             names=["c", "s", "e", "b", "n"], chunksize=2_000_000,
                             dtype={"c": str, "s": np.int64, "e": np.int64, "b": str, "n": np.int64}):
        barcodes.update(chunk["b"].unique().tolist())
        n_frag += int(chunk["n"].sum())
        for c, g in chunk.groupby("c"):
            if c not in counts:
                continue
            arr = counts[c]; lim = sizes[c]
            for pos in (g["s"].to_numpy() + 4, g["e"].to_numpy() - 5):
                ok = (pos >= 0) & (pos < lim)
                np.add.at(arr, pos[ok] // BIN_BP, g["n"].to_numpy()[ok].astype(np.int32))
    return counts, {"file": str(Path(path).resolve()), "n_fragments": n_frag,
                    "n_barcodes": len(barcodes), "sha256": sha256(path)}


def call_hotspots(rep_counts, genome, chroms, blacklist, window_bp, min_cuts, min_rep_cuts,
                  which="pooled"):
    """Pilot hotspot logic (smoothed pooled insertions, find_peaks), genome-wide.
    which = 'pooled' | 'rep<k>': candidate summits from that signal; pooled calls
    also require >= min_rep_cuts in every replicate (replicate-supported)."""
    from scipy.ndimage import uniform_filter1d
    from scipy.signal import find_peaks
    wb = math.ceil(window_bp / BIN_BP)
    out = []
    for c in chroms:
        arrs = [rc[c] for rc in rep_counts]
        sig = np.sum(arrs, axis=0) if which == "pooled" else arrs[int(which[3:]) - 1]
        smooth = uniform_filter1d(sig.astype(float), size=wb, mode="constant") * wb
        summits, _ = find_peaks(smooth, distance=wb + 2, height=min_cuts, prominence=max(2, min_cuts / 3))
        if not len(summits):
            continue
        prefixes = [np.r_[0, np.cumsum(a, dtype=np.int64)] for a in arrs]
        seq = genome.get(c); bl = sorted(blacklist.get(c, []))
        bls = np.array([b[0] for b in bl]); ble = np.array([b[1] for b in bl])
        for k in summits:
            center = int(k * BIN_BP + BIN_BP // 2); s = center - window_bp // 2; e = s + window_bp
            if s < 0 or e > len(seq):
                continue
            if len(bl):
                j = np.searchsorted(bls, e)
                if j and (ble[:j] > s).any():
                    continue
            sub = np.asarray(seq[s:e])
            if (sub == 4).sum() > window_bp * 0.05:
                continue
            lo, hi = s // BIN_BP, math.ceil(e / BIN_BP)
            per = [int(p[hi] - p[lo]) for p in prefixes]
            if which == "pooled" and min(per) < min_rep_cuts:
                continue
            gc = float(((sub == 1) | (sub == 2)).sum() / max(1, (sub < 4).sum()))
            out.append({"chrom": c, "start": s, "end": e, "score": float(smooth[k]),
                        "n_cuts": sum(per), "rep_cuts": per, "gc": round(gc, 4)})
    return out


def windows_from_fragments(fragments, genome_dir, blacklist, out, window_bp=300, max_windows=20000,
                           min_cuts=6, min_rep_cuts=2, chroms=None, label=None):
    g = Genome(genome_dir); chroms = chroms or g.chroms(); bl = read_blacklist(blacklist)
    reps, info = [], []
    for f in fragments:
        c, i = insertion_bins(f, g, chroms); reps.append(c); info.append(i)
    sets = {"pooled": call_hotspots(reps, g, chroms, bl, window_bp, min_cuts, min_rep_cuts, "pooled")}
    if len(reps) >= 2:
        for k in range(len(reps)):
            sets[f"rep{k+1}"] = call_hotspots(reps, g, chroms, bl, window_bp,
                                              max(3, min_cuts // 2), 0, f"rep{k+1}")
    rows, key2i, membership = [], {}, defaultdict(set)
    for name, ws in sets.items():
        cap = max_windows if name == "pooled" else max_windows // 2
        ws = sorted(ws, key=lambda w: (-w["score"], w["chrom"], w["start"]))[:cap]
        for w in ws:
            key = (w["chrom"], w["start"])
            if key not in key2i:
                key2i[key] = len(rows); rows.append(w)
            membership[key2i[key]].add(name)
    for i, w in enumerate(rows):
        w["sets"] = ",".join(sorted(membership[i]))
    write_windows(out, rows, len(reps), g)
    return {"label": label, "fragments": info, "n_windows": {k: sum(1 for m in membership.values() if k in m) for k in sets},
            "params": dict(window_bp=window_bp, max_windows=max_windows, min_cuts=min_cuts,
                           min_rep_cuts=min_rep_cuts, chroms=chroms)}


def write_windows(out, rows, n_reps, genome=None):
    if genome is not None:
        for w in rows:
            rf = genome.repeat_frac(w["chrom"], int(w["start"]), int(w["end"]))
            w["repeat_frac"] = "" if rf != rf else round(rf, 3)
    rows = sorted(rows, key=lambda w: (MAIN_CHROMS.index(w["chrom"]) if w["chrom"] in MAIN_CHROMS else 99, w["start"]))
    fields = ["window_id", "chrom", "start", "end", "n_cuts", "gc", "repeat_frac", "sets"] + [f"rep{k+1}_cuts" for k in range(n_reps)]
    for w in rows:
        w["window_id"] = f"{w['chrom']}:{w['start']}-{w['end']}"
        for k in range(n_reps):
            w[f"rep{k+1}_cuts"] = w["rep_cuts"][k]
    write_tsv(out, rows, fields)


# ---------------------------------------------------------------- windows from Corpus cell x cCRE
def windows_from_corpus(h5ad, genome_dir, blacklist, outdir, celltype_col="Cell_type (HSC)",
                        window_bp=300, max_windows=20000, min_cells=50, seed=7, label=""):
    """Per cell type: pseudobulk counts per cCRE, keep cCREs ranked by CPM (>=3 cells open);
    two random half-pseudoreplicates give rep1/rep2 window sets (replicate-split test)."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    a = ad.read_h5ad(h5ad)
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    names = np.asarray(a.var_names)
    m = np.array([re.match(r"^(chr[0-9XY]+)[:_-](\d+)[-_](\d+)$", n) is not None for n in names])
    parsed = [re.match(r"^(chr[0-9XY]+)[:_-](\d+)[-_](\d+)$", n) for n in names]
    chrom = np.array([p.group(1) if p else "" for p in parsed])
    mid = np.array([(int(p.group(2)) + int(p.group(3))) // 2 if p else -1 for p in parsed])
    g = Genome(genome_dir); bl = read_blacklist(blacklist); have = set(g.chroms())
    ok = m & np.isin(chrom, list(have))
    rng = np.random.default_rng(seed); outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    if celltype_col in a.obs.columns:
        cts = a.obs[celltype_col].astype(str)
    else:  # single-population sample (e.g. a cell line): treat all cells as one group
        cts = pd.Series(["all_cells"] * a.n_obs, index=a.obs_names)
    results = []
    for ct, n in cts.value_counts().items():
        if n < min_cells:
            continue
        idx = np.flatnonzero((cts == ct).to_numpy()); rng.shuffle(idx)
        halves = [idx[: len(idx) // 2], idx[len(idx) // 2:]]
        sub = X[idx]; tot = np.asarray(sub.sum(0)).ravel(); ncell = np.asarray((sub > 0).sum(0)).ravel()
        hs = [np.asarray(X[h].sum(0)).ravel() for h in halves]
        sets = {}
        for name, vec in (("pooled", tot), ("rep1", hs[0]), ("rep2", hs[1])):
            cand = np.flatnonzero(ok & (ncell >= 3) & (vec > 0))
            order = cand[np.argsort(-vec[cand], kind="stable")]
            cap = max_windows if name == "pooled" else max_windows // 2
            chosen = []
            for j in order:
                c = chrom[j]; s = int(mid[j]) - window_bp // 2; e = s + window_bp
                if s < 0 or any(b0 < e and s < b1 for b0, b1 in bl.get(c, [])):
                    continue
                seq = np.asarray(g.get(c)[s:e])
                if len(seq) < window_bp or (seq == 4).sum() > window_bp * 0.05:
                    continue
                chosen.append(j)
                if len(chosen) >= cap:
                    break
            sets[name] = chosen
        rows, key2i, mem = [], {}, defaultdict(set)
        for name, js in sets.items():
            for j in js:
                if j not in key2i:
                    c = chrom[j]; s = int(mid[j]) - window_bp // 2
                    seq = np.asarray(g.get(c)[s:s + window_bp])
                    gc = float(((seq == 1) | (seq == 2)).sum() / max(1, (seq < 4).sum()))
                    key2i[j] = len(rows)
                    rows.append({"chrom": c, "start": s, "end": s + window_bp, "n_cuts": int(tot[j]),
                                 "rep_cuts": [int(hs[0][j]), int(hs[1][j])], "gc": round(gc, 4), "ccre": names[j]})
                mem[key2i[j]].add(name)
        for i, w in enumerate(rows):
            w["sets"] = ",".join(sorted(mem[i]))
        slug = re.sub(r"[^A-Za-z0-9]+", "_", ct).strip("_")
        write_windows(outdir / f"{slug}.windows.tsv", rows, 2, g)
        results.append({"cell_type": ct, "slug": slug, "n_cells": int(n),
                        "n_windows": {k: len(v) for k, v in sets.items()}})
    meta = {"h5ad": str(h5ad), "h5ad_sha256": sha256(h5ad), "label": label, "celltype_col": celltype_col,
            "cell_types": results, "window_bp": window_bp, "max_windows": max_windows,
            "pseudoreplicates": "random half-splits of cells (seed %d)" % seed}
    (outdir / "corpus_windows.json").write_text(json.dumps(meta, indent=1))
    return meta


# ---------------------------------------------------------------- motifs
def read_meme(path, ids=None):
    motifs, lines, i = [], Path(path).read_text().splitlines(), 0
    while i < len(lines):
        if not lines[i].startswith("MOTIF "):
            i += 1; continue
        f = lines[i].split(); mid = f[1]; name = f[2] if len(f) > 2 else mid
        i += 1
        while i < len(lines) and not lines[i].startswith("letter-probability"):
            i += 1
        w = int(lines[i].split("w=")[1].split()[0])
        pwm = np.array([[float(x) for x in lines[i + 1 + j].split()[:4]] for j in range(w)])
        i += w + 1
        if ids is None or mid in ids or name in ids:
            motifs.append({"id": mid, "name": name, "pwm": pwm})
    return motifs


def score_threshold(weight, bg, pvalue):
    """Exact score threshold for P(score >= t) <= pvalue under an iid background (DP)."""
    scale = 100.0
    iw = np.round(weight * scale).astype(np.int64)
    lo = int(iw.min(1).sum()); hi = int(iw.max(1).sum())
    dist = np.array([1.0])  # index = score - cur_min
    cur_min = 0
    for row in iw:
        rmin = int(row.min()); new = np.zeros(len(dist) + int(row.max() - rmin))
        for b in range(4):
            sh = int(row[b] - rmin); new[sh:sh + len(dist)] += dist * bg[b]
        dist = new; cur_min += rmin
    tail = np.cumsum(dist[::-1])[::-1]
    k = np.flatnonzero(tail <= pvalue)
    t_int = (k[0] + cur_min) if len(k) else (len(dist) - 1 + cur_min)
    return t_int / scale


def build_weights(motifs, bg, pvalue, pseudo=0.01):
    L = max(m["pwm"].shape[0] for m in motifs)
    W = np.full((len(motifs), L, 5), 0.0); lens = np.zeros(len(motifs), np.int64); th = np.zeros(len(motifs))
    for k, m in enumerate(motifs):
        p = (m["pwm"] + pseudo) / (m["pwm"] + pseudo).sum(1, keepdims=True)
        w = np.log2(p / bg[None, :])
        W[k, :len(w), :4] = w; W[k, :len(w), 4] = -1000.0; lens[k] = len(w)
        th[k] = score_threshold(w, bg, pvalue)
    return W, lens, th


def _scan_kernel():
    import numba as nb

    @nb.njit(parallel=True, cache=True)
    def scan(codes, W, lens, th, count_only, out_w, out_m, out_p, out_s, out_sc, offsets):
        n, wl = codes.shape; M = W.shape[0]
        counts = np.zeros(M, np.int64)
        for m in nb.prange(M):
            L = lens[m]; c = 0; o = offsets[m]
            for i in range(n):
                for p in range(wl - L + 1):
                    sf = 0.0; sr = 0.0
                    for l in range(L):
                        b = codes[i, p + l]
                        sf += W[m, l, b]
                        rb = codes[i, p + L - 1 - l]
                        sr += W[m, l, 4 if rb == 4 else 3 - rb]
                    best = sf if sf >= sr else sr
                    if best >= th[m]:
                        if not count_only:
                            out_w[o + c] = i; out_m[o + c] = m; out_p[o + c] = p
                            out_s[o + c] = 1 if sf >= sr else -1; out_sc[o + c] = best
                        c += 1
            counts[m] = c
        return counts
    return scan


def scan_windows(windows_tsv, genome_dir, meme, out_npz, pvalue=1e-4, panel=None):
    rows = read_tsv(windows_tsv); g = Genome(genome_dir)
    wl = int(rows[0]["end"]) - int(rows[0]["start"])
    codes = np.stack([np.asarray(g.get(r["chrom"])[int(r["start"]):int(r["end"])]) for r in rows]).astype(np.int64)
    comp = np.array([(codes == b).sum() for b in range(4)], float); comp /= comp.sum()
    motifs = read_meme(meme, set(panel) if panel else None)
    W, lens, th = build_weights(motifs, comp, pvalue)
    scan = _scan_kernel(); M = len(motifs); dummy = np.zeros(1, np.int64)
    counts = scan(codes, W, lens, th, True, dummy, dummy, dummy, dummy, np.zeros(1), np.zeros(M, np.int64))
    offsets = np.r_[0, np.cumsum(counts)[:-1]]; tot = int(counts.sum())
    ow, om, op, os_ = (np.zeros(tot, np.int64) for _ in range(4)); osc = np.zeros(tot)
    scan(codes, W, lens, th, False, ow, om, op, os_, osc, offsets)
    np.savez_compressed(out_npz, window=ow.astype(np.int32), motif=om.astype(np.int16), pos=op.astype(np.int16),
                        strand=os_.astype(np.int8), score=osc.astype(np.float32),
                        motif_ids=np.array([m["id"] for m in motifs]), motif_names=np.array([m["name"] for m in motifs]),
                        motif_len=lens, thresholds=th, composition=comp, pvalue=pvalue)
    return {"n_windows": len(rows), "n_motifs": M, "n_hits": tot}


# ---------------------------------------------------------------- statistics
def bh(p):
    p = np.asarray(p, float); n = len(p)
    if not n:
        return p
    o = np.argsort(p); r = p[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(r, 1); return out


def _qbin(x, k):
    if k <= 1:
        return np.zeros(len(x), int)
    return np.searchsorted(np.unique(np.quantile(x, np.arange(1, k) / k)), x, side="left")


def strata(gc, cuts, rep=None, density=None, n_gc=5):
    """Matched-background strata: GC quintiles x accessibility tertiles (pilot used 3x3)
    x motif-density tertiles (number of distinct motifs hit; controls for generically
    'motif-rich'/AT-rich windows) x repeat bin (>=50% soft-masked) when available."""
    st = _qbin(gc, n_gc) * 3 + _qbin(cuts, 3)
    k = n_gc * 3
    if density is not None:
        st = st * 3 + _qbin(density, 3); k *= 3
    if rep is not None and np.isfinite(rep).all():
        st = st + k * (rep >= 0.5)
    return st


def stratified_pair_stats(P, st):
    """P: n x M bool presence. Returns obs, expected, var (sum of per-stratum hypergeometrics)."""
    M = P.shape[1]; obs = np.zeros((M, M)); exp = np.zeros((M, M)); var = np.zeros((M, M))
    for s in np.unique(st):
        Q = P[st == s].astype(np.float32); n = Q.shape[0]
        if n < 2:
            continue
        na = Q.sum(0).astype(float)
        obs += (Q.T @ Q)
        e = np.outer(na, na) / n; exp += e
        var += e * np.outer(n - na, n - na) / (n * (n - 1))
    return obs, exp, var


def pwm_similarity(a, b, min_overlap=6):
    """Max Pearson correlation of aligned PWM columns over offsets and both strands."""
    best = 0.0
    for bb in (b, b[::-1, ::-1]):
        for off in range(-(len(bb) - min_overlap), len(a) - min_overlap + 1):
            i0, j0 = max(0, off), max(0, -off); L = min(len(a) - i0, len(bb) - j0)
            if L < min_overlap:
                continue
            x = a[i0:i0 + L].ravel() - 0.25; y = bb[j0:j0 + L].ravel() - 0.25
            d = np.sqrt((x * x).sum() * (y * y).sum())
            if d > 0:
                best = max(best, float((x * y).sum() / d))
    return best


def analyse(windows_tsv, hits_npz, meme, outdir, label, dataset, context, min_windows_with_motif=20,
            permutations=999, top_perm=300, max_gap=150, q_cut=0.05, min_lift=1.25, seed=42, min_lift_matched=1.05,
            max_set_motifs=40, n_examples=25):
    """Pairs + replicate split + spacing + ambiguity + triplets + modules for one sample."""
    from scipy.stats import norm, poisson
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(windows_tsv); H = np.load(hits_npz)
    n_all = len(rows); names = [str(x) for x in H["motif_names"]]; ids = [str(x) for x in H["motif_ids"]]
    M = len(names)
    P_all = np.zeros((n_all, M), bool); P_all[H["window"], H["motif"]] = True
    insets = {k: np.array([k in r["sets"].split(",") for r in rows]) for k in ("pooled", "rep1", "rep2")}
    gc = np.array([float(r["gc"]) for r in rows]); cuts = np.array([float(r["n_cuts"]) for r in rows])
    rep = np.array([float(r["repeat_frac"]) if r.get("repeat_frac") not in (None, "") else np.nan for r in rows])
    rep_cut_cols = [k for k in rows[0] if re.match(r"rep\d+_cuts", k)]
    res = {}
    for k, msk in insets.items():
        if msk.sum() < 50:
            continue
        P = P_all[msk]; st = strata(gc[msk], cuts[msk], rep[msk], P.sum(1) if M >= 50 else None)
        obs, exp, var = stratified_pair_stats(P, st)
        z = (obs - exp - 0.5) / np.sqrt(np.maximum(var, 1e-9))
        res[k] = dict(P=P, st=st, obs=obs, exp=exp, z=z, p=norm.sf(z), n=int(msk.sum()))
    pooled = res["pooled"]; P = pooled["P"]; na = P.sum(0)
    # pilot-style background (GC x accessibility x repeat, no motif-density matching) for effect size
    st_basic = strata(gc[insets["pooled"]], cuts[insets["pooled"]], rep[insets["pooled"]])
    ob_b, ex_b, va_b = stratified_pair_stats(P, st_basic)
    keep = np.flatnonzero(na >= min_windows_with_motif)
    iu = [(a, b) for ai, a in enumerate(keep) for b in keep[ai + 1:]]
    A = np.array([x[0] for x in iu]); B = np.array([x[1] for x in iu])
    p = pooled["p"][A, B]; q = bh(p)
    lift = pooled["obs"][A, B] / np.maximum(pooled["exp"][A, B], 1e-9)
    lift_b = ob_b[A, B] / np.maximum(ex_b[A, B], 1e-9)
    z_b = (ob_b[A, B] - ex_b[A, B] - 0.5) / np.sqrt(np.maximum(va_b[A, B], 1e-9)); q_b = bh(norm.sf(z_b))
    # stratified permutations (pilot-compatible) for the strongest pairs
    rng = np.random.default_rng(seed); order = np.argsort(p)[:top_perm]
    perm_p = np.full(len(A), np.nan); st = pooled["st"]; groups = [np.flatnonzero(st == s) for s in np.unique(st)]
    if permutations and len(order):
        cols = sorted(set(A[order]) | set(B[order])); ci = {c: i for i, c in enumerate(cols)}
        Pc = P[:, cols].astype(np.float32); ge = np.zeros(len(order))
        obs_o = pooled["obs"][A[order], B[order]]
        ia = np.array([ci[a] for a in A[order]]); ib = np.array([ci[b] for b in B[order]])
        for _ in range(permutations):
            Sh = Pc.copy()
            for gidx in groups:
                Sh[gidx] = Sh[rng.permutation(gidx)]  # shuffle rows within stratum ...
            # ... for the B motif only: co-occurrence of A(original) with B(shuffled)
            co = (Pc[:, ia] * Sh[:, ib]).sum(0)
            ge += co >= obs_o
        perm_p[order] = (ge + 1) / (permutations + 1)
    # hit lists for spacing / examples
    win, mot, pos, ln = H["window"], H["motif"], H["pos"], H["motif_len"]
    pooled_idx = np.flatnonzero(insets["pooled"]); remap = -np.ones(n_all, np.int64); remap[pooled_idx] = np.arange(len(pooled_idx))
    by = defaultdict(list)
    sel = remap[win] >= 0
    for w_, m_, p_ in zip(remap[win][sel], mot[sel], pos[sel]):
        by[(int(w_), int(m_))].append(int(p_))
    rep_pooled = rep[insets["pooled"]]
    pwms = {m["id"]: m["pwm"] for m in read_meme(meme, set(ids))}
    fam = {k: tf_family(names[k]) for k in keep}
    out_rows = []
    cand = np.flatnonzero(((q <= 0.25) & (lift >= 1.0)) | ((q_b <= 0.05) & (lift_b >= min_lift)))
    cand_set = set(cand.tolist())
    for j in range(len(A)):
        a, b = int(A[j]), int(B[j])
        r = {"dataset": dataset, "sample": label, "context": context,
             "motif_a": ids[a], "motif_b": ids[b], "tf_a": names[a], "tf_b": names[b],
             "n_windows": pooled["n"], "n_a": int(na[a]), "n_b": int(na[b]),
             "n_both": int(pooled["obs"][a, b]), "expected_both": round(float(pooled["exp"][a, b]), 3),
             "lift": round(float(lift[j]), 4), "z": round(float(pooled["z"][a, b]), 3),
             "p": float(p[j]), "q": float(q[j]),
             "lift_basic": round(float(lift_b[j]), 4), "q_basic": float(q_b[j]),
             "perm_p": "" if np.isnan(perm_p[j]) else float(perm_p[j])}
        for k in ("rep1", "rep2"):
            if k in res:
                rr = res[k]; e = rr["exp"][a, b]
                r[f"{k}_lift"] = round(float(rr["obs"][a, b] / e), 4) if e > 0 else ""
                r[f"{k}_p"] = float(rr["p"][a, b])
        r["replicate_support"] = ("yes" if all(isinstance(r.get(f"{k}_lift"), float) and r[f"{k}_lift"] > 1 and r[f"{k}_p"] < 0.05 for k in ("rep1", "rep2"))
                                  else "no") if "rep1" in res and "rep2" in res else "n/a"
        if j in cand_set:
            close, gaps, overlap = 0, [], 0
            both = np.flatnonzero(P[:, a] & P[:, b])
            for w_ in both:
                pa, pb = by[(int(w_), a)], by[(int(w_), b)]
                la, lb = int(ln[a]), int(ln[b])
                gs = [max(x, y) - min(x + la, y + lb) for x in pa for y in pb]
                acc = [g_ for g_ in gs if 0 <= g_ <= max_gap]
                if acc:
                    close += 1; gaps.append(min(acc))
                if any(g_ < 0 for g_ in gs):
                    overlap += 1
            sim = pwm_similarity(pwms[ids[a]], pwms[ids[b]])
            shared_fam = sorted(set(fam[a]) & set(fam[b]))
            r.update(n_close_nonoverlap=close, median_gap_bp=float(np.median(gaps)) if gaps else "",
                     frac_overlapping_sites=round(overlap / max(1, len(both)), 3), pwm_similarity=round(sim, 3))
            flags = []
            if shared_fam: flags.append("same_family:" + "/".join(shared_fam))
            if sim >= 0.75: flags.append(f"similar_pwm:{sim:.2f}")
            if r["frac_overlapping_sites"] >= 0.5: flags.append("overlapping_sites")
            if "::" in names[a] or "::" in names[b]: flags.append("dimer_motif")
            rj = rep_pooled[both] if len(both) else np.array([np.nan])
            r["joint_repeat_frac"] = "" if not np.isfinite(rj).all() else round(float(np.mean(rj >= 0.5)), 3)
            if r["joint_repeat_frac"] != "" and r["joint_repeat_frac"] >= 0.5: flags.append("repeat_rich")
            r["ambiguity"] = ";".join(flags)
            distinct = not any(f.split(":")[0] in ("same_family", "similar_pwm", "overlapping_sites", "repeat_rich") for f in flags)
            passes = (q[j] <= q_cut and lift[j] >= min_lift_matched and lift_b[j] >= min_lift
                      and r["n_both"] >= 10 and close >= 5)
            r["status"] = ("priority_distinct" if passes and distinct else
                           "family_level" if passes else "exploratory")
        else:
            r.update(n_close_nonoverlap="", median_gap_bp="", frac_overlapping_sites="", pwm_similarity="",
                     ambiguity="", joint_repeat_frac="", status="not_enriched")
        r["families_a"] = "/".join(fam[a]); r["families_b"] = "/".join(fam[b])
        out_rows.append(r)
    rank = {"priority_distinct": 0, "family_level": 1, "exploratory": 2, "not_enriched": 3}
    out_rows.sort(key=lambda r: (rank[r["status"]], r["q"], -r["lift"]))
    fields = ["dataset", "sample", "context", "tf_a", "tf_b", "motif_a", "motif_b", "families_a", "families_b",
              "n_windows", "n_a", "n_b", "n_both", "expected_both", "lift", "z", "p", "q", "lift_basic", "q_basic", "perm_p",
              "rep1_lift", "rep1_p", "rep2_lift", "rep2_p", "replicate_support", "n_close_nonoverlap",
              "median_gap_bp", "frac_overlapping_sites", "pwm_similarity", "joint_repeat_frac", "ambiguity", "status"]
    write_tsv(outdir / "pairs.tsv", [r for r in out_rows if r["status"] != "not_enriched"], fields)
    # site examples
    ex = []
    prow = [r for r in out_rows if r["status"] in ("priority_distinct", "family_level")][:60]
    idmap = {i: k for k, i in enumerate(ids)}
    for r in prow:
        a, b = idmap[r["motif_a"]], idmap[r["motif_b"]]; la, lb = int(ln[a]), int(ln[b]); got = []
        for w_ in np.flatnonzero(P[:, a] & P[:, b]):
            gs = [(max(x, y) - min(x + la, y + lb), x, y) for x in by[(int(w_), a)] for y in by[(int(w_), b)]]
            acc = [t for t in gs if 0 <= t[0] <= max_gap]
            if acc:
                gap, x, y = min(acc); wr = rows[pooled_idx[w_]]; s0 = int(wr["start"])
                got.append({"tf_a": r["tf_a"], "tf_b": r["tf_b"], "window_id": wr["window_id"], "chrom": wr["chrom"],
                            "a_start": s0 + x, "a_end": s0 + x + la, "b_start": s0 + y, "b_end": s0 + y + lb,
                            "gap_bp": gap, "n_cuts": int(float(wr["n_cuts"])),
                            "rep_cuts": "/".join(wr[c] for c in rep_cut_cols)})
        got.sort(key=lambda t: -t["n_cuts"]); ex += got[:n_examples]
    write_tsv(outdir / "pair_site_examples.tsv", ex,
              ["tf_a", "tf_b", "window_id", "chrom", "a_start", "a_end", "b_start", "b_end", "gap_bp", "n_cuts", "rep_cuts"])
    # ---- higher-order sets: triplets beyond best pair, among motifs of significant pairs
    sig = [r for r in out_rows if r["status"] in ("priority_distinct", "family_level")]
    mset = []
    for r in sig:
        if r["status"] != "priority_distinct" or r["replicate_support"] == "no":
            continue
        for mid in (r["motif_a"], r["motif_b"]):
            k = idmap[mid]
            if k in mset or len(mset) >= max_set_motifs:
                continue
            # keep one representative per family / PWM-similarity cluster so sets are not homolog stacks
            if any((set(fam[k]) & set(fam[o])) or pwm_similarity(pwms[ids[k]], pwms[ids[o]]) >= 0.75 for o in mset):
                continue
            mset.append(k)
    trip = []
    st = pooled["st"]
    for a, b, c in combinations(sorted(mset), 3):
        pa, pb, pc = P[:, a], P[:, b], P[:, c]
        o = int((pa & pb & pc).sum())
        if o < 5:
            continue
        exps = []
        for x, y, zz in ((pa, pb, pc), (pa, pc, pb), (pb, pc, pa)):
            e = 0.0
            for s in np.unique(st):
                m_ = st == s; ns = m_.sum()
                e += (x[m_] & y[m_]).sum() * zz[m_].sum() / ns
            exps.append(e)
        e_ind = sum(((pa[st == s].mean() * pb[st == s].mean() * pc[st == s].mean()) * (st == s).sum()) for s in np.unique(st))
        e_cons = max(exps)
        pv = float(poisson.sf(o - 1, e_cons)) if e_cons > 0 else 1.0
        trip.append({"dataset": dataset, "sample": label, "tfs": "+".join(sorted([names[a], names[b], names[c]])),
                     "motifs": "+".join([ids[a], ids[b], ids[c]]), "n_all3": o,
                     "expected_independent": round(e_ind, 3), "lift_vs_independent": round(o / e_ind, 3) if e_ind else "",
                     "expected_given_best_pair": round(e_cons, 3), "synergy_lift": round(o / e_cons, 3) if e_cons else "",
                     "p_synergy": pv})
    tq = bh([t["p_synergy"] for t in trip]) if trip else []
    for t, qq in zip(trip, tq):
        t["q_synergy"] = float(qq)
        t["status"] = "synergistic_set" if qq <= 0.05 and t["synergy_lift"] and t["synergy_lift"] >= 1.2 else "additive_set"
    trip.sort(key=lambda t: (t["status"] != "synergistic_set", t["q_synergy"], -(t["lift_vs_independent"] or 0)))
    write_tsv(outdir / "motif_sets.tsv", trip, ["dataset", "sample", "tfs", "motifs", "n_all3", "expected_independent",
              "lift_vs_independent", "expected_given_best_pair", "synergy_lift", "p_synergy", "q_synergy", "status"])
    # modules: maximal cliques (size >= 3) in the graph of distinct significant pairs
    adj = defaultdict(set)
    top_d = [r for r in sig if r["status"] == "priority_distinct" and r["replicate_support"] != "no"][:400]
    for r in top_d:
        adj[r["tf_a"]].add(r["tf_b"]); adj[r["tf_b"]].add(r["tf_a"])
    cliques = []
    stack = [(frozenset(), frozenset(adj), frozenset())]
    while stack and len(cliques) < 500:  # iterative Bron-Kerbosch with pivot
        R_, P_, X_ = stack.pop()
        if not P_ and not X_:
            if len(R_) >= 3: cliques.append(sorted(R_))
            continue
        piv = max(P_ | X_, key=lambda v: len(adj[v] & P_))
        for v in list(P_ - adj[piv]):
            stack.append((R_ | {v}, P_ & adj[v], X_ & adj[v])); P_ = P_ - {v}; X_ = X_ | {v}
    mods = []
    for cl in sorted(cliques, key=lambda c: -len(c))[:100]:
        cols = [names.index(t) for t in cl]
        nall = int(np.all(P[:, cols], axis=1).sum())
        mods.append({"dataset": dataset, "sample": label, "module": "+".join(cl), "size": len(cl), "n_windows_all": nall})
    write_tsv(outdir / "modules.tsv", mods, ["dataset", "sample", "module", "size", "n_windows_all"])
    summ = {"dataset": dataset, "sample": label, "context": context, "n_windows_total": n_all,
            "n_windows": {k: v["n"] for k, v in res.items()}, "n_motifs_scanned": M,
            "n_motifs_tested": int(len(keep)), "n_hits": int(len(win)), "n_pairs_tested": int(len(A)),
            "n_priority_distinct": sum(r["status"] == "priority_distinct" for r in out_rows),
            "n_family_level": sum(r["status"] == "family_level" for r in out_rows),
            "n_replicate_supported_sig": sum(r["status"] in ("priority_distinct", "family_level") and r["replicate_support"] == "yes" for r in out_rows),
            "n_sets_tested": len(trip), "n_synergistic_sets": sum(t["status"] == "synergistic_set" for t in trip),
            "n_modules": len(mods), "params": dict(min_windows_with_motif=min_windows_with_motif, permutations=permutations,
            top_perm=top_perm, max_gap=max_gap, q_cut=q_cut, min_lift=min_lift, min_lift_matched=min_lift_matched,
            background="lift/q: GC quintile x accessibility tertile x motif-density tertile x repeat; lift_basic/q_basic: same without motif density", seed=seed, motif_pvalue=float(H["pvalue"])),
            "label_note": "Motif-sequence co-occurrence hypotheses; not ChromBPNet, not occupancy, not physical interaction."}
    (outdir / "summary.json").write_text(json.dumps(summ, indent=1))
    return summ
