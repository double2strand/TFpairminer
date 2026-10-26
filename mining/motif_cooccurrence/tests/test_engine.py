import gzip, itertools, json
from pathlib import Path

import numpy as np
import pytest

from motif_cooccurrence import engine as E
from motif_cooccurrence import aggregate as AG

MEME = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF M1 AAA1
letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 0
""" + "\n".join(" ".join("0.97" if i == k else "0.01" for i in range(4)) for k in [0, 1, 2, 3, 0, 1, 2, 3]) + """

MOTIF M2 BBB2
letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 0
""" + "\n".join(" ".join("0.97" if i == k else "0.01" for i in range(4)) for k in [3, 3, 1, 0, 2, 2, 0, 1]) + """

MOTIF M3 CCC3
letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 0
""" + "\n".join(" ".join("0.97" if i == k else "0.01" for i in range(4)) for k in [2, 0, 3, 1, 1, 3, 0, 2]) + "\n"

S1, S2, S3 = "ACGTACGT", "TTCAGGAC", "GATCCTAG"


def make_genome(tmp, n=600, seed=0):
    rng = np.random.default_rng(seed); L = n * 1000 + 2000
    seq = rng.choice(list("ACGT"), L)
    peaks = []
    for i in range(n):
        c = 1000 + i * 1000; peaks.append(c)
        if i % 3 == 0:   # planted M1+M2 pair
            seq[c - 40:c - 32] = list(S1); seq[c + 20:c + 28] = list(S2)
        elif i % 3 == 1:  # M3 alone
            seq[c:c + 8] = list(S3)
    fa = tmp / "g.fa"
    fa.write_text(">chr1\n" + "".join(seq) + "\n")
    E.prep_genome(fa, tmp / "g", chroms=["chr1"])
    for rep in (1, 2):
        with gzip.open(tmp / f"r{rep}.tsv.gz", "wt") as fh:
            for c in peaks:
                for k in range(8):
                    s = int(c + rng.integers(-60, 60)); fh.write(f"chr1\t{s}\t{s + 80}\tBC{k}\t1\n")
    (tmp / "m.meme").write_text(MEME)
    return peaks


def test_threshold_matches_bruteforce():
    w = np.log2(np.array([[0.7, 0.1, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1], [0.25, 0.25, 0.25, 0.25]]) / 0.25)
    t = E.score_threshold(w, np.full(4, 0.25), 0.1)
    scores = sorted((sum(w[i, b] for i, b in enumerate(x)) for x in itertools.product(range(4), repeat=3)), reverse=True)
    frac = np.mean(np.array(scores) >= t - 1e-9)
    assert frac <= 0.1 + 1e-9 and frac > 0.02


def test_bh_monotone():
    q = E.bh(np.array([0.01, 0.04, 0.03, 0.5]))
    assert np.all(q >= np.array([0.01, 0.04, 0.03, 0.5]) - 1e-12) and q.max() <= 1


def test_strata_shapes():
    rng = np.random.default_rng(1)
    st = E.strata(rng.random(500), rng.random(500), rng.random(500), rng.integers(0, 50, 500))
    assert st.min() >= 0 and len(np.unique(st)) > 20


def test_end_to_end_planted_pair(tmp_path):
    make_genome(tmp_path)
    meta = E.windows_from_fragments([tmp_path / "r1.tsv.gz", tmp_path / "r2.tsv.gz"], tmp_path / "g", "",
                                    tmp_path / "w.tsv", window_bp=300, max_windows=1000, min_cuts=4, min_rep_cuts=1)
    assert meta["n_windows"]["pooled"] > 300
    E.scan_windows(tmp_path / "w.tsv", tmp_path / "g", tmp_path / "m.meme", tmp_path / "h.npz", 1e-4)
    s = E.analyse(tmp_path / "w.tsv", tmp_path / "h.npz", tmp_path / "m.meme", tmp_path / "o", "T", "TEST", "synthetic",
                  min_windows_with_motif=10, permutations=49)
    rows = E.read_tsv(tmp_path / "o" / "pairs.tsv")
    top = {(r["tf_a"], r["tf_b"]) for r in rows if r["status"] == "priority_distinct"}
    assert ("AAA1", "BBB2") in top or ("BBB2", "AAA1") in top
    assert not any("CCC3" in k for k in top)
    assert json.loads((tmp_path / "o" / "summary.json").read_text())["n_priority_distinct"] >= 1
    # aggregation into site JSON
    d = tmp_path / "run" / "samples" / "T"; d.mkdir(parents=True)
    for f in (tmp_path / "o").iterdir():
        (d / f.name).write_bytes(f.read_bytes())
    idx = AG.aggregate(tmp_path / "run", tmp_path / "site")
    assert idx["counts"]["pairs_listed"] >= 1 and (tmp_path / "site" / "tf" / "AAA1.json").exists()
