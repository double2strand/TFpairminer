import collections, json
import numpy as np
import pytest
from motif_cooccurrence import chrombpnet as CB


def test_dinuc_shuffle_preserves_dinucleotides():
    rng = np.random.default_rng(1)
    for _ in range(100):
        s = rng.integers(0, 4, rng.integers(3, 400))
        y = CB.dinuc_shuffle(CB.onehot(s), rng).argmax(1)
        d = lambda q: collections.Counter(zip(q[:-1], q[1:]))
        assert d(list(s)) == d(list(y)) and y[0] == s[0]


def test_onehot_and_shard_slice():
    x = CB.onehot(np.array([0, 1, 2, 3, 4]))
    assert x.sum() == 4 and x[4].sum() == 0
    n = 103; got = [i for k in range(4) for i in range(n)[CB.shard_slice(n, k, 4)]]
    assert got == list(range(n))


def test_win_offset_centres_window():
    r = {"start": 1000, "end": 1300}
    assert CB.win_offset(r) == CB.INLEN // 2 - 150


def test_merge_bh_and_labels(tmp_path):
    res = [{"tfs": "A+B", "motifs": "m1+m2", "k": 2, "status": "tested", "synergy": 0.2, "synergy_p": 1e-6,
            "contribution_supported": "yes"},
           {"tfs": "C+D", "motifs": "m3+m4", "k": 2, "status": "tested", "synergy": -0.1, "synergy_p": 0.5,
            "contribution_supported": "no"},
           {"tfs": "E+F", "motifs": "m5+m6", "k": 2, "status": "too_few_windows"}]
    (tmp_path / "pairs_cbp.shard0.json").write_text(json.dumps({"results": res, "tracks": []}))
    CB.merge(tmp_path)
    s = json.loads((tmp_path / "cbp_summary.json").read_text())
    assert s == {"n_candidates": 3, "n_tested": 2, "n_contribution_supported": 1, "n_synergistic": 1, "n_subadditive": 0}


def test_tiny_model_gradients(tmp_path):
    tf = pytest.importorskip("tensorflow")
    inp = tf.keras.Input((CB.INLEN, 4)); h = tf.keras.layers.Conv1D(4, 9, activation="relu")(inp)
    prof = tf.keras.layers.Flatten()(tf.keras.layers.Conv1D(1, 1)(h[:, :1000]))
    cnt = tf.keras.layers.Dense(1)(tf.keras.layers.GlobalAveragePooling1D()(h))
    m = tf.keras.Model(inp, [prof, cnt]); p = tmp_path / "chrombpnet_nobias.h5"; m.save(p)
    km, grad, counts = CB.load_model(tmp_path)
    X = CB.onehot(np.random.default_rng(0).integers(0, 4, CB.INLEN))[None].repeat(3, 0)
    g = CB.batched(grad, X, 2); c = CB.batched(counts, X, 2)
    assert g.shape == X.shape and np.isfinite(g).all() and np.abs(g).sum() > 0 and c.shape == (3,)
