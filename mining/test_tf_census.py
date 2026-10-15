import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse

import tf_census as m


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.features = [
            {"ensembl_id": "ENSG1", "symbol": "TF1", "hgnc_id": "HGNC:1", "category": "TF", "is_tf": True},
            {"ensembl_id": "ENSG2", "symbol": "TF2", "hgnc_id": "HGNC:2", "category": "TF", "is_tf": True},
            {"ensembl_id": "ENSG3", "symbol": "CF", "hgnc_id": "HGNC:3", "category": "cofactor", "is_tf": False},
        ]
        self.meta = {"dataset_id": "dataset1", "collection_id": "collection1", **dict(zip(m.KEY_COLUMNS, ["CL:1", "cell", "UBERON:1", "tissue", "PATO:1", "normal", "EFO:1", "assay"]))}

    def test_register_reports_missing_ensembl(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "register.tsv"
            fields = ["Ensembl_gene_ID", "Symbol", "HGNC_ID", "Primary_category", "Strict_TF_set", "Exploratory_TF_set", "Expanded_regulator_set"]
            m.write_tsv(path, [dict(zip(fields, ["ENSG1", "TF1", "HGNC:1", "TF", "Yes", "Yes", "Yes"])),
                               dict(zip(fields, ["", "DUX1", "HGNC:2", "TF", "Yes", "Yes", "Yes"]))], fields)
            genes, missing = m.load_register(path, "exploratory")
            self.assertEqual(len(genes), 1)
            self.assertEqual(missing[0]["symbol"], "DUX1")

    def test_expression_and_pair_statistics(self):
        acc = m.Accumulator(self.features)
        x = sparse.csr_matrix([[1, 1, 0], [2, 2, 1], [0, 0, 0], [0, 0, 1]], dtype=float)
        acc.add(x[:2], [10, 20], ["d1", "d1"])
        acc.add(x[2:], [10, 10], ["d2", "d2"])
        profiles = list(acc.profiles(self.meta))
        self.assertEqual(profiles[0]["n_detected"], 2)
        self.assertAlmostEqual(profiles[0]["detection_fraction"], 0.5)
        self.assertAlmostEqual(profiles[0]["mean_log1p_cptt"], np.log1p(1000) / 2)
        self.assertEqual(profiles[2]["n_detected"], 2)
        pairs = list(acc.pairs(self.meta, 0.4, 0.5, 1.1, 2, 10))
        self.assertEqual(len(pairs), 1)
        self.assertEqual((pairs[0]["tf1_symbol"], pairs[0]["tf2_symbol"]), ("TF1", "TF2"))
        self.assertAlmostEqual(pairs[0]["cell_pearson_r"], 1.0)
        self.assertAlmostEqual(pairs[0]["co_detection_lift"], 2.0)

    def test_writer_roundtrip_and_aggregate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            m.write_tsv(root / "datasets.tsv", [{"dataset_id": f"dataset{i}"} for i in range(2)], ["dataset_id"])
            m.write_json_atomic(root / "manifest.json", {"dataset_selection": "all"})
            m.write_tsv(root / "features.tsv", self.features, ["ensembl_id", "symbol", "hgnc_id", "category", "is_tf"])
            priors = root / "interaction_priors.tsv"
            m.write_tsv(priors, [{"A_HGNC_ID": "HGNC:1", "B_HGNC_ID": "HGNC:2", "Source_pair_class": "TF-TF",
                                  "Best_source_evidence": "Direct-interaction annotation", "PMIDs": "12345"}],
                        ["A_HGNC_ID", "B_HGNC_ID", "Source_pair_class", "Best_source_evidence", "PMIDs"])
            for i in range(2):
                folder = root / "datasets" / f"dataset{i}"
                row_meta = {**self.meta, "dataset_id": f"dataset{i}", "collection_id": f"collection{i}"}
                profile = {**row_meta, "ensembl_id": "ENSG1", "symbol": "TF1", "hgnc_id": "HGNC:1", "category": "TF",
                           "n_cells": 100, "n_donors": 2, "n_detected": 60, "detection_fraction": 0.6,
                           "mean_log1p_cptt": 1.0 + i, "mean_raw_count": 2.0}
                pair = {**row_meta, "tf1_ensembl": "ENSG1", "tf1_symbol": "TF1", "tf2_ensembl": "ENSG2", "tf2_symbol": "TF2",
                        "n_cells": 100, "n_donors": 2, "tf1_detection_fraction": 0.6, "tf2_detection_fraction": 0.6,
                        "co_detected_cells": 45, "co_detection_lift": 1.25, "cell_pearson_r": 0.4 + 0.2 * i, "donor_pearson_r": ""}
                m.write_tsv(folder / "profiles.tsv.gz", [profile], m.PROFILE_COLUMNS)
                m.write_tsv(folder / "pairs.tsv.gz", [pair], m.PAIR_COLUMNS)
                (folder / "_SUCCESS").write_text('{"run_config": {"min_cells": 100}}')
            m.command_aggregate(SimpleNamespace(out=tmp, interaction_priors=str(priors)))
            with open(root / "aggregate" / "cell_type_profiles.tsv", newline="") as fh:
                profiles = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(profiles), 1)
            self.assertEqual(int(profiles[0]["n_collections"]), 2)
            self.assertAlmostEqual(float(profiles[0]["mean_log1p_cptt"]), 1.5)
            with open(root / "aggregate" / "cell_type_pair_candidates.tsv", newline="") as fh:
                pairs = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(pairs), 1)
            self.assertEqual(int(pairs[0]["n_collections_with_hit"]), 2)
            self.assertAlmostEqual(float(pairs[0]["cell_r_weighted_by_cells"]), 0.5)
            self.assertEqual(pairs[0]["prior_best_evidence"], "Direct-interaction annotation")
            self.assertIn('"complete": true', (root / "aggregate" / "summary.json").read_text())
            self.assertIn('"atlas_wide": true', (root / "aggregate" / "summary.json").read_text())
            self.assertEqual(int(pairs[0]["n_datasets_in_meta"]), 0)  # v1 rows carry no df
            self.assertTrue((root / "aggregate" / "positive_controls.tsv").exists())


class StatisticsTests(unittest.TestCase):
    def test_bh_matches_reference(self):
        p = np.array([0.01, 0.04, 0.03, 0.005, np.nan])
        q = m.bh_qvalues(p)
        # reference: R p.adjust(c(.01,.04,.03,.005), "BH") = 0.02 0.04 0.04 0.02
        np.testing.assert_allclose(q[:4], [0.02, 0.04, 0.04, 0.02])
        self.assertTrue(np.isnan(q[4]))

    def test_corr_pvalue_matches_scipy(self):
        rng = np.random.default_rng(0)
        a = rng.normal(size=50); b = a + rng.normal(size=50)
        r, p = m.sstats.pearsonr(a, b)
        self.assertAlmostEqual(float(m.corr_pvalues(np.array([r]), 48)[0]), p, places=10)

    def _sim(self, seed=1, donors=6, per=60):
        rng = np.random.default_rng(seed)
        feats = [{"ensembl_id": f"E{i}", "symbol": f"T{i}", "hgnc_id": f"H{i}", "category": "TF", "is_tf": True} for i in range(3)]
        rows, libs, dons = [], [], []
        for d in range(donors):
            shift = rng.normal(scale=2)             # donor effect shared by T0,T1 (confounder)
            for _ in range(per):
                lib = rng.integers(2000, 8000)
                latent = rng.normal()
                lam = np.exp([1 + shift + 0.8 * latent, 1 + shift + 0.8 * latent, 1 + rng.normal()]) * lib / 5000
                rows.append(rng.poisson(lam)); libs.append(lib); dons.append(f"d{d}")
        return feats, np.array(rows, dtype=float), np.array(libs, dtype=float), dons

    def test_resid_cell_equals_explicit_regression(self):
        feats, x, lib, dons = self._sim()
        acc = m.Accumulator(feats, metacell_size=10)
        for s in range(0, len(x), 70):
            acc.add(sparse.csr_matrix(x[s:s + 70]), lib[s:s + 70], dons[s:s + 70])
        y = np.log1p(x * (10000 / lib)[:, None])
        r_explicit, df_explicit = m.grouped_partial_corr(y, np.log(lib), np.array(dons))
        r_stream, df_stream = acc.resid_cell_corr()
        np.testing.assert_allclose(r_stream, r_explicit, atol=1e-8)
        self.assertEqual(df_stream, df_explicit)
        self.assertEqual(df_stream, len(x) - 6 - 1 - 1)

    def test_modes_fdr_and_positive_signal(self):
        feats, x, lib, dons = self._sim()
        acc = m.Accumulator(feats, metacell_size=10)
        acc.add(sparse.csr_matrix(x), lib, dons)
        meta = dict(zip(m.KEY_COLUMNS, ["CL:1"] * 8), dataset_id="d", collection_id="c")
        for mode in ["metacell", "pseudobulk", "resid_cell", "cell_raw"]:
            rows = list(acc.pairs(meta, 0.05, 0.2, 0.0, 0, 0, primary=mode, max_q=0.05))
            hit = {(r["tf1_symbol"], r["tf2_symbol"]) for r in rows if r["passed"]}
            self.assertIn(("T0", "T1"), hit, mode)          # true co-regulated pair recovered
            self.assertNotIn(("T0", "T2"), hit, mode)       # independent TF not called
            self.assertEqual(rows[0]["primary_stat"], mode)
            self.assertEqual(rows[0]["n_tested_pairs"], 3)
            for r in rows:
                self.assertLessEqual(r["p_value"], r["q_value_bh"] + 1e-12)

    def test_null_fdr_calibration(self):
        rng = np.random.default_rng(3)
        feats = [{"ensembl_id": f"E{i}", "symbol": f"T{i}", "hgnc_id": "", "category": "TF", "is_tf": True} for i in range(12)]
        x = rng.poisson(2.0, size=(400, 12)).astype(float)
        acc = m.Accumulator(feats, metacell_size=10)
        acc.add(sparse.csr_matrix(x), np.full(400, 5000.0), [f"d{i % 8}" for i in range(400)])
        meta = dict(zip(m.KEY_COLUMNS, ["x"] * 8), dataset_id="d", collection_id="c")
        rows = [r for r in acc.pairs(meta, 0.05, -1, 0, 0, 0, primary="resid_cell", max_q=0.05) if r["passed"]]
        self.assertLessEqual(len(rows), 2)

    def test_meta_analysis(self):
        k, rf, pf, rr, pr, _, tau2, i2 = m.meta_analyse([(0.5, 100), (0.5, 100)])
        self.assertEqual(k, 2); self.assertAlmostEqual(rf, 0.5); self.assertAlmostEqual(rr, 0.5)
        self.assertAlmostEqual(tau2, 0.0); self.assertLess(pf, 1e-6)
        k, rf, pf, rr, pr, _, tau2, i2 = m.meta_analyse([(0.8, 200), (-0.2, 200)])
        self.assertGreater(tau2, 0); self.assertGreater(i2, 0.9); self.assertGreater(pr, pf)
        self.assertEqual(m.meta_analyse([(0.3, None)])[0], 0)


if __name__ == "__main__":
    unittest.main()
