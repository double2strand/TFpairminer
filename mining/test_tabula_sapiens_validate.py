import tempfile
import unittest
from pathlib import Path

from scipy import sparse

import tabula_sapiens_validate as ts
import tf_census as base


class TabulaSapiensValidationTests(unittest.TestCase):
    def test_candidate_key_is_pair_order_invariant(self):
        row = {"cell_type_id": "CL:1", "tissue_id": "UBERON:1", "disease_id": "PATO:1",
               "assay_id": "EFO:1", "tf1_ensembl": "ENSG2", "tf2_ensembl": "ENSG1"}
        flipped = {**row, "tf1_ensembl": "ENSG1", "tf2_ensembl": "ENSG2"}
        self.assertEqual(ts.candidate_id(row), ts.candidate_id(flipped))

    def test_donor_supported_pair_and_low_expression(self):
        features = [{"ensembl_id": "ENSG1", "is_tf": True}, {"ensembl_id": "ENSG2", "is_tf": True}]
        acc = base.Accumulator(features, track_donor_covariance=True)
        x = sparse.csr_matrix([[1, 1], [2, 2], [3, 3], [4, 4],
                               [1, 1], [2, 2], [3, 3], [4, 4]])
        acc.add(x, [100] * 8, ["donor1"] * 4 + ["donor2"] * 4)
        thresholds = {"min_detection": 0.05, "min_rho": 0.25, "min_lift": 0.9,
                      "min_codetected": 4, "min_donor_cells": 4, "min_donors": 2,
                      "min_positive_donor_fraction": 0.75}
        metrics, status, _ = ts.pair_metrics(acc, 0, 1, thresholds)
        self.assertEqual(status, "replicated_donor_supported")
        self.assertEqual(metrics["ts_donors_tested"], 2)
        self.assertAlmostEqual(metrics["ts_cell_r"], 1.0)
        self.assertEqual(metrics["ts_positive_donor_fraction"], 1.0)
        _, low_status, _ = ts.pair_metrics(acc, 0, 1, {**thresholds, "min_codetected": 9})
        self.assertEqual(low_status, "inconclusive_low_expression")

    def test_completed_discovery_overlap_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base.write_tsv(root / "datasets.tsv", [{"dataset_id": "ts-all"}], ["dataset_id"])
            marker = root / "datasets" / "ts-all" / "_SUCCESS"
            base.write_json_atomic(marker, {"n_profile_rows": 10})
            with self.assertRaises(RuntimeError):
                ts.assert_no_discovery_overlap(root, {"ts-all"})
            base.write_json_atomic(marker, {"n_profile_rows": 0})
            self.assertEqual(ts.assert_no_discovery_overlap(root, {"ts-all"}), 1)


if __name__ == "__main__":
    unittest.main()
