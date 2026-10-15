#!/usr/bin/env python3
"""Small deterministic integration tests; no CATlas download or FIMO required."""

import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from catlas_pipeline import (command_aggregate, command_analyze, command_compile,
                             command_plan, read_tsv, write_json, write_tsv)


class PipelineTest(unittest.TestCase):
    def test_reviewed_exact_pair_screen(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            h5ad_path = root / "catlas.h5ad"
            ontology = root / "ontology.tsv"
            meme = root / "motifs.meme"
            candidates = root / "candidates.tsv"
            crosswalk = root / "crosswalk.tsv"
            out = root / "run"
            X = np.zeros((1, 100), dtype=np.uint8)
            X[0, :20] = 1
            obs = pd.DataFrame({"cell type": ["Test Cell"]}, index=["cluster_0"])
            var = pd.DataFrame({"chrom": ["chr1"] * 100,
                                "start": np.arange(100) * 200,
                                "end": np.arange(100) * 200 + 100,
                                "Class": ["Distal"] * 100},
                               index=[f"peak_{i}" for i in range(100)])
            ad.AnnData(X=X, obs=obs, var=var).write_h5ad(h5ad_path)
            write_tsv(ontology, [{"Cell type": "Test Cell", "Cell Ontology ID": "CL:1234567"}],
                      ["Cell type", "Cell Ontology ID"])
            meme.write_text("MEME version 5\n\nALPHABET= ACGT\n\n"
                            "MOTIF M1 TFA\nletter-probability matrix: alength= 4 w= 1\n0.25 0.25 0.25 0.25\n\n"
                            "MOTIF M2 TFB\nletter-probability matrix: alength= 4 w= 1\n0.25 0.25 0.25 0.25\n")
            write_tsv(candidates, [{"candidate_id": "pair1", "tf1_ensembl": "ENSGA", "tf1_symbol": "TFA",
                                    "tf2_ensembl": "ENSGB", "tf2_symbol": "TFB", "cell_type_id": "CL:1234567",
                                    "cell_type": "Test Cell", "tissue_id": "UBERON:1", "tissue": "Test Tissue",
                                    "validation_status": "replicated_donor_supported"}],
                      ["candidate_id", "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
                       "cell_type_id", "cell_type", "tissue_id", "tissue", "validation_status"])
            write_tsv(crosswalk, [{"catlas_cluster": "Test Cell", "cell_type_id": "CL:1234567",
                                   "tissue_id": "UBERON:1", "stage": "adult", "mapping_tier": "exact_reviewed",
                                   "reviewed": "Yes"}],
                      ["catlas_cluster", "cell_type_id", "tissue_id", "stage", "mapping_tier", "reviewed"])
            command_plan(Namespace(out=out, candidate_summary=candidates, catlas_h5ad=h5ad_path,
                                   ontology=ontology, meme=meme, cluster_map=crosswalk,
                                   include_status=["replicated_donor_supported"]))
            scan_dir = out / "scan"
            write_tsv(scan_dir / "0000.hits.tsv.gz",
                      [{"motif_id": mid, "peak_index": i} for mid in ("M1", "M2") for i in range(20)],
                      ["motif_id", "peak_index"])
            write_tsv(scan_dir / "0000.gc.tsv.gz",
                      [{"peak_index": i, "gc": 0.5} for i in range(100)], ["peak_index", "gc"])
            write_json(scan_dir / "0000.done.json", {"shard_index": 0, "num_shards": 1,
                                                      "fimo_p": 1e-4, "max_stored_scores": 5000000,
                                                      "n_distal_peaks": 100,
                                                      "n_fimo_sites": 40})
            command_compile(Namespace(out=out))
            command_analyze(Namespace(out=out, catlas_h5ad=h5ad_path, num_shards=1,
                                      shard_index=0, min_matched=10, min_motif_hits=5))
            command_aggregate(Namespace(out=out, q=0.05, min_joint=5))
            result = list(read_tsv(out / "aggregate" / "catlas_pair_results.tsv.gz"))
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0]["catlas_result"], "catlas_coenriched")
            self.assertEqual(int(result[0]["joint_fg_hits"]), 20)
            self.assertEqual(int(result[0]["joint_bg_hits"]), 0)


if __name__ == "__main__":
    unittest.main()
