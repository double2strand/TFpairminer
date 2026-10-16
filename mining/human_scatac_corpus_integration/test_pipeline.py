#!/usr/bin/env python3
"""Synthetic RNA-first corpus test; no portal download or FIMO needed."""

import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from catlas_tf_integration.catlas_pipeline import read_tsv, write_json, write_tsv
from human_scatac_corpus_integration.pipeline import (
    command_aggregate, command_analyze, command_compile, command_index, command_plan,
)


class CorpusPipelineTest(unittest.TestCase):
    def test_independent_source_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data_dir = root / "data"
            data_dir.mkdir()
            catalog = root / "catalog.tsv"
            cell_index = root / "celltypes.tsv"
            crosswalk = root / "crosswalk.tsv"
            candidates = root / "candidates.tsv"
            meme = root / "motifs.meme"
            out = root / "run"
            x = np.zeros((50, 100), dtype=np.uint8)
            x[:, :20] = 1
            obs = pd.DataFrame({"Cell_type (HSC)": ["Test Cell"] * 50,
                                "Donor (HSC)": ["donor1"] * 25 + ["donor2"] * 25,
                                "Batch (HSC)": ["batch1"] * 50},
                               index=[f"cell_{i}" for i in range(50)])
            var = pd.DataFrame(index=[f"chr1:{i * 200}-{i * 200 + 100}" for i in range(100)])
            ad.AnnData(X=x, obs=obs, var=var).write_h5ad(data_dir / "1.h5ad")
            write_tsv(catalog, [{"sample_id": "1", "sample": "Independent2025-test",
                                 "dataset": "Independent2025", "tissue": "Test Tissue", "organ": "Test Organ",
                                 "pmid": "123", "overlaps_catlas_source": "No"}],
                      ["sample_id", "sample", "dataset", "tissue", "organ", "pmid",
                       "overlaps_catlas_source"])
            command_index(Namespace(catalog=catalog, data_dir=data_dir, out=cell_index))
            write_tsv(crosswalk, [{"cell_type_id": "CL:1234567", "tissue_id": "UBERON:1",
                                   "sample_id": "1", "hsc_cell_type": "Test Cell",
                                   "mapping_tier": "exact_reviewed", "reviewed": "Yes"}],
                      ["cell_type_id", "tissue_id", "sample_id", "hsc_cell_type",
                       "mapping_tier", "reviewed"])
            write_tsv(candidates, [{"candidate_id": "pair1", "tf1_ensembl": "ENSGA", "tf1_symbol": "TFA",
                                    "tf2_ensembl": "ENSGB", "tf2_symbol": "TFB", "cell_type_id": "CL:1234567",
                                    "cell_type": "Test Cell", "tissue_id": "UBERON:1", "tissue": "Test Tissue",
                                    "validation_status": "replicated_donor_supported"}],
                      ["candidate_id", "tf1_ensembl", "tf1_symbol", "tf2_ensembl", "tf2_symbol",
                       "cell_type_id", "cell_type", "tissue_id", "tissue", "validation_status"])
            meme.write_text("MEME version 5\n\nALPHABET= ACGT\n\n"
                            "MOTIF M1 TFA\nletter-probability matrix: alength= 4 w= 1\n0.25 0.25 0.25 0.25\n\n"
                            "MOTIF M2 TFB\nletter-probability matrix: alength= 4 w= 1\n0.25 0.25 0.25 0.25\n")
            command_plan(Namespace(candidate_summary=candidates, index=cell_index,
                                   crosswalk=crosswalk, meme=meme, out=out,
                                   include_status=["replicated_donor_supported"]))
            scan = out / "scan"
            write_tsv(scan / "0000.hits.tsv.gz",
                      [{"motif_id": motif, "peak_index": i} for motif in ("M1", "M2")
                       for i in range(20)], ["motif_id", "peak_index"])
            write_tsv(scan / "0000.gc.tsv.gz",
                      [{"peak_index": i, "gc": 0.5} for i in range(100)], ["peak_index", "gc"])
            write_json(scan / "0000.done.json", {"shard_index": 0, "num_shards": 1,
                                                   "fimo_p": 1e-4, "max_stored_scores": 5000000,
                                                   "n_distal_peaks": 100, "n_fimo_sites": 40})
            command_compile(Namespace(out=out))
            command_analyze(Namespace(out=out, data_dir=data_dir, num_shards=1, shard_index=0,
                                      min_cells=20, min_open_cells=5, min_open_fraction=0.05,
                                      min_matched=10, min_motif_hits=5, cell_batch_size=10))
            command_aggregate(Namespace(out=out, q=0.05, min_joint=5))
            result = list(read_tsv(out / "aggregate" / "corpus_pair_results.tsv.gz"))
            self.assertEqual(result[0]["corpus_result"], "corpus_coenriched_independent_source")
            self.assertEqual(int(result[0]["joint_fg_hits"]), 20)
            rollup = list(read_tsv(out / "aggregate" / "candidate_summary.tsv"))
            self.assertEqual(rollup[0]["corpus_summary"], "independent_source_support")


if __name__ == "__main__":
    unittest.main()
