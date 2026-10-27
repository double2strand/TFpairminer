import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from tf_pair_evidence.integrate import integrate, read_tsv


CTX = {"cell_type_id": "CL:0001", "tissue_id": "UBERON:0001", "disease_id": "PATO:0000461"}
GENES = {"TFA": "ENSG00000000001", "TFB": "ENSG00000000002", "TFC": "ENSG00000000003",
         "TFD": "ENSG00000000004", "TFE": "ENSG00000000005"}


def tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = (gzip.open(path, "wt", newline="", encoding="utf-8") if path.suffix == ".gz"
              else path.open("w", newline="", encoding="utf-8"))
    with handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class IntegrationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.register = self.root / "gene_register.tsv"
        tsv(self.register, ("Symbol", "Ensembl_gene_ID"),
            [{"Symbol": s, "Ensembl_gene_ID": e} for s, e in GENES.items()])
        digest = hashlib.sha256(self.register.read_bytes()).hexdigest()
        self.rna1 = self.make_rna("rna1", "dataset-one", digest,
                                  [("TFA", "TFB"), ("TFA", "TFC"), ("TFB", "TFC"), ("TFB", "TFD")])
        self.rna2 = self.make_rna("rna2", "dataset-two", digest,
                                  [("TFA", "TFB"), ("TFA", "TFD"), ("TFB", "TFC"), ("TFB", "TFD")])
        cbp = self.root / "pairs_cbp.tsv"
        cbp_rows = []
        for a, b, support in (("TFA", "TFB", "yes"), ("TFA", "TFC", "yes"),
                               ("TFA", "TFD", "yes"), ("TFA", "TFE", "yes"),
                               ("TFB", "TFC", "no")):
            cbp_rows.append({"tfs": f"{a}+{b}", "k": "2", "status": "tested",
                             "contribution_supported": support, "synergy_q": "0.02",
                             "synergy": "0.3", "cooperative": "synergistic"})
        tsv(cbp, tuple(cbp_rows[0]), cbp_rows)
        self.atac_manifest = self.root / "atac_manifest.tsv"
        tsv(self.atac_manifest, ("atac_dataset", "atac_context", "pairs_cbp_path"),
            [{"atac_dataset": "atlas-one", "atac_context": "cell-A", "pairs_cbp_path": str(cbp)}])
        self.context_map = self.root / "context_map.tsv"
        tsv(self.context_map, ("atac_dataset", "atac_context", *CTX, "mapping_level"),
            [{"atac_dataset": "atlas-one", "atac_context": "cell-A", **CTX, "mapping_level": "exact"}])
        self.out = self.root / "results"

    def tearDown(self):
        self.temp.cleanup()

    def make_rna(self, name, dataset, digest, hits):
        root = self.root / name
        root.mkdir()
        (root / "manifest.json").write_text(json.dumps({"resolved_census_version": "2026-01-01",
            "register_sha256": digest, "dataset_selection": "all" if name == "rna1" else "explicit_ids",
            "excluded_collection_names": ["Tabula Sapiens"] if name == "rna1" else []}), encoding="utf-8")
        tsv(root / "datasets.tsv", ("dataset_id", "collection_id", "collection_name", "dataset_title"),
            [{"dataset_id": dataset, "collection_id": name,
              "collection_name": "Other atlas" if name == "rna1" else "Tabula Sapiens",
              "dataset_title": "Other cells" if name == "rna1" else "Tabula Sapiens - All Cells"}])
        (root / "aggregate").mkdir()
        (root / "aggregate" / "summary.json").write_text(json.dumps({"complete": True,
            "atlas_wide_except_heldout_collections": name == "rna1",
            "run_config": {"max_strata": 0, "top_k": 0, "primary_stat": "metacell",
                           "max_q": 0.05, "min_cells": 100, "min_detection": 0.05}}), encoding="utf-8")
        profiles = [{**CTX, "assay_id": "EFO:1", "ensembl_id": gene,
                     "n_cells": "200", "detection_fraction": "0.5"} for gene in GENES.values()]
        tsv(root / "aggregate" / "cell_type_profiles.tsv", tuple(profiles[0]), profiles)
        tsv(root / "datasets" / dataset / "profiles.tsv.gz", tuple(profiles[0]), profiles)
        candidates = [{**CTX, "assay_id": "EFO:1", "tf1_ensembl": GENES[a],
                       "tf2_ensembl": GENES[b], "primary_stat": "metacell",
                       "best_primary_r": "0.4", "min_q_within_stratum": "0.01",
                       "mean_co_detection_lift": "1.2", "co_detected_cells": "50"} for a, b in hits]
        fields = (*CTX, "assay_id", "tf1_ensembl", "tf2_ensembl", "primary_stat",
                  "best_primary_r", "min_q_within_stratum", "mean_co_detection_lift", "co_detected_cells")
        tsv(root / "aggregate" / "cell_type_pair_candidates.tsv", fields, candidates)
        return root

    def run_pipeline(self):
        integrate(SimpleNamespace(atac_manifest=self.atac_manifest, context_map=self.context_map,
                                  register=self.register, rna1=self.rna1, rna2=self.rna2, out=self.out))

    def test_four_views_and_unassessed(self):
        self.run_pipeline()
        expected = {"01_atac_only.tsv": ("TFA", "TFE"),
                    "02_atac_plus_rna1.tsv": ("TFA", "TFC"),
                    "03_atac_plus_rna2.tsv": ("TFA", "TFD"),
                    "04_rna_replicated_cbp_no_support.tsv": ("TFB", "TFC")}
        for file, pair in expected.items():
            pairs = {(r["tf1_symbol"], r["tf2_symbol"]) for r in read_tsv(self.out / file)}
            self.assertIn(pair, pairs)
        both = list(read_tsv(self.out / "06_atac_plus_both_rna.tsv"))
        self.assertEqual([(r["tf1_symbol"], r["tf2_symbol"]) for r in both], [("TFA", "TFB")])
        modules = list(read_tsv(self.out / "modules.tsv"))
        self.assertTrue(any(r["view"] == "atac_rna1" and r["n_tfs"] == "3"
                            and r["tf_symbols"] == "TFA;TFB;TFC" for r in modules))
        unassessed = list(read_tsv(self.out / "05_unassessed.tsv"))
        self.assertTrue(any((r["tf1_symbol"], r["tf2_symbol"]) == ("TFB", "TFD") for r in unassessed))

    def test_incomplete_screen_fails_closed(self):
        path = self.rna2 / "aggregate" / "summary.json"
        summary = json.loads(path.read_text())
        summary["complete"] = False
        path.write_text(json.dumps(summary))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            self.run_pipeline()

    def test_proxy_context_excluded(self):
        tsv(self.context_map, ("atac_dataset", "atac_context", *CTX, "mapping_level"),
            [{"atac_dataset": "atlas-one", "atac_context": "cell-A", **CTX, "mapping_level": "proxy"}])
        self.run_pipeline()
        for file in ("01_atac_only.tsv", "02_atac_plus_rna1.tsv", "03_atac_plus_rna2.tsv",
                     "04_rna_replicated_cbp_no_support.tsv"):
            self.assertEqual(list(read_tsv(self.out / file)), [])

    def test_missing_same_dataset_coverage_is_not_atac_only(self):
        # RNA2 sees TFA but not TFE in the only source dataset. An aggregate-
        # level profile could conceal this if other datasets supplied TFE.
        path = self.rna2 / "datasets" / "dataset-two" / "profiles.tsv.gz"
        rows = [r for r in read_tsv(path) if r["ensembl_id"] != GENES["TFE"]]
        tsv(path, tuple(rows[0]), rows)
        self.run_pipeline()
        self.assertFalse(any(r["tf2_symbol"] == "TFE" for r in read_tsv(self.out / "01_atac_only.tsv")))
        row = next(r for r in read_tsv(self.out / "master.tsv") if r["tf2_symbol"] == "TFE")
        self.assertEqual(row["rna2_status"], "unassessed")

    def test_cell_raw_hit_is_not_significant(self):
        path = self.rna1 / "aggregate" / "cell_type_pair_candidates.tsv"
        rows = list(read_tsv(path))
        target = next(r for r in rows if r["tf1_ensembl"] == GENES["TFA"] and r["tf2_ensembl"] == GENES["TFC"])
        target["primary_stat"] = "cell_raw"
        target["min_q_within_stratum"] = ""
        tsv(path, tuple(rows[0]), rows)
        self.run_pipeline()
        self.assertFalse(any(r["tf2_symbol"] == "TFC" for r in read_tsv(self.out / "02_atac_plus_rna1.tsv")))


if __name__ == "__main__":
    unittest.main()
