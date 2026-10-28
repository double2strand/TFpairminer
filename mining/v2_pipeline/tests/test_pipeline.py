"""Unit tests for the v2 cluster pipeline (no network, no GPU).  Run:
   python -m pytest mining/v2_pipeline/tests -q      (needs pyyaml; numpy/pandas/anndata for the RNA test)
The full end-to-end check is the smoke run (RUNBOOK 'Smoke test')."""
import csv, io, os, subprocess, sys, tempfile, unittest
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "mining")); sys.path.insert(0, str(ROOT / "mining/tf_pair_evidence"))
import pipeline as P  # noqa: E402


def cfg_with(tmp, **over):
    import yaml
    c = yaml.safe_load((HERE / "config.yaml").read_text()); c["out"] = str(tmp)
    for k, v in over.items():
        c[k] = v
    p = Path(tmp) / "cfg.yaml"; p.write_text(yaml.safe_dump(c)); return P.Cfg(p)


class TestConfig(unittest.TestCase):
    def test_env_expansion(self):
        os.environ["TFPM_V2P_OUT"] = "/x/y"
        self.assertEqual(P._expand("${TFPM_V2P_OUT:-/d}/o"), "/x/y/o")
        self.assertEqual(P._expand("${TFPM_NOPE_ZZ:-/d}/o"), "/d/o")

    def test_stage_selection(self):
        self.assertEqual(P.expand_stages("all"), P.STAGES)
        self.assertEqual(P.expand_stages("05"), ["05", "05b"])
        self.assertEqual(P.expand_stages("03-04"), ["03", "03b", "04"])
        self.assertEqual(P.expand_stages("07a,09"), ["07a", "09"])

    def test_resources_every_stage(self):
        with tempfile.TemporaryDirectory() as t:
            c = cfg_with(t)
            for s in P.STAGES:
                r = c.res(s); self.assertTrue(r["cpus"] >= 1 and r["mem"] >= 1, s)
            self.assertIn("--misc", P.wr_res(c, "04")); self.assertNotIn("--misc", P.wr_res(c, "03"))

    def test_secs(self):
        self.assertEqual(P._secs("30m"), 1800); self.assertEqual(P._secs("2d"), 172800); self.assertEqual(P._secs("90"), 90)


class TestManifest(unittest.TestCase):
    def test_samples_tsv(self):
        rows = list(csv.DictReader(open(HERE / "samples.tsv"), delimiter="\t"))
        ids = [r["profile_id"] for r in rows]
        self.assertEqual(len(ids), len(set(ids)))
        ccla = [r for r in rows if r["atlas"].startswith("CCLA")]; corp = [r for r in rows if "Corpus" in r["atlas"]]
        self.assertEqual(len(ccla) + len(corp), len(rows)); self.assertGreaterEqual(len(ccla), 39); self.assertGreaterEqual(len(corp), 77)
        for r in rows:
            self.assertTrue(r["inputs"].startswith("http"), r["profile_id"])
            self.assertTrue(all(s.isdigit() for s in r["sizes"].split(",") if s), r["profile_id"])
        hep = next(r for r in rows if r["profile_id"] == "CCLA_HepG2")
        self.assertIn("ENCSR291GJU", hep["cbp_pretrained"]); self.assertEqual(len(hep["cbp_pretrained_sha256"]), 64)


class TestContextMap(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory(); self.c = cfg_with(self.t.name)

    def tearDown(self):
        self.t.cleanup()

    def ctx(self, profile, celltype="", sample="S"):
        r = P.rna_context(self.c, {"profile": profile, "celltype": celltype, "sample": sample})
        return (r["rna_context"], r["mapping_level"]) if r else None

    def test_lines(self):
        self.assertEqual(self.ctx("CCLA_HepG2"), ("hepatocyte", "proxy"))
        self.assertEqual(self.ctx("CCLA_K562"), ("erythroid_myeloid_progenitor", "proxy"))
        self.assertIsNone(self.ctx("CCLA_A549"))

    def test_corpus(self):
        self.assertEqual(self.ctx("X__Hepatocyte", "Hepatocyte"), ("hepatocyte", "exact"))
        self.assertEqual(self.ctx("X__m", "Macrophage_Liver"), ("macrophage", "proxy"))
        self.assertEqual(self.ctx("X__c", "CD8 T cell"), ("CD8_T", "exact"))
        self.assertIsNone(self.ctx("X__w", "Weird_cluster_7"))

    def test_rules_valid(self):
        for r in P.context_rules(self.c):
            self.assertIn(r["mapping_level"], ("exact", "proxy", "unmapped"))
            if r["mapping_level"] != "unmapped":
                self.assertTrue(r["cl_ids"].startswith("CL:"), r)


class TestSafety(unittest.TestCase):
    def test_fetch_size_mismatch(self):
        with tempfile.TemporaryDirectory() as t:
            src = Path(t) / "a.txt"; src.write_text("hello"); log = Path(t) / "l.log"
            P.fetch_verified(str(src), Path(t) / "d" / "a.txt", log, size="5", lock=Path(t) / "lock.tsv")
            with self.assertRaises(RuntimeError):
                P.fetch_verified(str(src), Path(t) / "e" / "a.txt", log, size="6")

    def test_publish_refuses_v1(self):
        with tempfile.TemporaryDirectory() as t:
            c = cfg_with(t, publish={"dir": "docs", "commit": False})
            u = P.st10(c)[0]
            with self.assertRaises(RuntimeError):
                u.fn({}, Path(t) / "l.log")

    def test_wr_dry_run_chain(self):
        with tempfile.TemporaryDirectory() as t:
            c = cfg_with(t); buf = io.StringIO()
            with redirect_stdout(buf):
                P.submit_wr(c, P.STAGES, dry=True)
            s = buf.getvalue()
            self.assertIn("-i tfpm-v2-exp-00", s); self.assertIn("-d tfpm-v2-exp-10", s)
            self.assertNotIn("sbatch", s)


    def test_farm22_container_and_gpu_flags(self):
        c = P.Cfg(P.HERE / "config.farm22.yaml")
        g = c.gpu_py(); m = c.gpu_py(gpu=False)
        self.assertIn("apptainer exec --nv -B", g); self.assertIn("chrombpnet.sif python", g)
        self.assertIn("module load ISG/apptainer/1.4.0", g); self.assertIn("/lustre/scratch126", g)
        self.assertNotIn("--nv", m)
        self.assertIn("PYTHONPATH=/lustre/scratch126/gengen/projects_v2/ncap_selex/tfpm/cbp_extra:", g)
        a4 = " ".join(P.wr_res(c, "04")); a5 = " ".join(P.wr_res(c, "05")); a3 = " ".join(P.wr_res(c, "03"))
        self.assertIn("--queue gpu-basement", a4); self.assertIn('-gpu "num=1"', a4)
        self.assertIn("--queue gpu ", a5 + " "); self.assertIn("--queue long", a3); self.assertNotIn("-gpu", a3)


    def test_write_tsv_parallel(self):
        import threading
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "m.tsv"; errs = []
            def w():
                try:
                    for _ in range(50):
                        P.write_tsv(p, [{"a": 1}], ["a"])
                except Exception as e:
                    errs.append(e)
            th = [threading.Thread(target=w) for _ in range(8)]; [x.start() for x in th]; [x.join() for x in th]
            self.assertEqual(errs, []); self.assertEqual(p.read_text(), "a\n1\n")


class TestTargetedRnaLocal(unittest.TestCase):
    def test_local_source_fills_missing_columns(self):
        try:
            import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
            import targeted_rna as T
        except Exception as e:  # pragma: no cover
            self.skipTest(f"deps missing: {e}")
        with tempfile.TemporaryDirectory() as t:
            obs = pd.DataFrame({"dataset_id": ["d1"] * 4, "cell_type_ontology_term_id": ["CL:1", "CL:1", "CL:2", "CL:1"],
                                "is_primary_data": [True] * 4}, index=list("abcd"))
            A = ad.AnnData(X=sp.csr_matrix(np.ones((4, 2), dtype=np.float32)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
            p = Path(t) / "x.h5ad"; A.write_h5ad(p)
            class A_: local_h5ad = str(p); census_version = "local"
            src = T.open_source(A_())
            o = src.obs("cell_type_ontology_term_id in ['CL:1'] and is_primary_data == True",
                        ["dataset_id", "tissue_general_ontology_term_id"])
            self.assertEqual(len(o), 3); self.assertTrue((o["tissue_general_ontology_term_id"] == "").all())


if __name__ == "__main__":
    unittest.main()


class TestDisabledRows(unittest.TestCase):
    def test_disabled_rows_skipped(self):
        with tempfile.TemporaryDirectory() as t:
            m = Path(t) / "m.tsv"
            m.write_text("profile_id\tatlas\tcontext\tkind\tinputs\tenabled\tdisabled_reason\n"
                         "A\tCorpus\tc\tcorpus\thttp://x/a.h5ad\t1\t\n"
                         "B\tCorpus\tc\tcorpus\thttp://x/b.h5ad\t0\tHTTP 404\n"
                         "C\tCorpus\tc\tcorpus\thttp://x/c.h5ad\t\t\n")
            cfg = cfg_with(t, manifest=str(m))
            self.assertEqual([r["profile_id"] for r in P.samples(cfg)], ["A", "C"])
            self.assertEqual([r["profile_id"] for r in P.disabled_samples(cfg)], ["B"])
            self.assertNotIn("B", [u.sample for u in P.stage_units(cfg, "01")])

    def test_shipped_manifest_rows_have_enabled(self):
        rows = P.read_tsv(HERE / "samples.tsv")
        self.assertEqual(len(rows), 116)
        for r in rows:
            self.assertIn(r["enabled"], ("0", "1"))
            if r["enabled"] == "0":
                self.assertTrue(r["disabled_reason"])


def test_host_slot_limits_and_curl_resume(tmp_path, monkeypatch):
    import pipeline as P, threading, time as _t
    monkeypatch.setitem(P.FETCH, "slot_dir", str(tmp_path / "slots"))
    monkeypatch.setitem(P.FETCH, "max_parallel_per_host", 2)
    log = tmp_path / "l.log"; cur = [0]; peak = [0]; lk = threading.Lock()
    _sl = _t.sleep; monkeypatch.setattr(P.time, "sleep", lambda s: _sl(0.01))
    def w():
        with P.host_slot("h", log):
            with lk: cur[0] += 1; peak[0] = max(peak[0], cur[0])
            _sl(0.05)
            with lk: cur[0] -= 1
    # flock is per open-file-description, so threads contend like separate processes
    ts = [threading.Thread(target=w) for _ in range(6)]; [t.start() for t in ts]; [t.join() for t in ts]
    assert peak[0] <= 2
    calls = []
    def fake_run(cmd, log):
        calls.append(cmd)
        if len(calls) < 3: raise RuntimeError("command failed (92): curl")
        Path(tmp_path / "f.bin.part").write_text("x")
    monkeypatch.setattr(P, "run", fake_run); monkeypatch.setitem(P.FETCH, "backoff", 0)
    P.curl("https://h/f.bin", tmp_path / "f.bin", log)
    assert len(calls) == 3 and "--http1.1" in calls[0] and "--retry-all-errors" in calls[0] and "-C -" in calls[0]
    assert (tmp_path / "f.bin").read_text() == "x"
