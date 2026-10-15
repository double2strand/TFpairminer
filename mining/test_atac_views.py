import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import atac_views as av


class AtacViewsTest(unittest.TestCase):
    def test_similarity_identity_and_revcomp(self):
        m = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]], float)
        self.assertAlmostEqual(av.similarity(m, m), 1.0)
        self.assertAlmostEqual(av.similarity(m, m[::-1, ::-1]), 1.0)
        other = np.array([[0.25] * 4] * 6) + np.eye(6, 4) * 0.1
        self.assertLess(av.similarity(m, other), 0.8)

    def test_wilson(self):
        p, lo, hi = av.wilson(50, 100)
        self.assertAlmostEqual(p, 0.5); self.assertTrue(lo < 0.5 < hi)
        self.assertEqual(av.wilson(0, 0), [None, None, None])


if __name__ == "__main__":
    unittest.main()
