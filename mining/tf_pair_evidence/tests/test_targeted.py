import unittest
from tf_pair_evidence.targeted_integrate import rna_status


class TestRnaStatus(unittest.TestCase):
    def test_missing_screen_is_unassessed(self):
        self.assertEqual(rna_status(None, [("A", "B")])["status"], "unassessed")

    def test_untested_subpair_makes_triplet_unassessed(self):
        per = {("A", "B"): [(True, .4, .01, "s")], ("A", "C"): [(True, .3, .02, "s")]}
        self.assertEqual(rna_status(per, [("A", "B"), ("A", "C"), ("B", "C")])["status"], "unassessed")

    def test_triplet_needs_all_subpairs(self):
        per = {("A", "B"): [(True, .4, .01, "s")], ("A", "C"): [(True, .3, .02, "s")], ("B", "C"): [(False, .1, .5, "s")]}
        self.assertEqual(rna_status(per, [("A", "B"), ("A", "C"), ("B", "C")])["status"], "screened_no_hit")
        per[("B", "C")] = [(True, .26, .04, "s")]
        r = rna_status(per, [("A", "B"), ("A", "C"), ("B", "C")])
        self.assertEqual((r["status"], r["best_r"]), ("hit", .26))


if __name__ == "__main__":
    unittest.main()
