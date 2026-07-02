import os
import tempfile
import unittest

from conftest import make_dataset
from seqvar import call_variants, vcf_compare
from seqvar.vcf import read_vcf, VcfRecord


class VcfCompareTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.info = make_dataset(cls.tmp)
        cls.truth = read_vcf(os.path.join(cls.tmp, "truth.vcf"))
        calls = call_variants.call_variants(
            os.path.join(cls.tmp, "alignment.sam"),
            os.path.join(cls.tmp, "ref.fa"),
        )
        cls.calls = calls

    def test_perfect_recovery(self):
        res = vcf_compare.compare(self.truth, self.calls)
        self.assertEqual(res.tp, 2)
        self.assertEqual(res.fp, 0)
        self.assertEqual(res.fn, 0)
        self.assertAlmostEqual(res.precision, 1.0)
        self.assertAlmostEqual(res.recall, 1.0)
        self.assertAlmostEqual(res.f1, 1.0)

    def test_normalize_key_trims_indel(self):
        # Two equivalent spellings of the same 1 bp deletion normalise equal.
        a = VcfRecord("chr1", 100, ".", "CT", ["C"], None, [], {})
        b = VcfRecord("chr1", 100, ".", "CTT", ["CT"], None, [], {})
        self.assertEqual(
            vcf_compare.normalize_key(a), vcf_compare.normalize_key(b)
        )

    def test_false_positive_and_negative(self):
        # Truth has a variant the test misses; test has one truth lacks.
        truth = read_vcf(os.path.join(self.tmp, "truth.vcf"))
        subset = read_vcf(os.path.join(self.tmp, "truth.vcf"))
        # Remove one truth record from the "test" set -> one FN.
        subset.records = subset.records[:1]
        # Add a spurious record -> one FP.
        subset.records.append(
            VcfRecord("chr_test", 10, ".", "A", ["C"], 50.0, ["PASS"], {})
        )
        res = vcf_compare.compare(truth, subset)
        self.assertEqual(res.tp, 1)
        self.assertEqual(res.fp, 1)
        self.assertEqual(res.fn, 1)

    def test_snps_only_filter(self):
        res = vcf_compare.compare(self.truth, self.calls, snps_only=True)
        self.assertEqual(res.tp, 1)  # only the SNP is compared


if __name__ == "__main__":
    unittest.main()
