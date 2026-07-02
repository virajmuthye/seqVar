import os
import tempfile
import unittest

from conftest import make_dataset
from seqvar import call_variants
from seqvar.call_variants import CallParams


class CallVariantsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.info = make_dataset(cls.tmp)
        cls.sam = os.path.join(cls.tmp, "alignment.sam")
        cls.ref = os.path.join(cls.tmp, "ref.fa")
        cls.vcf = call_variants.call_variants(cls.sam, cls.ref)

    def test_calls_both_variants(self):
        # Exactly the SNP and the deletion should be recovered.
        self.assertEqual(len(self.vcf), 2)

    def test_snp_call(self):
        snp_pos, ref_base, alt_base = self.info["snp"]
        snps = [r for r in self.vcf if r.is_snp]
        self.assertEqual(len(snps), 1)
        rec = snps[0]
        self.assertEqual(rec.pos, snp_pos)
        self.assertEqual(rec.ref, ref_base)
        self.assertEqual(rec.alt, [alt_base])
        self.assertGreater(rec.qual, 20)
        self.assertEqual(rec.samples[0]["GT"], "1/1")

    def test_deletion_call(self):
        del_pos, del_ref, del_alt = self.info["deletion"]
        dels = [r for r in self.vcf if r.is_indel]
        self.assertEqual(len(dels), 1)
        rec = dels[0]
        self.assertEqual(rec.pos, del_pos)
        self.assertEqual(rec.ref, del_ref)
        self.assertEqual(rec.alt, [del_alt])

    def test_depth_threshold_suppresses_calls(self):
        # Requiring absurd depth should yield no calls.
        strict = CallParams(min_depth=1000)
        vcf = call_variants.call_variants(self.sam, self.ref, params=strict)
        self.assertEqual(len(vcf), 0)

    def test_binomial_tail_monotonic(self):
        # More alt observations -> smaller error probability -> higher QUAL.
        p = CallParams()
        q_low = call_variants._phred_qual(2, 20, p)
        q_high = call_variants._phred_qual(15, 20, p)
        self.assertGreater(q_high, q_low)


if __name__ == "__main__":
    unittest.main()
