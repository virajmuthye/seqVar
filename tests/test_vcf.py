import os
import tempfile
import unittest

from seqvar import vcf


SAMPLE_VCF = """\
##fileformat=VCFv4.2
##source=test
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP\t0/1:30
chr1\t200\t.\tC\tT\t10\tLowQual\tDP=5\tGT:DP\t0/1:5
chr1\t300\t.\tCT\tC\t60\tPASS\tDP=40\tGT:DP\t1/1:40
chr1\t400\t.\tA\tT\t80\tPASS\tDP=50\tGT:DP\t1/1:50
"""


class VcfTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.path = os.path.join(cls.tmp, "sample.vcf")
        with open(cls.path, "w") as fh:
            fh.write(SAMPLE_VCF)
        cls.vcf = vcf.read_vcf(cls.path)

    def test_parse(self):
        self.assertEqual(len(self.vcf), 4)
        self.assertEqual(self.vcf.samples, ["S1"])
        rec = self.vcf.records[0]
        self.assertEqual(rec.chrom, "chr1")
        self.assertEqual(rec.pos, 100)
        self.assertEqual(rec.ref, "A")
        self.assertEqual(rec.alt, ["G"])
        self.assertEqual(rec.info["DP"], "30")
        self.assertEqual(rec.samples[0]["GT"], "0/1")

    def test_snp_indel_classification(self):
        self.assertTrue(self.vcf.records[0].is_snp)
        self.assertFalse(self.vcf.records[0].is_indel)
        self.assertTrue(self.vcf.records[2].is_indel)  # CT -> C
        self.assertFalse(self.vcf.records[2].is_snp)

    def test_roundtrip(self):
        out = os.path.join(self.tmp, "roundtrip.vcf")
        self.vcf.write(out)
        again = vcf.read_vcf(out)
        self.assertEqual(len(again), len(self.vcf))
        for a, b in zip(again.records, self.vcf.records):
            self.assertEqual(a.to_line(), b.to_line())

    def test_filter_min_qual(self):
        out = vcf.filter_vcf(self.vcf, min_qual=20)
        # Drops the QUAL=10 record.
        self.assertEqual(len(out), 3)

    def test_filter_pass_only(self):
        out = vcf.filter_vcf(self.vcf, pass_only=True)
        self.assertEqual(len(out), 3)

    def test_filter_snps_only(self):
        out = vcf.filter_vcf(self.vcf, snps_only=True)
        self.assertEqual(len(out), 3)  # three SNPs, one indel excluded

    def test_stats_ts_tv(self):
        st = vcf.compute_stats(self.vcf)
        self.assertEqual(st.snps, 3)
        self.assertEqual(st.indels, 1)
        # A>G and C>T are transitions; A>T is a transversion.
        self.assertEqual(st.transitions, 2)
        self.assertEqual(st.transversions, 1)
        self.assertAlmostEqual(st.ts_tv, 2.0)


if __name__ == "__main__":
    unittest.main()
