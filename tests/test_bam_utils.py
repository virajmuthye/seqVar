import os
import tempfile
import unittest

from conftest import make_dataset
from seqvar import bam_utils


class BamUtilsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.info = make_dataset(cls.tmp)
        cls.sam = os.path.join(cls.tmp, "alignment.sam")
        cls.contig = cls.info["contig"]

    def test_flagstat(self):
        fs = bam_utils.flagstat(self.sam)
        self.assertEqual(fs.total, self.info["n_reads"])
        self.assertEqual(fs.mapped, self.info["n_reads"])
        self.assertEqual(fs.unmapped, 0)
        self.assertAlmostEqual(fs.mapped_fraction, 1.0)

    def test_coverage_positive(self):
        depths = bam_utils.coverage(self.sam, self.contig, start=1, end=60)
        self.assertEqual(len(depths), 60)
        self.assertTrue(all(d >= 1 for d in depths))
        # Interior positions get tiled coverage from several reads.
        self.assertGreater(max(depths), 1)

    def test_filter_by_mapq(self):
        out = os.path.join(self.tmp, "filtered.sam")
        # All reads have MAPQ 60; a threshold of 61 should keep none.
        kept = bam_utils.filter_records(self.sam, out, min_mapq=61)
        self.assertEqual(kept, 0)
        kept = bam_utils.filter_records(self.sam, out, min_mapq=60)
        self.assertEqual(kept, self.info["n_reads"])

    def test_filter_exclude_flags(self):
        out = os.path.join(self.tmp, "filtered2.sam")
        # Excluding unmapped (0x4) keeps everything since all are mapped.
        kept = bam_utils.filter_records(self.sam, out, exclude_flags=0x4)
        self.assertEqual(kept, self.info["n_reads"])


if __name__ == "__main__":
    unittest.main()
