import os
import tempfile
import unittest

from conftest import make_dataset
from seqvar import fastq_qc


class FastqQcTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.info = make_dataset(cls.tmp)
        cls.fastq = os.path.join(cls.tmp, "reads.fastq")

    def test_parse_counts(self):
        records = list(fastq_qc.parse_fastq(self.fastq))
        self.assertEqual(len(records), self.info["n_reads"])
        name, seq, qual = records[0]
        self.assertEqual(len(seq), len(qual))
        self.assertTrue(name.startswith("read"))

    def test_report_metrics(self):
        rep = fastq_qc.qc_fastq(self.fastq)
        self.assertEqual(rep.n_reads, self.info["n_reads"])
        self.assertEqual(rep.min_length, 60)
        self.assertEqual(rep.max_length, 60)
        self.assertAlmostEqual(rep.mean_quality, 40.0, places=5)
        self.assertEqual(rep.n_content, 0.0)
        self.assertEqual(len(rep.per_position_quality), 60)
        self.assertTrue(0.0 < rep.gc_content < 1.0)

    def test_malformed_raises(self):
        bad = os.path.join(self.tmp, "bad.fastq")
        with open(bad, "w") as fh:
            fh.write("@r1\nACGT\n+\nII\n")  # seq/qual length mismatch
        with self.assertRaises(ValueError):
            list(fastq_qc.parse_fastq(bad))

    def test_empty_file(self):
        empty = os.path.join(self.tmp, "empty.fastq")
        open(empty, "w").close()
        rep = fastq_qc.qc_fastq(empty)
        self.assertEqual(rep.n_reads, 0)


if __name__ == "__main__":
    unittest.main()
