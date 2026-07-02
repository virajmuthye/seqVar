import unittest

from seqvar import sam
from seqvar.sam import parse_cigar, parse_sam_line


class CigarTest(unittest.TestCase):
    def test_parse_cigar(self):
        self.assertEqual(parse_cigar("60M"), [(60, "M")])
        self.assertEqual(parse_cigar("10M2D48M"),
                         [(10, "M"), (2, "D"), (48, "M")])
        self.assertEqual(parse_cigar("*"), [])

    def test_malformed_cigar(self):
        with self.assertRaises(ValueError):
            parse_cigar("10Z")

    def test_reference_length(self):
        rec = parse_sam_line(
            "r\t0\tchr\t100\t60\t10M2D48M\t*\t0\t0\t"
            + "A" * 58 + "\t" + "I" * 58
        )
        # 10M + 2D + 48M consume 60 reference bases.
        self.assertEqual(rec.reference_length, 60)
        self.assertEqual(rec.reference_end, 160)

    def test_flags(self):
        rec = parse_sam_line(
            "r\t147\tchr\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII"
        )
        self.assertTrue(rec.flag & sam.FLAG_PAIRED)
        self.assertTrue(rec.is_reverse)
        self.assertTrue(rec.is_proper_pair)
        self.assertTrue(rec.is_mapped)
        self.assertFalse(rec.is_secondary)

    def test_aligned_pairs_with_deletion(self):
        rec = parse_sam_line(
            "r\t0\tchr\t10\t60\t2M2D2M\t*\t0\t0\tACGT\tIIII"
        )
        pairs = list(rec.aligned_pairs())
        # 2 matches, 2 deletions (query None), 2 matches
        self.assertEqual(pairs[0], (0, 10))
        self.assertEqual(pairs[1], (1, 11))
        self.assertEqual(pairs[2], (None, 12))
        self.assertEqual(pairs[3], (None, 13))
        self.assertEqual(pairs[4], (2, 14))
        self.assertEqual(pairs[5], (3, 15))

    def test_short_line_raises(self):
        with self.assertRaises(ValueError):
            parse_sam_line("too\tfew\tfields")


if __name__ == "__main__":
    unittest.main()
