"""Compare two VCF call sets.

Given a "truth" VCF and a "test" VCF, compute concordance: true positives
(shared), false positives (test-only) and false negatives (truth-only), plus
precision, recall and F1. This is the everyday task when benchmarking a caller
or checking two pipelines against each other.

Matching is by normalised (CHROM, POS, REF, ALT). Records are left-normalised
so that, for example, indels written with different but equivalent
representations still line up when they share the same trimmed coordinates.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field

from .vcf import VcfFile, VcfRecord, read_vcf


def normalize_key(rec: VcfRecord) -> tuple[str, int, str, str]:
    """Left-normalise a biallelic record to a canonical (chrom,pos,ref,alt) key.

    Trims common trailing then leading bases (parsimony), which collapses
    equivalent indel spellings such as REF=CT/ALT=C and REF=CTT/ALT=CT onto a
    consistent representation. Only the first ALT is considered.
    """
    chrom = rec.chrom
    pos = rec.pos
    ref = rec.ref.upper()
    alt = (rec.alt[0].upper() if rec.alt else ".")

    # Trim shared suffix (keep at least one base on each allele).
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    # Trim shared prefix, advancing the position accordingly.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return (chrom, pos, ref, alt)


@dataclass
class ComparisonResult:
    true_positives: list[VcfRecord] = field(default_factory=list)
    false_positives: list[VcfRecord] = field(default_factory=list)
    false_negatives: list[VcfRecord] = field(default_factory=list)

    @property
    def tp(self) -> int:
        return len(self.true_positives)

    @property
    def fp(self) -> int:
        return len(self.false_positives)

    @property
    def fn(self) -> int:
        return len(self.false_negatives)

    @property
    def precision(self) -> float:
        denom = self.tp + self.fp
        return self.tp / denom if denom else 0.0

    @property
    def recall(self) -> float:
        denom = self.tp + self.fn
        return self.tp / denom if denom else 0.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if (p + r) else 0.0

    def summary(self) -> str:
        return (
            f"true positives  (shared)     : {self.tp}\n"
            f"false positives (test-only)  : {self.fp}\n"
            f"false negatives (truth-only) : {self.fn}\n"
            f"precision                    : {self.precision:.4f}\n"
            f"recall (sensitivity)         : {self.recall:.4f}\n"
            f"F1 score                     : {self.f1:.4f}"
        )


def compare(
    truth: VcfFile,
    test: VcfFile,
    snps_only: bool = False,
    indels_only: bool = False,
) -> ComparisonResult:
    """Compare two loaded VCFs and return a :class:`ComparisonResult`."""

    def keep(rec: VcfRecord) -> bool:
        if snps_only and not rec.is_snp:
            return False
        if indels_only and not rec.is_indel:
            return False
        return True

    truth_map = {normalize_key(r): r for r in truth if keep(r)}
    test_map = {normalize_key(r): r for r in test if keep(r)}

    truth_keys = set(truth_map)
    test_keys = set(test_map)

    result = ComparisonResult()
    for k in test_keys & truth_keys:
        result.true_positives.append(test_map[k])
    for k in test_keys - truth_keys:
        result.false_positives.append(test_map[k])
    for k in truth_keys - test_keys:
        result.false_negatives.append(truth_map[k])

    # Deterministic ordering makes output and tests reproducible.
    for lst in (result.true_positives, result.false_positives,
                result.false_negatives):
        lst.sort(key=lambda r: (r.chrom, r.pos, r.ref, r.alt[0] if r.alt else ""))
    return result


def compare_files(
    truth_path: str,
    test_path: str,
    snps_only: bool = False,
    indels_only: bool = False,
) -> ComparisonResult:
    return compare(
        read_vcf(truth_path), read_vcf(test_path),
        snps_only=snps_only, indels_only=indels_only,
    )


def _build_parser():
    p = argparse.ArgumentParser(
        prog="seqvar-vcfcompare",
        description="Compare a test VCF against a truth VCF (precision/recall).",
    )
    p.add_argument("truth", help="Truth / gold-standard VCF.")
    p.add_argument("test", help="Test VCF to evaluate.")
    p.add_argument("--snps-only", action="store_true")
    p.add_argument("--indels-only", action="store_true")
    p.add_argument("--fp-out", help="Optional: write false positives to this VCF.")
    p.add_argument("--fn-out", help="Optional: write false negatives to this VCF.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    truth = read_vcf(args.truth)
    test = read_vcf(args.test)
    result = compare(
        truth, test, snps_only=args.snps_only, indels_only=args.indels_only
    )
    print(result.summary())

    if args.fp_out:
        out = VcfFile(meta=test.meta[:], samples=test.samples[:])
        out.records = result.false_positives
        out.write(args.fp_out)
    if args.fn_out:
        out = VcfFile(meta=truth.meta[:], samples=truth.samples[:])
        out.records = result.false_negatives
        out.write(args.fn_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
