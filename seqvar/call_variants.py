"""Pileup-based variant caller (pure Python).

Builds a per-position pileup from aligned reads and calls SNVs and short
indels against a reference. The model is intentionally simple and transparent
so you can read exactly how a call is made, rather than being a black box:

1. Walk every alignment and tally, at each reference position, the observed
   bases (for substitutions), plus insertions and deletions.
2. At a position, the alternate allele is the most common non-reference
   observation. It is emitted as a variant when it clears three thresholds:
   minimum depth, minimum supporting-allele count, and minimum allele
   fraction.
3. A Phred-scaled QUAL is derived from a binomial tail probability of seeing
   at least this many alt observations under a background error rate. This is
   a heuristic, not a full genotype-likelihood model like GATK/bcftools, and
   is documented as such.

Genotype is a simple heuristic: allele fraction >= 0.8 is called homozygous
alternate (1/1), otherwise heterozygous (0/1).

Input is SAM or BAM (via ``bam_utils``) plus a reference FASTA. Output is VCF.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass

from . import bam_utils
from .fasta import read_fasta
from .vcf import VcfFile, VcfRecord

# Phred-33 base-quality offset for read qualities.
_QUAL_OFFSET = 33


@dataclass
class CallParams:
    min_depth: int = 4
    min_alt_count: int = 2
    min_alt_fraction: float = 0.2
    min_base_quality: int = 13  # ~ P(error) 0.05; skip low-quality bases
    min_mapq: int = 1
    error_rate: float = 0.01  # background per-base sequencing error for QUAL
    hom_fraction: float = 0.8  # >= this alt fraction -> homozygous alt
    max_qual: float = 200.0


class _PositionTally:
    """Accumulates evidence at one reference position."""

    __slots__ = ("depth", "bases", "insertions", "deletions")

    def __init__(self) -> None:
        self.depth = 0
        self.bases: dict[str, int] = defaultdict(int)
        # inserted sequence -> count (insertion occurs after this position)
        self.insertions: dict[str, int] = defaultdict(int)
        # length of deletion starting here -> count
        self.deletions: dict[int, int] = defaultdict(int)


def build_pileup(
    alignment_path: str,
    reference_name: str,
    params: CallParams,
) -> dict[int, _PositionTally]:
    """Build {ref_pos(1-based): _PositionTally} for one reference sequence."""
    pileup: dict[int, _PositionTally] = defaultdict(_PositionTally)

    with bam_utils.open_alignments(alignment_path) as records:
        for rec in records:
            if rec.is_unmapped or rec.rname != reference_name:
                continue
            if rec.is_secondary or rec.is_supplementary or rec.is_duplicate:
                continue
            if rec.mapq < params.min_mapq:
                continue
            _tally_read(rec, pileup, params)
    return pileup


def _base_quality(rec, qi: int, params: CallParams) -> int:
    """Return the Phred quality of query base qi, or a pass-through if absent."""
    if rec.qual and rec.qual != "*" and qi < len(rec.qual):
        return ord(rec.qual[qi]) - _QUAL_OFFSET
    # No qualities available (e.g. '*'): treat as acceptable.
    return params.min_base_quality


def _tally_read(rec, pileup: dict[int, _PositionTally], params: CallParams) -> None:
    """Walk one alignment's CIGAR and accumulate evidence into the pileup."""
    ops = rec.cigar_ops()
    qi = 0
    rp = rec.pos
    seq = rec.seq
    for length, op in ops:
        if op in ("M", "=", "X"):
            for _ in range(length):
                if qi < len(seq):
                    base = seq[qi].upper()
                    if _base_quality(rec, qi, params) >= params.min_base_quality:
                        tally = pileup[rp]
                        tally.depth += 1
                        tally.bases[base] += 1
                qi += 1
                rp += 1
        elif op == "I":
            # Insertion is anchored to the previous reference position.
            inserted = seq[qi:qi + length].upper()
            anchor = rp - 1
            if anchor >= rec.pos:
                pileup[anchor].insertions[inserted] += 1
            qi += length
        elif op == "D":
            pileup[rp].deletions[length] += 1
            # Deleted reference positions still count toward depth.
            for _ in range(length):
                pileup[rp].depth += 1
                rp += 1
        elif op == "N":
            rp += length
        elif op == "S":
            qi += length
        elif op in ("H", "P"):
            continue


def _phred_qual(alt_count: int, depth: int, params: CallParams) -> float:
    """Phred-scaled confidence that alt evidence is not pure sequencing error.

    Uses the binomial tail P(X >= alt_count | depth, error_rate) and converts
    to a Phred score -10*log10(p). Capped at ``max_qual``.
    """
    p = _binom_tail_ge(alt_count, depth, params.error_rate)
    if p <= 0:
        return params.max_qual
    q = -10.0 * math.log10(p)
    return min(q, params.max_qual)


def _binom_tail_ge(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p), computed directly (n is small here)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    total = 0.0
    # Sum the upper tail; n is read depth, typically tens to low hundreds.
    for i in range(k, n + 1):
        total += math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
    return min(1.0, total)


def call_variants(
    alignment_path: str,
    reference_fasta: str,
    params: CallParams | None = None,
    sample_name: str = "SAMPLE",
) -> VcfFile:
    """Call SNVs and short indels; return a populated :class:`VcfFile`."""
    params = params or CallParams()
    reference = read_fasta(reference_fasta)

    vcf = VcfFile(samples=[sample_name])
    vcf.meta = [
        "##fileformat=VCFv4.2",
        "##source=seqvar.call_variants",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
        '##INFO=<ID=AC,Number=1,Type=Integer,Description="Alt allele count">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Alt allele fraction">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="SNP, INS or DEL">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Ref,Alt depths">',
    ]

    for ref_name, ref_seq in reference.items():
        pileup = build_pileup(alignment_path, ref_name, params)
        for pos in sorted(pileup):
            rec = _call_position(ref_name, ref_seq, pos, pileup[pos], params)
            if rec is not None:
                vcf.records.append(rec)
    return vcf


def _call_position(ref_name, ref_seq, pos, tally, params):
    """Decide whether a variant exists at one position; return a record or None."""
    depth = tally.depth
    if depth < params.min_depth:
        return None

    ref_base = ref_seq[pos - 1].upper() if 0 < pos <= len(ref_seq) else "N"

    # Assemble candidate alt alleles. Each candidate is a dict so that indels
    # can carry their own VCF POS (deletions anchor on the *preceding* base,
    # per the VCF spec) independent of the pileup position.
    candidates: list[dict] = []

    # SNV candidate: most common base that differs from the reference.
    for base, count in tally.bases.items():
        if base != ref_base and base in "ACGT":
            candidates.append(
                {"pos": pos, "ref": ref_base, "alt": base,
                 "count": count, "type": "SNP"}
            )

    # Insertion candidate: anchored on this base (ALT = base + inserted seq).
    for inserted, count in tally.insertions.items():
        candidates.append(
            {"pos": pos, "ref": ref_base, "alt": ref_base + inserted,
             "count": count, "type": "INS"}
        )

    # Deletion candidate: VCF anchors on the base BEFORE the deleted region.
    # REF = anchor + deleted bases; ALT = anchor. Requires a preceding base.
    for dlen, count in tally.deletions.items():
        end = pos - 1 + dlen  # 0-based exclusive end of deleted region
        if pos >= 2 and end <= len(ref_seq):
            anchor = ref_seq[pos - 2].upper()
            deleted_ref = anchor + ref_seq[pos - 1:end].upper()
            candidates.append(
                {"pos": pos - 1, "ref": deleted_ref, "alt": anchor,
                 "count": count, "type": "DEL"}
            )

    if not candidates:
        return None

    best = max(candidates, key=lambda c: c["count"])
    vcf_pos = best["pos"]
    ref_allele = best["ref"]
    alt_allele = best["alt"]
    alt_count = best["count"]
    vtype = best["type"]

    if alt_count < params.min_alt_count:
        return None
    af = alt_count / depth
    if af < params.min_alt_fraction:
        return None

    qual = _phred_qual(alt_count, depth, params)
    gt = "1/1" if af >= params.hom_fraction else "0/1"
    ref_support = tally.bases.get(ref_base, 0)

    info = {
        "DP": depth,
        "AC": alt_count,
        "AF": round(af, 4),
        "TYPE": vtype,
    }
    sample = {
        "GT": gt,
        "DP": str(depth),
        "AD": f"{ref_support},{alt_count}",
    }
    return VcfRecord(
        chrom=ref_name,
        pos=vcf_pos,
        id=".",
        ref=ref_allele,
        alt=[alt_allele],
        qual=round(qual, 1),
        filter=["PASS"] if qual >= 20 else ["LowQual"],
        info=info,
        format_keys=["GT", "DP", "AD"],
        samples=[sample],
    )


def _build_parser():
    p = argparse.ArgumentParser(
        prog="seqvar-call",
        description="Call SNVs and short indels from a SAM/BAM against a reference.",
    )
    p.add_argument("alignment", help="Input SAM or BAM (coordinate order not required).")
    p.add_argument("-r", "--reference", required=True, help="Reference FASTA.")
    p.add_argument("-o", "--output", required=True, help="Output VCF.")
    p.add_argument("-s", "--sample", default="SAMPLE", help="Sample name.")
    p.add_argument("--min-depth", type=int, default=4)
    p.add_argument("--min-alt-count", type=int, default=2)
    p.add_argument("--min-alt-fraction", type=float, default=0.2)
    p.add_argument("--min-base-quality", type=int, default=13)
    p.add_argument("--min-mapq", type=int, default=1)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    params = CallParams(
        min_depth=args.min_depth,
        min_alt_count=args.min_alt_count,
        min_alt_fraction=args.min_alt_fraction,
        min_base_quality=args.min_base_quality,
        min_mapq=args.min_mapq,
    )
    vcf = call_variants(
        args.alignment, args.reference, params=params, sample_name=args.sample
    )
    vcf.write(args.output)
    print(f"Called {len(vcf)} variants -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
