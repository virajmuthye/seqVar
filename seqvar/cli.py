"""Unified command-line entry point for the seqvar toolkit.

Dispatches to the individual tool modules so the whole pipeline is reachable
from one command:

    seqvar fastqc   reads.fastq
    seqvar align    ref.fa r1.fq r2.fq -o aln.bam
    seqvar bam      flagstat aln.bam
    seqvar call     aln.sam -r ref.fa -o calls.vcf
    seqvar vcf      stats calls.vcf
    seqvar compare  truth.vcf calls.vcf
"""

from __future__ import annotations

import argparse
import sys

from . import __version__
from . import fastq_qc, align, bam_utils, call_variants, vcf, vcf_compare


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="seqvar",
        description="Sequencer-to-variant toolkit: QC, align, BAM ops, call, VCF.",
    )
    parser.add_argument("--version", action="version",
                        version=f"seqvar {__version__}")
    sub = parser.add_subparsers(dest="tool", required=True)

    sub.add_parser("fastqc", add_help=False, help="FASTQ quality control.")
    sub.add_parser("align", add_help=False, help="Align reads (bwa + samtools).")
    sub.add_parser("bam", add_help=False, help="SAM/BAM stats, filter, coverage.")
    sub.add_parser("call", add_help=False, help="Call variants -> VCF.")
    sub.add_parser("vcf", add_help=False, help="Filter / summarise VCF.")
    sub.add_parser("compare", add_help=False, help="Compare two VCFs.")

    if not argv:
        parser.print_help()
        return 1

    tool, rest = argv[0], argv[1:]
    dispatch = {
        "fastqc": fastq_qc.main,
        "align": align.main,
        "bam": bam_utils.main,
        "call": call_variants.main,
        "vcf": vcf.main,
        "compare": vcf_compare.main,
    }
    if tool in ("-h", "--help"):
        parser.print_help()
        return 0
    if tool == "--version":
        print(f"seqvar {__version__}")
        return 0
    if tool not in dispatch:
        parser.print_help()
        return 1
    return dispatch[tool](rest)


if __name__ == "__main__":
    raise SystemExit(main())
