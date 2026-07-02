"""FASTQ quality-control metrics.

Pure Python, no third-party dependencies. Reads plain or gzipped FASTQ and
computes the summary statistics you typically want before alignment: read
counts, length distribution, per-position mean Phred quality, GC content and
the fraction of ambiguous (N) bases.

Phred encoding is assumed to be Sanger / Illumina 1.8+ (offset 33), which has
been the standard since 2011.

Example
-------
    from seqvar.fastq_qc import qc_fastq
    report = qc_fastq("reads.fastq")
    print(report.mean_quality)

CLI
---
    python -m seqvar.fastq_qc reads.fastq
    python -m seqvar.fastq_qc reads.fastq --json report.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from dataclasses import dataclass, field, asdict
from typing import Iterator, TextIO

PHRED_OFFSET = 33


def _open(path: str) -> TextIO:
    """Open a FASTQ file transparently, whether gzipped or plain text."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_fastq(path: str) -> Iterator[tuple[str, str, str]]:
    """Yield (name, sequence, quality) tuples from a FASTQ file.

    Raises ValueError on malformed records (wrong number of lines, mismatched
    sequence/quality length, or a missing '@'/'+' marker), which is more useful
    than silently producing garbage downstream.
    """
    handle = _open(path)
    try:
        while True:
            header = handle.readline()
            if not header:
                break
            header = header.rstrip("\n")
            seq = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            qual = handle.readline().rstrip("\n")
            if not header.startswith("@"):
                raise ValueError(f"Expected '@' header, got: {header!r}")
            if not plus.startswith("+"):
                raise ValueError(f"Expected '+' separator, got: {plus!r}")
            if len(seq) != len(qual):
                raise ValueError(
                    f"Sequence/quality length mismatch in record {header!r}: "
                    f"{len(seq)} vs {len(qual)}"
                )
            yield header[1:].split()[0], seq, qual
    finally:
        if handle is not sys.stdin:
            handle.close()


@dataclass
class FastqReport:
    """Summary QC metrics for a single FASTQ file."""

    n_reads: int = 0
    total_bases: int = 0
    min_length: int = 0
    max_length: int = 0
    mean_length: float = 0.0
    gc_content: float = 0.0
    n_content: float = 0.0
    mean_quality: float = 0.0
    # Per-position mean quality, index 0 == first base of the read.
    per_position_quality: list[float] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        lines = [
            f"reads            : {self.n_reads:,}",
            f"total bases      : {self.total_bases:,}",
            f"read length      : {self.min_length}-{self.max_length} "
            f"(mean {self.mean_length:.1f})",
            f"GC content       : {self.gc_content * 100:.2f}%",
            f"N content        : {self.n_content * 100:.4f}%",
            f"mean quality     : Q{self.mean_quality:.1f}",
        ]
        return "\n".join(lines)


def qc_fastq(path: str) -> FastqReport:
    """Compute QC metrics for a FASTQ file in a single streaming pass."""
    n_reads = 0
    total_bases = 0
    gc = 0
    n_bases = 0
    qual_sum = 0
    min_len = None
    max_len = 0

    # Accumulate per-position quality sums; grows as needed for variable lengths.
    pos_qual_sum: list[int] = []
    pos_count: list[int] = []

    for _name, seq, qual in parse_fastq(path):
        n_reads += 1
        length = len(seq)
        total_bases += length
        min_len = length if min_len is None else min(min_len, length)
        max_len = max(max_len, length)

        upper = seq.upper()
        gc += upper.count("G") + upper.count("C")
        n_bases += upper.count("N")

        if length > len(pos_qual_sum):
            pos_qual_sum.extend([0] * (length - len(pos_qual_sum)))
            pos_count.extend([0] * (length - len(pos_count)))

        for i, ch in enumerate(qual):
            q = ord(ch) - PHRED_OFFSET
            qual_sum += q
            pos_qual_sum[i] += q
            pos_count[i] += 1

    if n_reads == 0:
        return FastqReport()

    per_position = [
        (pos_qual_sum[i] / pos_count[i]) if pos_count[i] else 0.0
        for i in range(len(pos_qual_sum))
    ]

    return FastqReport(
        n_reads=n_reads,
        total_bases=total_bases,
        min_length=min_len or 0,
        max_length=max_len,
        mean_length=total_bases / n_reads,
        gc_content=gc / total_bases if total_bases else 0.0,
        n_content=n_bases / total_bases if total_bases else 0.0,
        mean_quality=qual_sum / total_bases if total_bases else 0.0,
        per_position_quality=per_position,
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="seqvar-fastqc",
        description="Compute quality-control metrics for a FASTQ file.",
    )
    p.add_argument("fastq", help="Input FASTQ file (plain or .gz; '-' for stdin).")
    p.add_argument("--json", metavar="PATH", help="Write full report as JSON to PATH.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    report = qc_fastq(args.fastq)
    print(report.summary())
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report.to_dict(), fh, indent=2)
        print(f"\nFull report written to {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
