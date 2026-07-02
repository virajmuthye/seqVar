"""Lightweight SAM parsing.

A pure-Python reader for the SAM text format (SAMv1 spec). This is deliberately
small: it covers the columns and CIGAR/flag logic needed to build pileups and
compute alignment statistics, without pulling in htslib.

For binary BAM input, see ``bam_utils``, which will use pysam when it is
installed and otherwise shell out to ``samtools view`` to stream SAM text
through this parser.

References
----------
SAM/BAM format specification: https://samtools.github.io/hts-specs/SAMv1.pdf
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterator, TextIO

# SAM FLAG bits (subset that we actually reason about).
FLAG_PAIRED = 0x1
FLAG_PROPER_PAIR = 0x2
FLAG_UNMAPPED = 0x4
FLAG_MATE_UNMAPPED = 0x8
FLAG_REVERSE = 0x10
FLAG_MATE_REVERSE = 0x20
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUPLICATE = 0x400
FLAG_SUPPLEMENTARY = 0x800

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# CIGAR operations that consume reference bases / query bases.
_CONSUMES_REF = set("MDN=X")
_CONSUMES_QUERY = set("MIS=X")


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """Return a list of (length, operation) tuples from a CIGAR string.

    An unavailable CIGAR ('*') yields an empty list.
    """
    if cigar == "*" or not cigar:
        return []
    ops = _CIGAR_RE.findall(cigar)
    # Sanity check: the regex must consume the whole string.
    if sum(len(n) + 1 for n, _ in ops) != len(cigar):
        raise ValueError(f"Malformed CIGAR string: {cigar!r}")
    return [(int(n), op) for n, op in ops]


@dataclass
class SamRecord:
    """A single SAM alignment record (11 mandatory fields + optional tags)."""

    qname: str
    flag: int
    rname: str
    pos: int  # 1-based leftmost mapping position, per the SAM spec
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: str
    qual: str
    tags: dict[str, str]

    # --- flag helpers -----------------------------------------------------
    @property
    def is_unmapped(self) -> bool:
        return bool(self.flag & FLAG_UNMAPPED)

    @property
    def is_mapped(self) -> bool:
        return not self.is_unmapped

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & FLAG_REVERSE)

    @property
    def is_secondary(self) -> bool:
        return bool(self.flag & FLAG_SECONDARY)

    @property
    def is_supplementary(self) -> bool:
        return bool(self.flag & FLAG_SUPPLEMENTARY)

    @property
    def is_duplicate(self) -> bool:
        return bool(self.flag & FLAG_DUPLICATE)

    @property
    def is_qcfail(self) -> bool:
        return bool(self.flag & FLAG_QCFAIL)

    @property
    def is_proper_pair(self) -> bool:
        return bool(self.flag & FLAG_PROPER_PAIR)

    @property
    def is_primary(self) -> bool:
        return not (self.is_secondary or self.is_supplementary)

    def cigar_ops(self) -> list[tuple[int, str]]:
        return parse_cigar(self.cigar)

    @property
    def reference_length(self) -> int:
        """Number of reference bases the alignment spans (M/D/N/=/X)."""
        return sum(n for n, op in self.cigar_ops() if op in _CONSUMES_REF)

    @property
    def reference_end(self) -> int:
        """1-based exclusive end coordinate on the reference."""
        return self.pos + self.reference_length

    def aligned_pairs(self) -> Iterator[tuple[int | None, int | None]]:
        """Yield (query_index, ref_pos) pairs walking the CIGAR.

        query_index is 0-based into ``seq``; ref_pos is 1-based on the
        reference. Either element is None for an insertion (no ref) or a
        deletion/skip (no query). Soft clips advance the query only; hard
        clips are skipped entirely.
        """
        qi = 0
        rp = self.pos
        for length, op in self.cigar_ops():
            if op in ("M", "=", "X"):
                for _ in range(length):
                    yield qi, rp
                    qi += 1
                    rp += 1
            elif op == "I":
                for _ in range(length):
                    yield qi, None
                    qi += 1
            elif op in ("D", "N"):
                for _ in range(length):
                    yield None, rp
                    rp += 1
            elif op == "S":
                qi += length
            elif op in ("H", "P"):
                continue


def parse_sam_line(line: str) -> SamRecord:
    """Parse a single non-header SAM line into a :class:`SamRecord`."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        raise ValueError(f"SAM line has {len(fields)} fields, expected >= 11")
    tags: dict[str, str] = {}
    for tag in fields[11:]:
        parts = tag.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return SamRecord(
        qname=fields[0],
        flag=int(fields[1]),
        rname=fields[2],
        pos=int(fields[3]),
        mapq=int(fields[4]),
        cigar=fields[5],
        rnext=fields[6],
        pnext=int(fields[7]),
        tlen=int(fields[8]),
        seq=fields[9],
        qual=fields[10],
        tags=tags,
    )


def iter_sam(handle: TextIO) -> Iterator[SamRecord]:
    """Iterate alignment records from an open SAM text stream, skipping headers."""
    for line in handle:
        if not line or line.startswith("@"):
            continue
        yield parse_sam_line(line)


def read_sam(path: str) -> Iterator[SamRecord]:
    """Iterate alignment records from a SAM file path."""
    with open(path) as fh:
        yield from iter_sam(fh)
