"""SAM/BAM statistics, filtering and coverage.

Works natively on SAM text (via :mod:`seqvar.sam`). For binary BAM it prefers
pysam if installed, and otherwise transparently streams the file through
``samtools view``. This keeps the core logic pure Python and dependency-free
while still handling real BAM files when the ecosystem tools are present.

Provided operations
--------------------
- ``flagstat``  : samtools-flagstat-style counts (total, mapped, duplicates...).
- ``filter_records`` : filter by MAPQ, flags, and mapped status.
- ``coverage``  : per-position depth over a reference region.
"""

from __future__ import annotations

import shutil
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from typing import Iterator

from . import sam
from .sam import SamRecord


@contextmanager
def open_alignments(path: str) -> Iterator[Iterator[SamRecord]]:
    """Yield an iterator of :class:`SamRecord` from a SAM or BAM path.

    - ``.sam`` (or anything non-.bam): parsed directly as text.
    - ``.bam``: uses pysam if available, else ``samtools view``.
    """
    if not path.endswith(".bam"):
        with open(path) as fh:
            yield sam.iter_sam(fh)
        return

    # BAM: try pysam first for a pure in-process path.
    try:
        import pysam  # type: ignore

        af = pysam.AlignmentFile(path, "rb")
        try:
            yield (_from_pysam(r, af) for r in af)
        finally:
            af.close()
        return
    except ImportError:
        pass

    # Fall back to samtools view streaming SAM text.
    if shutil.which("samtools") is None:
        raise RuntimeError(
            f"Reading BAM {path!r} requires either the 'pysam' package or "
            "'samtools' on PATH. Neither was found. Convert to SAM first, or "
            "install one of them."
        )
    proc = subprocess.Popen(
        ["samtools", "view", path], stdout=subprocess.PIPE, text=True
    )
    try:
        yield sam.iter_sam(proc.stdout)
    finally:
        proc.stdout.close()
        proc.wait()


def _from_pysam(read, af) -> SamRecord:
    """Adapt a pysam AlignedSegment to our SamRecord dataclass."""
    return SamRecord(
        qname=read.query_name or "*",
        flag=read.flag,
        rname=(af.get_reference_name(read.reference_id)
               if read.reference_id >= 0 else "*"),
        pos=read.reference_start + 1 if read.reference_start is not None else 0,
        mapq=read.mapping_quality,
        cigar=read.cigarstring or "*",
        rnext="*",
        pnext=0,
        tlen=read.template_length or 0,
        seq=read.query_sequence or "*",
        qual="*",
        tags={},
    )


@dataclass
class FlagStat:
    total: int = 0
    mapped: int = 0
    unmapped: int = 0
    secondary: int = 0
    supplementary: int = 0
    duplicates: int = 0
    paired: int = 0
    properly_paired: int = 0
    read1: int = 0
    read2: int = 0

    @property
    def mapped_fraction(self) -> float:
        return self.mapped / self.total if self.total else 0.0

    def to_dict(self) -> dict:
        d = asdict(self)
        d["mapped_fraction"] = self.mapped_fraction
        return d

    def summary(self) -> str:
        pct = self.mapped_fraction * 100
        return (
            f"{self.total} total\n"
            f"{self.mapped} mapped ({pct:.2f}%)\n"
            f"{self.unmapped} unmapped\n"
            f"{self.secondary} secondary\n"
            f"{self.supplementary} supplementary\n"
            f"{self.duplicates} duplicates\n"
            f"{self.paired} paired in sequencing\n"
            f"{self.properly_paired} properly paired\n"
            f"{self.read1} read1\n"
            f"{self.read2} read2"
        )


def flagstat(path: str) -> FlagStat:
    """Compute samtools-flagstat-style counts for a SAM/BAM file."""
    fs = FlagStat()
    with open_alignments(path) as records:
        for r in records:
            fs.total += 1
            if r.is_unmapped:
                fs.unmapped += 1
            else:
                fs.mapped += 1
            if r.is_secondary:
                fs.secondary += 1
            if r.is_supplementary:
                fs.supplementary += 1
            if r.is_duplicate:
                fs.duplicates += 1
            if r.flag & sam.FLAG_PAIRED:
                fs.paired += 1
                if r.is_proper_pair:
                    fs.properly_paired += 1
                if r.flag & sam.FLAG_READ1:
                    fs.read1 += 1
                if r.flag & sam.FLAG_READ2:
                    fs.read2 += 1
    return fs


def filter_records(
    in_path: str,
    out_sam: str,
    min_mapq: int = 0,
    exclude_flags: int = 0,
    require_flags: int = 0,
    mapped_only: bool = False,
) -> int:
    """Filter alignments and write a headerless SAM. Returns count kept.

    exclude_flags / require_flags mirror samtools view -F / -f.
    """
    kept = 0
    with open_alignments(in_path) as records, open(out_sam, "w") as out:
        for r in records:
            if mapped_only and r.is_unmapped:
                continue
            if r.mapq < min_mapq:
                continue
            if exclude_flags and (r.flag & exclude_flags):
                continue
            if require_flags and (r.flag & require_flags) != require_flags:
                continue
            out.write(_to_sam_line(r) + "\n")
            kept += 1
    return kept


def _to_sam_line(r: SamRecord) -> str:
    core = [
        r.qname, str(r.flag), r.rname, str(r.pos), str(r.mapq), r.cigar,
        r.rnext, str(r.pnext), str(r.tlen), r.seq, r.qual,
    ]
    tags = [f"{k}:{'i' if v.lstrip('-').isdigit() else 'Z'}:{v}"
            for k, v in r.tags.items()]
    return "\t".join(core + tags)


def coverage(
    path: str,
    reference_name: str,
    start: int = 1,
    end: int | None = None,
    min_mapq: int = 0,
) -> list[int]:
    """Per-position read depth over ``reference_name`` in [start, end].

    Coordinates are 1-based inclusive. If ``end`` is None it is inferred from
    the rightmost aligned base observed. Positions with no coverage are 0.
    Secondary/supplementary/unmapped reads are ignored.
    """
    depths: dict[int, int] = {}
    max_ref = start
    with open_alignments(path) as records:
        for r in records:
            if r.is_unmapped or r.rname != reference_name:
                continue
            if r.is_secondary or r.is_supplementary:
                continue
            if r.mapq < min_mapq:
                continue
            for _qi, rp in r.aligned_pairs():
                if rp is None:
                    continue
                if rp < start or (end is not None and rp > end):
                    continue
                depths[rp] = depths.get(rp, 0) + 1
                max_ref = max(max_ref, rp)
    last = end if end is not None else max_ref
    return [depths.get(pos, 0) for pos in range(start, last + 1)]


def _build_parser():
    import argparse

    p = argparse.ArgumentParser(
        prog="seqvar-bam",
        description="SAM/BAM statistics, filtering and coverage.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    fs = sub.add_parser("flagstat", help="Print alignment flag statistics.")
    fs.add_argument("input", help="SAM or BAM file.")

    fl = sub.add_parser("filter", help="Filter alignments to a SAM file.")
    fl.add_argument("input")
    fl.add_argument("-o", "--output", required=True)
    fl.add_argument("-q", "--min-mapq", type=int, default=0)
    fl.add_argument("-F", "--exclude-flags", type=int, default=0)
    fl.add_argument("-f", "--require-flags", type=int, default=0)
    fl.add_argument("--mapped-only", action="store_true")

    cov = sub.add_parser("coverage", help="Per-position depth over a region.")
    cov.add_argument("input")
    cov.add_argument("-r", "--reference", required=True)
    cov.add_argument("--start", type=int, default=1)
    cov.add_argument("--end", type=int, default=None)
    cov.add_argument("-q", "--min-mapq", type=int, default=0)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.cmd == "flagstat":
        print(flagstat(args.input).summary())
    elif args.cmd == "filter":
        n = filter_records(
            args.input, args.output,
            min_mapq=args.min_mapq,
            exclude_flags=args.exclude_flags,
            require_flags=args.require_flags,
            mapped_only=args.mapped_only,
        )
        print(f"Kept {n} records -> {args.output}")
    elif args.cmd == "coverage":
        depths = coverage(
            args.input, args.reference,
            start=args.start, end=args.end, min_mapq=args.min_mapq,
        )
        for i, d in enumerate(depths, start=args.start):
            print(f"{args.reference}\t{i}\t{d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
