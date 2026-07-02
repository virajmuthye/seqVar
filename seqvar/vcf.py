"""VCF reading, writing, filtering and summary statistics.

Pure Python implementation of a practical subset of VCF 4.2. It parses header
metadata and data records, exposes typed access to the INFO and FORMAT/sample
fields, and can round-trip a file. This is enough for filtering, summarising
and comparing call sets, which is what most downstream analysis needs.

References
----------
VCF v4.2 specification: https://samtools.github.io/hts-specs/VCFv4.2.pdf
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Iterator, Optional


@dataclass
class VcfRecord:
    """A single VCF data line."""

    chrom: str
    pos: int  # 1-based
    id: str
    ref: str
    alt: list[str]
    qual: Optional[float]
    filter: list[str]
    info: dict[str, object]
    format_keys: list[str] = field(default_factory=list)
    samples: list[dict[str, str]] = field(default_factory=list)

    @property
    def is_snp(self) -> bool:
        return len(self.ref) == 1 and all(len(a) == 1 and a != "." for a in self.alt)

    @property
    def is_indel(self) -> bool:
        return any(len(a) != len(self.ref) for a in self.alt if a not in (".", "*"))

    @property
    def is_biallelic(self) -> bool:
        return len(self.alt) == 1

    def key(self) -> tuple[str, int, str, str]:
        """Identity used for comparison: (chrom, pos, ref, first-alt)."""
        return (self.chrom, self.pos, self.ref, self.alt[0] if self.alt else ".")

    def to_line(self) -> str:
        info_str = _format_info(self.info)
        qual_str = "." if self.qual is None else _fmt_number(self.qual)
        filt = ";".join(self.filter) if self.filter else "."
        cols = [
            self.chrom, str(self.pos), self.id or ".", self.ref,
            ",".join(self.alt) if self.alt else ".",
            qual_str, filt, info_str,
        ]
        if self.format_keys:
            cols.append(":".join(self.format_keys))
            for sample in self.samples:
                cols.append(":".join(sample.get(k, ".") for k in self.format_keys))
        return "\t".join(cols)


def _fmt_number(x: float) -> str:
    """Format a number without a trailing '.0' for integers."""
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x)


def _parse_info(info_str: str) -> dict[str, object]:
    info: dict[str, object] = {}
    if info_str in (".", ""):
        return info
    for entry in info_str.split(";"):
        if "=" in entry:
            k, v = entry.split("=", 1)
            info[k] = v
        else:
            info[entry] = True  # flag
    return info


def _format_info(info: dict[str, object]) -> str:
    if not info:
        return "."
    parts = []
    for k, v in info.items():
        if v is True:
            parts.append(k)
        else:
            parts.append(f"{k}={v}")
    return ";".join(parts)


def _parse_qual(q: str) -> Optional[float]:
    if q == ".":
        return None
    val = float(q)
    return val


@dataclass
class VcfFile:
    """An in-memory VCF: header lines, sample names and records."""

    meta: list[str] = field(default_factory=list)  # '##...' lines, no newline
    samples: list[str] = field(default_factory=list)
    records: list[VcfRecord] = field(default_factory=list)

    def __iter__(self) -> Iterator[VcfRecord]:
        return iter(self.records)

    def __len__(self) -> int:
        return len(self.records)

    def write(self, path: str) -> None:
        with open(path, "w") as fh:
            fh.write(self.header_text())
            for rec in self.records:
                fh.write(rec.to_line() + "\n")

    def header_text(self) -> str:
        meta = self.meta[:]
        if not any(m.startswith("##fileformat") for m in meta):
            meta.insert(0, "##fileformat=VCFv4.2")
        cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if self.samples:
            cols += ["FORMAT", *self.samples]
        return "\n".join(meta) + "\n" + "\t".join(cols) + "\n"


def read_vcf(path: str) -> VcfFile:
    """Parse a VCF file (plain text) into a :class:`VcfFile`."""
    vcf = VcfFile()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                vcf.meta.append(line)
            elif line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) > 9:
                    vcf.samples = cols[9:]
            else:
                vcf.records.append(_parse_record(line, vcf.samples))
    return vcf


def _parse_record(line: str, sample_names: list[str]) -> VcfRecord:
    cols = line.split("\t")
    if len(cols) < 8:
        raise ValueError(f"VCF data line has {len(cols)} columns, expected >= 8")
    chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
    format_keys: list[str] = []
    samples: list[dict[str, str]] = []
    if len(cols) > 8:
        format_keys = cols[8].split(":")
        for sample_col in cols[9:]:
            values = sample_col.split(":")
            samples.append(dict(zip(format_keys, values)))
    return VcfRecord(
        chrom=chrom,
        pos=int(pos),
        id=vid,
        ref=ref,
        alt=[] if alt == "." else alt.split(","),
        qual=_parse_qual(qual),
        filter=[] if filt in (".", "") else filt.split(";"),
        info=_parse_info(info),
        format_keys=format_keys,
        samples=samples,
    )


# --- operations -----------------------------------------------------------

def filter_vcf(
    vcf: VcfFile,
    min_qual: float | None = None,
    pass_only: bool = False,
    snps_only: bool = False,
    indels_only: bool = False,
    min_depth: int | None = None,
) -> VcfFile:
    """Return a new VcfFile keeping records that satisfy all given criteria."""
    out = VcfFile(meta=vcf.meta[:], samples=vcf.samples[:])
    for rec in vcf.records:
        if min_qual is not None and (rec.qual is None or rec.qual < min_qual):
            continue
        if pass_only and rec.filter not in ([], ["PASS"]):
            continue
        if snps_only and not rec.is_snp:
            continue
        if indels_only and not rec.is_indel:
            continue
        if min_depth is not None:
            dp = rec.info.get("DP")
            if dp is None or int(dp) < min_depth:
                continue
        out.records.append(rec)
    return out


@dataclass
class VcfStats:
    total: int = 0
    snps: int = 0
    indels: int = 0
    transitions: int = 0
    transversions: int = 0
    multiallelic: int = 0

    @property
    def ts_tv(self) -> float:
        return self.transitions / self.transversions if self.transversions else 0.0

    def summary(self) -> str:
        return (
            f"records        : {self.total}\n"
            f"SNPs           : {self.snps}\n"
            f"indels         : {self.indels}\n"
            f"multiallelic   : {self.multiallelic}\n"
            f"transitions    : {self.transitions}\n"
            f"transversions  : {self.transversions}\n"
            f"Ts/Tv          : {self.ts_tv:.2f}"
        )


_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def compute_stats(vcf: VcfFile) -> VcfStats:
    """Summarise a VCF: SNP/indel counts and the Ts/Tv ratio.

    Ts/Tv (transitions over transversions) is a standard sanity check; whole-
    genome human call sets typically sit around 2.0-2.1, and a value far from
    that hints at false positives.
    """
    st = VcfStats()
    for rec in vcf.records:
        st.total += 1
        if not rec.is_biallelic:
            st.multiallelic += 1
        if rec.is_snp:
            st.snps += 1
            for a in rec.alt:
                if (rec.ref, a) in _TRANSITIONS:
                    st.transitions += 1
                else:
                    st.transversions += 1
        elif rec.is_indel:
            st.indels += 1
    return st


def _build_parser():
    p = argparse.ArgumentParser(
        prog="seqvar-vcf",
        description="Filter and summarise VCF files.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("stats", help="Print SNP/indel/Ts-Tv statistics.")
    s.add_argument("input")

    f = sub.add_parser("filter", help="Filter records and write a new VCF.")
    f.add_argument("input")
    f.add_argument("-o", "--output", required=True)
    f.add_argument("--min-qual", type=float, default=None)
    f.add_argument("--pass-only", action="store_true")
    f.add_argument("--snps-only", action="store_true")
    f.add_argument("--indels-only", action="store_true")
    f.add_argument("--min-depth", type=int, default=None)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    vcf = read_vcf(args.input)
    if args.cmd == "stats":
        print(compute_stats(vcf).summary())
    elif args.cmd == "filter":
        out = filter_vcf(
            vcf,
            min_qual=args.min_qual,
            pass_only=args.pass_only,
            snps_only=args.snps_only,
            indels_only=args.indels_only,
            min_depth=args.min_depth,
        )
        out.write(args.output)
        print(f"Kept {len(out)}/{len(vcf)} records -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
