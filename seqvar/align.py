"""Read alignment (hybrid: orchestrates industry-standard command-line tools).

Aligning short reads to a reference is a solved problem, and re-implementing a
Burrows-Wheeler aligner in Python would be slower and less correct than the
tools everyone already trusts. So this module is deliberately a thin, well-
behaved orchestration layer over ``bwa`` and ``samtools``: it builds the
reference index, runs the aligner, and produces a coordinate-sorted, indexed
BAM ready for variant calling.

It checks that the external tools exist before running, streams their stderr
back on failure, and returns structured results, which is the part that
actually saves you time in a pipeline.

If bwa/samtools are not installed, the functions raise ``ToolNotFoundError``
with an actionable message. The rest of the toolkit (QC, variant calling, VCF
handling) works on SAM/BAM produced by any aligner, so this module is optional.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass


class ToolNotFoundError(RuntimeError):
    """Raised when a required external executable is not on PATH."""


class AlignmentError(RuntimeError):
    """Raised when an external alignment/indexing command fails."""


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise ToolNotFoundError(
            f"'{tool}' was not found on your PATH. Install it, e.g. "
            f"`conda install -c bioconda {tool}` or `apt-get install {tool}`."
        )
    return path


def _run(cmd: list[str], stdout=None) -> subprocess.CompletedProcess:
    """Run a command, raising AlignmentError with captured stderr on failure."""
    result = subprocess.run(
        cmd,
        stdout=stdout,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise AlignmentError(
            f"Command failed ({' '.join(cmd)}):\n{result.stderr.strip()}"
        )
    return result


def tools_available() -> dict[str, bool]:
    """Report which external aligner tools are available on PATH."""
    return {tool: shutil.which(tool) is not None for tool in ("bwa", "samtools")}


def index_reference(reference: str) -> str:
    """Build a bwa index for ``reference`` (no-op if it already exists)."""
    _require("bwa")
    if not os.path.exists(reference):
        raise FileNotFoundError(reference)
    if os.path.exists(reference + ".bwt"):
        return reference
    _run(["bwa", "index", reference])
    return reference


@dataclass
class AlignmentResult:
    bam: str
    reference: str
    n_input_files: int
    sorted: bool = True
    indexed: bool = True


def align_reads(
    reference: str,
    reads: list[str],
    output_bam: str,
    read_group: str | None = None,
    threads: int = 1,
    sort: bool = True,
    index: bool = True,
) -> AlignmentResult:
    """Align FASTQ reads to ``reference`` and produce a sorted, indexed BAM.

    Parameters
    ----------
    reference : path to the reference FASTA (will be bwa-indexed if needed).
    reads : one FASTQ (single-end) or two FASTQs (paired-end).
    output_bam : path of the BAM to write.
    read_group : optional @RG string, e.g. '@RG\\tID:s1\\tSM:sample1\\tPL:ILLUMINA'.
    threads : threads to pass to bwa mem and samtools sort.

    Returns an :class:`AlignmentResult`. Requires bwa and samtools on PATH.
    """
    _require("bwa")
    _require("samtools")
    if len(reads) not in (1, 2):
        raise ValueError("Provide one FASTQ (single-end) or two (paired-end).")
    for r in reads:
        if not os.path.exists(r):
            raise FileNotFoundError(r)

    index_reference(reference)

    mem_cmd = ["bwa", "mem", "-t", str(threads)]
    if read_group:
        mem_cmd += ["-R", read_group]
    mem_cmd += [reference, *reads]

    # bwa mem -> samtools sort (or view) via a pipe, so we never materialise SAM.
    bwa_proc = subprocess.Popen(
        mem_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False
    )
    if sort:
        down_cmd = ["samtools", "sort", "-@", str(threads), "-o", output_bam, "-"]
    else:
        down_cmd = ["samtools", "view", "-b", "-o", output_bam, "-"]

    down_proc = subprocess.run(
        down_cmd, stdin=bwa_proc.stdout, stderr=subprocess.PIPE
    )
    _, bwa_err = bwa_proc.communicate()
    if bwa_proc.returncode != 0:
        raise AlignmentError(f"bwa mem failed:\n{bwa_err.decode(errors='replace')}")
    if down_proc.returncode != 0:
        raise AlignmentError(
            f"samtools failed:\n{down_proc.stderr.decode(errors='replace')}"
        )

    if index and sort:
        _run(["samtools", "index", output_bam])

    return AlignmentResult(
        bam=output_bam,
        reference=reference,
        n_input_files=len(reads),
        sorted=sort,
        indexed=index and sort,
    )


def _build_parser():
    import argparse

    p = argparse.ArgumentParser(
        prog="seqvar-align",
        description="Align FASTQ reads to a reference (bwa mem + samtools sort).",
    )
    p.add_argument("reference", help="Reference FASTA.")
    p.add_argument("reads", nargs="+", help="One or two FASTQ files.")
    p.add_argument("-o", "--output", required=True, help="Output sorted BAM.")
    p.add_argument("-R", "--read-group", help="@RG read-group string.")
    p.add_argument("-t", "--threads", type=int, default=1)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    result = align_reads(
        args.reference,
        args.reads,
        args.output,
        read_group=args.read_group,
        threads=args.threads,
    )
    print(f"Wrote {result.bam} (sorted={result.sorted}, indexed={result.indexed})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
