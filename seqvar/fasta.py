"""Minimal FASTA reader.

Loads a reference genome into memory as a dict of {name: sequence}. Fine for
the small references used in tests and teaching; for chromosome-scale
references you would want an indexed (.fai) random-access reader instead.
"""

from __future__ import annotations

from typing import Dict


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file into an ordered dict of {sequence_name: sequence}.

    The sequence name is the first whitespace-delimited token after '>'.
    """
    sequences: Dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def write_fasta(path: str, sequences: Dict[str, str], width: int = 60) -> None:
    """Write {name: sequence} to a FASTA file wrapped at ``width`` columns."""
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
