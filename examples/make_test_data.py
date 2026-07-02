"""Generate a tiny, self-consistent synthetic dataset for seqvar.

Produces, under the target directory (default: ../data):

    ref.fa        reference sequence (one contig, ~240 bp)
    reads.fastq   simulated reads (FASTQ)
    alignment.sam aligned reads with correct POS/CIGAR (no external aligner)
    truth.vcf     the ground-truth variants that were spiked in

The variants spiked in are a SNP and a 2 bp deletion. Reads are generated
directly against the reference with correct CIGAR strings, so the whole
pipeline (call -> compare against truth) runs without bwa/samtools installed
and yields a known, checkable answer.

Everything is deterministic (fixed RNG seed) so tests are reproducible.
"""

from __future__ import annotations

import os
import random

# Make the package importable when run as a plain script.
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from seqvar.fasta import write_fasta  # noqa: E402

CONTIG = "chr_test"
REF_LEN = 240
READ_LEN = 60
STEP = 6  # start-position spacing between tiled reads -> ~10x coverage
SEED = 42

# 1-based variant definitions spiked into the sampled reads.
SNP_POS = 80          # reference position of the SNP
DEL_START = 150       # first deleted reference base (1-based)
DEL_LEN = 2           # number of bases deleted


def make_reference(rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(REF_LEN))


def _snp_alt(ref_base: str, rng: random.Random) -> str:
    return rng.choice([b for b in "ACGT" if b != ref_base])


def generate(out_dir: str) -> dict:
    rng = random.Random(SEED)
    ref = make_reference(rng)

    ref_snp_base = ref[SNP_POS - 1]
    alt_base = _snp_alt(ref_snp_base, rng)

    # --- write reference ---
    os.makedirs(out_dir, exist_ok=True)
    write_fasta(os.path.join(out_dir, "ref.fa"), {CONTIG: ref})

    # --- write truth VCF (standard, preceding-base-anchored deletion) ---
    del_anchor_pos = DEL_START - 1
    del_ref = ref[del_anchor_pos - 1: DEL_START - 1 + DEL_LEN]  # anchor + deleted
    del_alt = ref[del_anchor_pos - 1]
    truth_lines = [
        "##fileformat=VCFv4.2",
        "##source=make_test_data",
        f"##contig=<ID={CONTIG},length={REF_LEN}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        f"{CONTIG}\t{SNP_POS}\t.\t{ref_snp_base}\t{alt_base}\t.\tPASS\tTYPE=SNP",
        f"{CONTIG}\t{del_anchor_pos}\t.\t{del_ref}\t{del_alt}\t.\tPASS\tTYPE=DEL",
    ]
    with open(os.path.join(out_dir, "truth.vcf"), "w") as fh:
        fh.write("\n".join(truth_lines) + "\n")

    # --- simulate reads + SAM ---
    fastq_records: list[tuple[str, str, str]] = []
    sam_records: list[str] = []
    read_id = 0

    for start in range(1, REF_LEN - READ_LEN + 2, STEP):  # 1-based read start
        end = start + READ_LEN - 1  # inclusive ref end for a pure-match read
        spans_del = start <= DEL_START and end >= DEL_START + DEL_LEN - 1

        if spans_del:
            # Read carries the deletion: left match, D, right match.
            left_len = DEL_START - start           # bases before deletion
            right_len = READ_LEN - left_len          # remaining read bases
            left = ref[start - 1: DEL_START - 1]
            right = ref[DEL_START - 1 + DEL_LEN:
                        DEL_START - 1 + DEL_LEN + right_len]
            seq = left + right
            cigar = f"{left_len}M{DEL_LEN}D{right_len}M"
        else:
            seq = list(ref[start - 1: end])
            if start <= SNP_POS <= end:
                seq[SNP_POS - start] = alt_base  # inject the SNP allele
            seq = "".join(seq)
            cigar = f"{len(seq)}M"

        qname = f"read{read_id}"
        qual = "I" * len(seq)  # Phred ~40, well above the caller's threshold
        fastq_records.append((qname, seq, qual))
        # FLAG 0 = mapped, forward strand, single-end.
        sam_records.append(
            "\t".join([qname, "0", CONTIG, str(start), "60", cigar,
                       "*", "0", "0", seq, qual])
        )
        read_id += 1

    with open(os.path.join(out_dir, "reads.fastq"), "w") as fh:
        for name, seq, qual in fastq_records:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")

    with open(os.path.join(out_dir, "alignment.sam"), "w") as fh:
        fh.write(f"@HD\tVN:1.6\tSO:coordinate\n")
        fh.write(f"@SQ\tSN:{CONTIG}\tLN:{REF_LEN}\n")
        for line in sam_records:
            fh.write(line + "\n")

    return {
        "contig": CONTIG,
        "ref_len": REF_LEN,
        "snp": (SNP_POS, ref_snp_base, alt_base),
        "deletion": (del_anchor_pos, del_ref, del_alt),
        "n_reads": read_id,
        "out_dir": out_dir,
    }


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(__file__), "..", "data"
    )
    info = generate(os.path.abspath(target))
    print("Generated synthetic dataset:")
    for k, v in info.items():
        print(f"  {k}: {v}")
