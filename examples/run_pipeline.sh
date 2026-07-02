#!/usr/bin/env bash
# End-to-end demo of the seqvar toolkit on the bundled synthetic dataset.
#
# It runs the full path from reads to an evaluated variant call set. No
# external tools (bwa/samtools) are required because the synthetic SAM is
# generated with correct coordinates already; the alignment step is shown as a
# commented example of how you would produce that SAM/BAM from real FASTQs.
set -euo pipefail

cd "$(dirname "$0")/.."

echo "==> 0. (re)generate the synthetic dataset"
python3 examples/make_test_data.py data

echo; echo "==> 1. FASTQ quality control"
python3 -m seqvar.cli fastqc data/reads.fastq

# echo; echo "==> 2. Align reads (requires bwa + samtools on PATH)"
# python3 -m seqvar.cli align data/ref.fa data/reads.fastq -o data/alignment.bam

echo; echo "==> 2. Alignment statistics"
python3 -m seqvar.cli bam flagstat data/alignment.sam

echo; echo "==> 3. Call variants"
python3 -m seqvar.cli call data/alignment.sam -r data/ref.fa -o data/calls.vcf

echo; echo "==> 4. Summarise the call set"
python3 -m seqvar.cli vcf stats data/calls.vcf

echo; echo "==> 5. Benchmark calls against the truth set"
python3 -m seqvar.cli compare data/truth.vcf data/calls.vcf

echo; echo "Done. Calls written to data/calls.vcf"
