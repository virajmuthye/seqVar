# seqvar — Quick Start

Copy-paste commands. No reading required.

## Install (once)

```bash
cd seq2variant
pip install -e .        # gives you the `seqvar` command, no dependencies
```

## Try it right now

```bash
bash examples/run_pipeline.sh
```

That generates demo data and runs the whole thing: QC → stats → call → compare. Done.

## The 5 commands you'll actually use

```bash
# 1. Check read quality
seqvar fastqc reads.fastq

# 2. Alignment stats (SAM or BAM)
seqvar bam flagstat aligned.bam

# 3. Call variants  ->  VCF
seqvar call aligned.bam -r reference.fa -o calls.vcf

# 4. Summarize a VCF (SNP/indel counts, Ts/Tv)
seqvar vcf stats calls.vcf

# 5. Compare your calls to a truth set (precision/recall/F1)
seqvar compare truth.vcf calls.vcf
```

## Starting from raw FASTQ? (needs bwa + samtools installed)

```bash
seqvar align reference.fa reads_R1.fq reads_R2.fq -o aligned.bam
# then go to step 3 above
```

## Common tweaks

```bash
# stricter variant calling
seqvar call aligned.bam -r reference.fa -o calls.vcf --min-depth 10 --min-alt-fraction 0.25

# keep only PASS, high-quality variants
seqvar vcf filter calls.vcf -o filtered.vcf --pass-only --min-qual 30

# filter a BAM: mapped reads only, MAPQ >= 20
seqvar bam filter aligned.bam -o clean.sam -q 20 --mapped-only

# save the mismatches from a comparison
seqvar compare truth.vcf calls.vcf --fp-out false_positives.vcf --fn-out missed.vcf
```

## In Python instead of the shell

```python
from seqvar.call_variants import call_variants
from seqvar.vcf import read_vcf
from seqvar.vcf_compare import compare

vcf = call_variants("aligned.bam", "reference.fa")
vcf.write("calls.vcf")

r = compare(read_vcf("truth.vcf"), read_vcf("calls.vcf"))
print(r.precision, r.recall, r.f1)
```

## Getting help on any command

```bash
seqvar call --help
seqvar bam --help
```

That's it. For the how-it-works details, see `README.md`.
