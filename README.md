# seqvar

A compact, readable Python toolkit for the core steps of a small-variant workflow: going from **sequencer reads (FASTQ)** through **alignment** and **SAM/BAM manipulation** to **variant calling** and **VCF handling and comparison**.

The design goal is clarity, not another wrapper around a monolithic pipeline. The variant caller and the VCF logic are implemented in **pure Python** so you can read exactly how a call is made and how two call sets are compared. Alignment, where re-implementing a Burrows-Wheeler aligner would be pointless, is a thin **hybrid** layer that orchestrates the tools everyone already trusts (`bwa` + `samtools`).

The core toolkit has **no third-party dependencies** (standard library only). It runs, and its full test suite passes, without installing anything.

---

## Why this exists

Most variant-calling code is either a black box or a shell script gluing binaries together. This repo shows the actual mechanics: how a pileup is built from CIGAR strings, how an alternate allele clears depth/fraction/quality thresholds, how a Phred quality is derived, and how you score a caller against a truth set with precision and recall. It is meant to be legible and correct on small data, and honest about where it is a teaching model rather than a production caller.

## The pipeline

```
 FASTQ            reference          SAM / BAM           VCF                VCF
  reads   ─QC─▶    +  ─align─▶      alignments  ─call─▶   calls   ─compare─▶  metrics
                (bwa+samtools)     (stats/filter/       (pileup       (precision,
   fastqc          align            coverage)            SNV+indel)    recall, F1)
                                     bam                  call          compare
```

Each stage is a module with a library API and a CLI subcommand.

| Stage | Module | Implementation | What it does |
|-------|--------|----------------|--------------|
| QC | `seqvar.fastq_qc` | pure Python | Read counts, length distribution, per-position mean Phred quality, GC and N content |
| Align | `seqvar.align` | hybrid (`bwa`+`samtools`) | Index reference, align FASTQ, emit sorted+indexed BAM |
| SAM parsing | `seqvar.sam` | pure Python | SAM records, FLAG helpers, CIGAR walking, aligned-pair iteration |
| BAM ops | `seqvar.bam_utils` | pure Python (pysam optional) | `flagstat`, flag/MAPQ filtering, per-position coverage |
| Call | `seqvar.call_variants` | pure Python | Pileup-based SNV and short-indel caller → VCF |
| VCF | `seqvar.vcf` | pure Python | Parse/write VCF, filter, SNP/indel/Ts-Tv stats |
| Compare | `seqvar.vcf_compare` | pure Python | Concordance of two VCFs: TP/FP/FN, precision/recall/F1 |

## Install

```bash
git clone https://github.com/virajmuthye/seq2variant.git
cd seq2variant
pip install -e .            # installs the `seqvar` command; core needs no deps
```

Optional extras:

```bash
pip install -e ".[bam]"     # pysam, for native BAM reading (else samtools is used)
pip install -e ".[test]"    # pytest
```

`bwa` and `samtools` are only needed for the `align` stage. Everything else operates on SAM/BAM produced by any aligner.

## Quick start

The repo ships a tiny synthetic dataset generator with two known variants (a SNP and a 2 bp deletion) spiked in, so you can run the whole thing end to end with a checkable answer and no external tools:

```bash
bash examples/run_pipeline.sh
```

Or step by step:

```bash
# 0. generate the demo data (ref.fa, reads.fastq, alignment.sam, truth.vcf)
python -m examples.make_test_data data

# 1. QC the reads
seqvar fastqc data/reads.fastq

# 2. alignment statistics
seqvar bam flagstat data/alignment.sam

# 3. call variants
seqvar call data/alignment.sam -r data/ref.fa -o data/calls.vcf

# 4. summarise the calls
seqvar vcf stats data/calls.vcf

# 5. benchmark against truth
seqvar compare data/truth.vcf data/calls.vcf
```

Step 5 prints a perfect recovery on the demo data:

```
true positives  (shared)     : 2
false positives (test-only)  : 0
false negatives (truth-only) : 0
precision                    : 1.0000
recall (sensitivity)         : 1.0000
F1 score                     : 1.0000
```

## Aligning real reads

With `bwa` and `samtools` installed, produce a sorted, indexed BAM from FASTQs:

```bash
seqvar align data/ref.fa sample_R1.fastq.gz sample_R2.fastq.gz \
    -o sample.bam -R '@RG\tID:s1\tSM:sample1\tPL:ILLUMINA' -t 4
seqvar call sample.bam -r data/ref.fa -o sample.vcf
```

## Library usage

```python
from seqvar.fastq_qc import qc_fastq
from seqvar.call_variants import call_variants, CallParams
from seqvar.vcf import read_vcf, compute_stats
from seqvar.vcf_compare import compare

print(qc_fastq("reads.fastq").summary())

vcf = call_variants(
    "alignment.sam", "ref.fa",
    params=CallParams(min_depth=8, min_alt_fraction=0.25),
)
vcf.write("calls.vcf")

result = compare(read_vcf("truth.vcf"), read_vcf("calls.vcf"))
print(result.precision, result.recall, result.f1)
```

## How the variant caller works

The caller is intentionally transparent:

1. **Pileup.** Every alignment is walked by its CIGAR. Matched bases are tallied per reference position (filtered by base quality and mapping quality); insertions are anchored to the preceding reference base; deletions are recorded by length.
2. **Candidate selection.** At each position the most-supported non-reference observation (substitution, insertion or deletion) becomes the candidate alternate allele.
3. **Thresholds.** A candidate is emitted only if it clears minimum depth, minimum supporting-allele count, and minimum allele fraction.
4. **Quality.** `QUAL` is the Phred-scaled binomial tail probability of seeing at least this many alternate observations under a background sequencing-error rate. Genotype is a simple heuristic: allele fraction ≥ 0.8 → homozygous (`1/1`), otherwise heterozygous (`0/1`).

Indels are written in standard VCF form (deletions anchored on the preceding base), so calls line up with real truth sets and with the normalization used by `vcf_compare`.

**Scope, honestly:** this is a legible pileup caller, not GATK/DeepVariant. It has no local realignment, no full genotype-likelihood/HWE model, no base-quality recalibration, and no population priors. It is meant for learning, small datasets, and as a clear reference implementation.

## VCF comparison

`vcf_compare` matches records by a **normalized** `(CHROM, POS, REF, ALT)` key. Alleles are left-normalized (shared prefix/suffix trimmed) so equivalent indel spellings collapse to the same key before matching. It reports true positives (shared), false positives (test-only), false negatives (truth-only), and the derived precision, recall and F1. False-positive and false-negative records can be written back out as VCFs (`--fp-out`, `--fn-out`) for inspection.

## Running the tests

The suite is written with the standard-library `unittest` framework, so it runs with or without `pytest`:

```bash
python -m pytest                                   # if pytest is installed
python -m unittest discover -s tests -p "test_*.py"  # stdlib only
```

30 tests cover CIGAR/FLAG parsing, FASTQ QC metrics, flagstat/coverage/filtering, the SNV and indel calls (against the spiked truth), VCF round-tripping and Ts/Tv statistics, and the precision/recall comparison logic. Test fixtures are regenerated deterministically from `examples/make_test_data.py`, so nothing binary is committed.

## Repository layout

```
seq2variant/
├── seqvar/
│   ├── fastq_qc.py       # FASTQ QC metrics (pure Python)
│   ├── fasta.py          # minimal FASTA reader/writer
│   ├── sam.py            # SAM record + CIGAR/FLAG logic (pure Python)
│   ├── align.py          # bwa + samtools orchestration (hybrid)
│   ├── bam_utils.py      # flagstat / filter / coverage (pysam optional)
│   ├── call_variants.py  # pileup-based SNV + indel caller (pure Python)
│   ├── vcf.py            # VCF parse/write/filter/stats (pure Python)
│   ├── vcf_compare.py    # truth-vs-test concordance (pure Python)
│   └── cli.py            # unified `seqvar` command
├── tests/                # unittest suite (runs under pytest too)
├── examples/
│   ├── make_test_data.py # deterministic synthetic dataset generator
│   └── run_pipeline.sh   # end-to-end demo
├── data/                 # generated demo inputs (ref, reads, SAM, truth VCF)
├── pyproject.toml
├── LICENSE               # MIT
└── README.md
```

## File formats

Implemented against the public HTS specifications: [SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) and [VCF v4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Coordinates follow the conventions of each format (SAM/VCF are 1-based).

## License

MIT — see [LICENSE](LICENSE).
