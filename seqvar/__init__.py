"""seqvar: a small pure-Python toolkit for going from sequencer reads to variant calls.

Modules
-------
fastq_qc      : FASTQ quality-control metrics (pure Python).
align         : reference indexing / read alignment (hybrid: wraps bwa + samtools).
sam           : lightweight SAM record parser and helpers (pure Python).
bam_utils     : SAM/BAM statistics, filtering and coverage (pysam optional).
call_variants : pileup-based SNV/indel caller (pure Python).
vcf           : VCF reader/writer and records (pure Python).
vcf_compare   : concordance between two VCF call sets (pure Python).
"""

__version__ = "0.1.0"

__all__ = [
    "fastq_qc",
    "align",
    "sam",
    "bam_utils",
    "call_variants",
    "vcf",
    "vcf_compare",
]
