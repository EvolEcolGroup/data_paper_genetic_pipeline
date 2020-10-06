# Genetic pipeline to process genetic data from FASTQ files to a pairwise pi matrix

There are four main steps to go from FASTQ files to a pairwise pi matrix:

1. FASTQ -> BAM.
2. BAM -> VCF.
3. Merge VCF files.
4. Calculate pairwise pi and within sample heterozygosity.

## 1. FASTQ to BAM

##### Dependencies
1. Python v3.
2. [SAMtools](http://www.htslib.org/download/) (v1.9).
3. [BWA](http://bio-bwa.sourceforge.net/) (v0.7.12).
4. [Cutadapt](http://cutadapt.readthedocs.io/en/stable/) (v1.9.1).
5. softclip.py (provided in [scripts](XXX) folder).
5. [Picard](https://broadinstitute.github.io/picard/) (v2.9.2).
6. [GenomeAnalysisToolkit](https://software.broadinstitute.org/gatk/) (v3.7).
