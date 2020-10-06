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
5. softclip.py (provided in [scripts](XXX) folder). **Add scripts into folder, copy scripts from cluster**
5. [Picard](https://broadinstitute.github.io/picard/) (v2.9.2).
6. [GenomeAnalysisToolkit](https://software.broadinstitute.org/gatk/) (v3.7).

##### Inputs and top level scripts

1. A sample information file (csv format).
2. A parameter and links file (csv format).
3. The alignment script ["make_aligment_script.py"](XXX) which makes an executable bash script. **Add scripts into folder, copy scripts from cluster**

##### Sample information file

The sample information file should have four columns of information in csv format: 
sample id, input gzipped fastq file, library index and library ID.

*Example sample information file*
```
sample_id,fastq_file,barcode,library
KK1,KK1_1_1.fastq.gz,AGACTCC,KK1_AGACTCC_hiseq17
KK1,KK1_2_1.fastq.gz,GTTACCG,KK1_GTTACCG_hiseq17
NE5,NE5_1_1.fastq.gz,ACGGCAG,NE5_ACGGCAG_hiseq17
NE5,NE5_2_1.fastq.gz,GACCGAT,NE5_GACCGAT_hiseq17
ZVEJ31,ZVEJ31_1_1.fastq.gz,GGTAACT,ZVEJ31_GGTAACT_hiseq17
```

##### Parameters and links file

This file should contain links to program paths and the reference genome. It should also contain parameter information for Cutadapt, BWA and the mapping quality threshold for read filtering.

*Example parameters and links file* **Copy one param file from the cluster**
```
paths:,
Reference:,/home/eppie/hg19_with_rCRS_ordered/hg19_with_rCRS_ordered.fa
Picard:,/home/eppie/programs/picard/build/libs/picard.jar
GATK:,/home/lara/Software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
softclip.py:,/home/eppie/hiseq_2017_alignment/programs/softclip.py
,
parameters:,
Cutadapt:,-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 34
BWA aln:,-l 1000 -t 12 -n0.01 -o2
mapping quality threshold:,20
```

##### Alignment script
This python script  ["make_aligment_script.py"](XXX) (**Add script in the folder**) creates an alignment script which should be run with bash. It works on gzipped fastq files which have already been split by library index. Adapters are trimmed from the end of reads with Cutadapt using the parameters set in the parameter file. The reads are aligned using BWA (parameters set in parameter file) and the files are sorted and indexed using SAMtools. All steps apart from the BWA aln step are run in parallel so you will need as many cores as input files. If this is an issue you should split your sample information file into smaller files and run commands in series. Only run one alignment bash script per folder (or else output statistics files will be overwritten). The number of threads to use for the BWA aln step should be set in the parameter file. 

Reads from the same sample are merged using Picard. Reads are filtering by a mapping quality threshold set in the parameter file. Duplicate reads are removed using SAMtools. Indels are realigned using The Genome Analysis Toolkit and 2bp are softclip from the start and ends of reads. 

Basic alignment statistics are provided in "alignment_stats" output files. 

The robusticity of this pipeline has not been tested extensively. It is recommended to use the template input files provided and to check the alignment bash script before running it.

**Running the python and bash scripts**

```
python make_aligment_script.py <sample_information_file.csv> <parameters_and_links_file.csv> <output_file.sh>
bash output_file.sh
```

### 2. BAM to VCF

Convert from BAM files to filtered VCF files.

##### Dependencies
1. Python (v2.7.12)
2. [SAMtools](http://www.htslib.org/download/) (v1.9).
3. [GenomeAnalysisToolkit](https://software.broadinstitute.org/gatk/) (v3.7).
4. [VCFtools](http://vcftools.sourceforge.net/perl_module.html) (v0.1.5)
5. filter_GC.py (provided in [scripts](XXX) folder) **Add script to the folder**
6. filter_vcf_minDepth_maxDepth.py (provided in [scripts](XXX) folder) **Add script to the folder**

##### Inputs and top level scripts
1. A parameter and links file (csv format).
2. The [filter_ancient_bams.sh](https://gitlab.com/manica-group/genetics_pipeline/tree/master/scripts) script.
