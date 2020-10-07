# Genetic pipeline to process genetic data from ancient and modern individuals

The pipeline is divided in two main sections for ancient and modern samples.

## Ancient samples

There are three main steps to go from FASTQ files to a filtered multi-sample vcf file:

1. FASTQ -> BAM.
2. BAM -> VCF.

### 1. FASTQ to BAM

##### Dependencies
1. Python (v3).
2. [SAMtools](http://www.htslib.org/download/) (v1.9).
3. [BWA](http://bio-bwa.sourceforge.net/) (v0.7.12).
4. [Cutadapt](http://cutadapt.readthedocs.io/en/stable/) (v1.9.1).
5. softclip.py (provided in [scripts](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) folder).
5. [Picard](https://broadinstitute.github.io/picard/) (v2.9.2).
6. [GenomeAnalysisToolkit](https://software.broadinstitute.org/gatk/) (v3.7).

##### Inputs and top level scripts

1. A sample information file (csv format).
2. A parameter and links file (csv format).
3. The alignment script ["make_aligment_script.py"](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) which makes an executable bash script.

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

*Example parameters and links file*
```
paths:,
Reference:,/home/pm604/ref_genome/hg19_with_rCRS_ordered/hg19_with_rCRS_ordered.fa
Picard:,/home/pm604/tools/picard/picard.jar
GATK:,/home/pm604/tools/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar
softclip.py:,/home/pm604/scripts/softclip.py
,
parameters:,
Cutadapt:,-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 34 -j 32
BWA aln:,-l 1000 -t 32 -n0.01 -o2
mapping quality threshold:,20
```
The options `-t` and `-j` for BWA and Cutadapt respectively are set to use 32 CPUs. They need to be adjusted depedning on the system.

##### Alignment script
This python script  ["make_aligment_script.py"](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) creates an alignment script which should be run with bash. It works on gzipped fastq files which have already been split by library index. Adapters are trimmed from the end of reads with Cutadapt using the parameters set in the parameter file. The reads are aligned using BWA (parameters set in parameter file) and the files are sorted and indexed using SAMtools. All steps apart from the BWA aln step are run in parallel so you will need as many cores as input files. If this is an issue you should split your sample information file into smaller files and run commands in series. Only run one alignment bash script per folder (or else output statistics files will be overwritten). The number of threads to use for the BWA aln step should be set in the parameter file. 

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
1. Python (v3)
2. [SAMtools](http://www.htslib.org/download/) (v1.9).
3. [GenomeAnalysisToolkit](https://software.broadinstitute.org/gatk/) (v3.7).
4. [VCFtools](http://vcftools.sourceforge.net/perl_module.html) (v0.1.5).
5. filter_GC.py (provided in [scripts](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) folder).
6. filter_vcf_minDepth_maxDepth.py (provided in [scripts](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) folder).

##### Inputs and top level scripts
1. A parameter and links file (csv format).
2. The [filter_ancient_bams.sh](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) script.

##### Parameters and links file

The parameter and links file should contain links to the reference genome and a link to lists of positions you want to call in your sample(s). This list should be split by chromosome. I am using **Gronau regions** (see [README](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/blob/main/Gronau_filters/Gronau_filters.md) ). It should also contain links to GATK, VCFtools and python scripts (provided in  [scripts](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts)). 

The parameter and links file should specify the minimum and maximum genotype depth for filtering. The depth can be written in absolute terms or an X can be used to specify the coverage relative to the average genotype coverage within the Gronau windows e.g. max_depth: 2X is 2 times the average coverage. 

The *clean_up* option specifies whether you want intermediate files to be deleted (Y/N), and *merge* specifies whether you want a single vcf file at the end (Y) or one vcf per chromosome (N).

*Example parameters and links file*

```
data,
snp_list:,/home/pm604/Gronau_filters/Gronau_regions_to_keep_per_chromosome/chr_in_chr_name
reference_genome:,/home/pm604/ref_genome/hg19_with_rCRS_ordered/hg19_with_rCRS_ordered.fa
,
links,
GATK:,/home/pm604/tools/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar
filter_GC.py:,/home/pm604/scripts/filter_GC.py
filter_vcf_minDepth_maxDepth.py:,/home/pm604/scripts/filter_vcf_minDepth_maxDepth.py
vcf-concat:,vcf-concat
,
options,
min_depth:,8
max_depth:,2X
clean_up:,Y
merge:,N
```

##### The filter\_ancient\_bams.sh script

You will need a machine which has at least 23 cores. BAM files are split by chromosome (autosomes only) and genotypes called using Genome Analysis Toolkit Unified Genotyper. The input priors are set to equal for homozygous reference and homozygous alternate genotypes (as done for the Simons Genome Diversity panel). Positions of interest are specified in the parameters and links file. Positions where there is a G with a C called immediately up or downstream are removed as are positions with a C called and a G immediately up or downstream of it.  C and G sites with a missing genotype immediately preceding or proceeding it are also removed. Genotypes are filtered by a minimum and maximum depth threshold as specified in the parameters and links file. The depth can be written in absolute terms or an X can be used to specify the coverage relative to the average genotype coverage within the Gronau windows e.g. max_depth: 2X is 2 times the average coverage. 

**Running the script**

```bash
bash filter_ancient_bams.sh <parameters_and_links.csv>
```


## Modern samples

Mondern samples were processed from bam files. All modern samples were downloaded from the [Cancer Genomics Cloud](https://cgc.sbgenomics.com/u/sevenbridges/simons-genome-diversity-project-sgdp/files) except for JHM06 which is described in [McColl H et al., 2018](https://science.sciencemag.org/content/361/6397/88?fbclid=IwAR2K_oBfFR7SUPCk7l8r4GBQO5XsIR01MJCxHz3vRKRlovQ2iQSgIuhVIo8). 

1. BAM -> VCF.

### 1. BAM to VCF

Convert BAM files to filtered VCF files.

##### Dependencies
1. Python (v3)
2. [SAMtools](http://www.htslib.org/download/) (v1.9).
3. [GenomeAnalysisToolkit](https://software.broadinstitute.org/gatk/) (v3.7).
4. [VCFtools](http://vcftools.sourceforge.net/perl_module.html) (v0.1.5).
5. filter_GC.py (provided in [scripts](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) folder).
6. filter_vcf_minDepth_maxDepth.py (provided in [scripts](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) folder).

##### Inputs and top level scripts
1. The [filter_modern_bams.sh](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) script.

##### The filter\_modern\_bams.sh script

You will need a machine which has at least 23 cores. BAM files are split by chromosome (autosomes only) and genotypes called using Genome Analysis Toolkit Unified Genotyper. Only sites within the **Gronau regions** will be called (see [README](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/blob/main/Gronau_filters/Gronau_filters.md) ). The input priors are set to equal for homozygous reference and homozygous alternate genotypes (as done for the Simons Genome Diversity panel). Positions where there is a G with a C called immediately up or downstream are removed as are positions with a C called and a G immediately up or downstream of it. C and G sites with a missing genotype immediately preceding or proceeding it are also removed. Genotypes are filtered by a minimum coverage of 20x and maximum depth defined as twice the average genotype coverage within the Gronau regions. Paths to specific tools/scripts, bed files containing the intervals of the Gronau regions and reference genome have to be manually changed within the filter\_modern\_bams.sh script 

**Running the script**

The script shoul be submitted within the folder containing the BAM file. 
```bash
bash filter_modern_bams.sh
```

## Modern and ancient samples

After VCF files were generated for both moderna dn ancient samples, they can be merged to create a final dataset

1. Merging VCF files.
2. Filtering final dataset.

### 1. Merging VCF files

VCF files are merged semi-manually. It is more efficient to merge per chromosome across samples (can be run in parallel) and then concatenate the chromosomes together (i.e. merge sample1_chr1 sample2_chr1 sample3_chr1  > all_samples_chr1. Then merge all_samples_chr1 all_samples_chr2 etc) than to merge across all chromosomes within a sample and then merge across samples. An example script "merge_vcf_example.sh" which merges files using VCFtools can be found in the [scripts](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/scripts) folder. If you already have vcf files per sample and you want to merge them across samples, you can use `bcftools merge -O z --threads 32 sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz > final_dataset.vcf.gz`.


### 2. Filtering final dataset

After merging, the dataset is filtered to remove all missing data and triallelic sites.

##### Removing missing data
All sites showing even just one missing genotype across the dataset are discarded using vcftools:

`vcftools --gzvcf final_dataset.vcf.gz" --max-missing 1 --recode-INFO-all --recode --stdout | bgzip -c -@ 32 > final_dataset_noMissing.vcf.gz`

##### Removing triallelic sites
Triallelic site are removed with bcftools:

`bcftools view final_dataset_noMissing.vcf.gz -O z -m1 -M2 --threads 32  -o final_dataset_noMissing_noTriallelic.vcf.gz`

##### Number of transitions and transversions
The number of transitions and transversions is calculate using bcftools:

`bcftools stats -s - final_dataset_noMissing_noTriallelic.vcf.gz > final_dataset_noMissing_noTriallelic_bcftools_stats`
