#!/bin/bash

for file in `ls *bam`
do
	file_we=${file%.bam}
	samtools index $file

	#split bam file up by chromosome
	for i in {1..22}
	do
		samtools view -bh $file $i > ${file_we}_chr${i}.bam & 
	done 
	wait

	#index bam files
	for i in {1..22}
	do
		samtools index ${file_we}_chr${i}.bam &
	done
	wait

	#call genotypes using GATK Unified Genotyper
	#input priors equal for homozgous ref and homzygous alt (as in Simons Diversity panel)
	#minimum base quality of 20
	for i in {1..22}
	do
		java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -I ${file_we}_chr${i}.bam -L ~/SGDP/Gronau_filters/Gronau_regions_to_keep_per_chromosome/autosomes_with_Gronau_filters_chr${i}.bed -R /home/eppie/hs37d5/hs37d5_nochr/hs37d5.fa -out_mode EMIT_ALL_SITES  -mbq 20 -inputPrior 0.0010 -inputPrior 0.4995 -o ${file_we}_Gronau_filters_chr${i}.vcf &
	done
	wait

	#filter out GC sites (and indels and missing sites)
	for i in {1..22}
	do
		python filter_GC.py ${file_we}_Gronau_filters_chr${i}.vcf ${file_we}_Gronau_GC_filters_chr${i}.vcf &
	done 
	wait

	#concatenate vcf files
	concat_list=$(ls -v1 *_Gronau_GC_filters_chr*vcf)
	vcf-concat $concat_list  > ${file_we}_Gronau_GC_filters_merge.vcf

	#filter by minimum depth of coverage 20 and maximum depth of coverage of (2 x average genotype coverage)
	num=$(grep -v '^#' ${file_we}_Gronau_GC_filters_merge.vcf  | cut -f10 | cut -f3 -d':' | sort | uniq -c | awk '{sum+=($1*$2)} END {print sum}')
	den=$(grep -v '^#' ${file_we}_Gronau_GC_filters_merge.vcf  | cut -f10 | cut -f3 -d':' | sort | uniq -c | awk '{sum+=$1} END {print sum}')
	twice_average_genotype_coverage=$(awk  'BEGIN { rounded = sprintf("%.0f", "'${num}'"/"'${den}'"); print rounded*2 }')
	python filter_vcf_minDepth_maxDepth.py ${file_we}_Gronau_GC_filters_merge.vcf ${file_we}_Gronau_GC_filters_minDP20_maxDP${twice_average_genotype_coverage}.vcf 20 ${twice_average_genotype_coverage}

	bgzip ${file_we}_Gronau_GC_filters_merge.vcf
	bgzip ${file_we}_Gronau_GC_filters_minDP20_maxDP${twice_average_genotype_coverage}.vcf
	rm *chr*
done
