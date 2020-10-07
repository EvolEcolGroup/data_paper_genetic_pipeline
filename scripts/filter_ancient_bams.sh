#!/bin/bash

#####functions#####

parse_csv() {

	#get variables from csv file
	snp_list=$(awk -F ',' '{if ($1=="snp_list:") print $2}' $csv_file)
	reference_genome=$(awk -F ',' '{if ($1=="reference_genome:") print $2}' $csv_file)
	GATK=$(awk -F ',' '{if ($1=="GATK:") print $2}' $csv_file)
	filter_GC=$(awk -F ',' '{if ($1=="filter_GC.py:") print $2}' $csv_file)
	filter_vcf_minDepth_maxDepth=$(awk -F ',' '{if ($1=="filter_vcf_minDepth_maxDepth.py:") print $2}' $csv_file)
	vcf_concat=$(awk -F ',' '{if ($1=="vcf-concat:") print $2}' $csv_file)
	min_depth=$(awk -F ',' '{if ($1=="min_depth:") print $2}' $csv_file)
	max_depth=$(awk -F ',' '{if ($1=="max_depth:") print $2}' $csv_file)
	clean_up=$(awk -F ',' '{if ($1=="clean_up:") print $2}' $csv_file)
	merge=$(awk -F ',' '{if ($1=="merge:") print $2}' $csv_file)

}


check_bam_split() {

	#check that the sum of the records in the split bam files is equal to the number of records (chr1-22) in the original bam file.
	n_records_big_bam=$(samtools idxstats $file | egrep -w "chr[1-9]|chr1[0-9]|chr2[0-2]" | awk '{sum+=$3} END {print sum}')

	counter=0
	for i in {1..22}
	do
		n_records=$(samtools flagstat ${file_we}_chr${i}.bam | head -n1 | cut -f1 -d ' ')
		counter=$((counter+n_records)) 
	done 

	wait 

 	if [ "$counter" != "$n_records_big_bam" ]
 	then
 	 	echo "ERROR: BAM file splitting. Is the BAM file corrupted?"
 		exit 1
 	else
 		echo "split BAM file check: passed"
 	fi

}


check_DP_flag() {

	echo "checking the DP flag location..."
 	#check the DP field is in the expected location in the vcf file for all records
 	for i in {1..22}
 	do
		DP_check=$(grep -v '#' ${file_we}_Gronau_filters_chr${i}.vcf | grep -v '\./\.' | cut -f9 | cut -f3 -d":" | sort | uniq)
 		if [ "$DP_check" != "DP" ]
		then
 			echo "ERROR: vcf file ${file_we}_Gronau_filters_chr${i}.vcf is not in the expected format for depth filtering."
 			exit 1
 		fi
	done

 	echo "DP field in vcf file check: passed"
}

	
average_genotype_coverage() {
	#calculate the average depth of coverage across all Gronau sites
	echo "calculating average depth of coverage across all Gronau-filtered sites..."
	num_counter=0
	den_counter=0

	for i in {1..22}
	do

		grep -v '^#' ${file_we}_Gronau_filters_chr${i}.vcf  | awk '{if ($9=="GT:DP") {split($10, a,":"); print a[2]} else if  ($9=="GT:AD:DP:GQ:PL") {split($10, a,":"); print a[3]} else if ($9=="GT")  print "0"}' - |  sort | uniq -c |grep -v './.'| awk '{sum+=($1*$2)} END {print sum}' > ${file_we}_chr${i}_numcount.txt &
	done 

	wait

	for i in {1..22}
	do
		grep -v '^#' ${file_we}_Gronau_filters_chr${i}.vcf  | awk '{if ($9=="GT:DP") {split($10, a,":"); print a[2]} else if  ($9=="GT:AD:DP:GQ:PL") {split($10, a,":"); print a[3]} else if ($9=="GT")  print "0"}' - | sort | uniq -c | awk '{sum+=$1} END {print sum}' > ${file_we}_chr${i}_dencount.txt &
	done 

	wait

	num_counter=$(cat ${file_we}_chr*_numcount.txt | awk '{sum+=$1} END {print sum}')
	den_counter=$(cat ${file_we}_chr*_dencount.txt | awk '{sum+=$1} END {print sum}')

	average_genotype_coverage=$(awk 'BEGIN { rounded = sprintf("%.10f", "'${num_counter}'"/"'${den_counter}'"); print rounded }')

	echo "average depth of coverage across all Gronau sites is:" $average_genotype_coverage

	#tidy up intermediate files
	rm ${file_we}_chr*_{numcount,dencount}.txt
}


specify_min_max_depth() {

	if [[ ($min_depth = *"X"*) ]]
	then
		min_depth=$(echo $min_depth | sed s/X//g - | sed s/x//g -)
		min_depth=$(awk  'BEGIN { rounded = sprintf("%.0f", "'${min_depth}'"*"'${average_genotype_coverage}'"); print rounded }')
	fi

	if [[ ($max_depth = *"X"*) ]]
	then
		max_depth=$(echo $max_depth | sed s/X//g - | sed s/x//g -)
		max_depth=$(awk  'BEGIN { rounded = sprintf("%.0f", "'${max_depth}'"*"'${average_genotype_coverage}'"); print rounded }')
	fi

}


check_min_max_depth(){

	if [[ ($min_depth = 0) ]]
	then
		min_depth=1
	fi

	if (($max_depth < $min_depth))
	then
		echo "ERROR: minimum depth of coverage is greater than maximum depth of coverage"
		exit 1
	fi
}


clean_up() {

	if [[ $clean_up =~ ^(Y|YES|Yes|yes|1)$ ]]
	then
		rm $files
	fi
}


###################

csv_file=$1

#parse csv file
parse_csv

#if clean_up option is TRUE, remove intermediate files as we go along
if [[ $clean_up =~ ^(Y|YES|Yes|yes|1|T|TRUE|True)$ ]] 
then
	echo "file clean-up initialized..."
else   
	echo "file clean-up turned off..."
fi

#work on each bam file in the current working directory
for file in $(ls *bam)
do

	echo "working on file "$file

#####main code##### 

	#index file if necessary
	file_we=${file%.bam}
	samtools index $file

	if [ ! -f ${file_we}.bam.bai ] 
	then 
		echo "indexing file..."
	 	samtools index $file
	 else
	 	echo "index present, skipping indexing step..."
	 fi

	#split bam file up by chromosome (chromosomes 1-22)
	#all instances of "chr" are stripped from the file
	echo "splitting BAM file by chromosome..."
	for i in {1..22}
	do
		samtools view -h $file chr$i | samtools view -Sb - > ${file_we}_chr${i}.bam 2> /dev/null &
		#samtools view -h $file chr$i | sed s/chr//g - | samtools view -Sb - > ${file_we}_chr${i}.bam 2> /dev/null & 
	done 

	wait
	
	#check that the bam file split correctly
	check_bam_split

	#index split bam files
	echo "indexing split BAM files..."
	for i in {1..22}
	do
		samtools index ${file_we}_chr${i}.bam &
	done

	wait

	echo "calling genotypes..."
	#call genotypes using GATK Unified Genotyper
	#input priors equal for homozgous ref and homzygous alt (as in Simons Diversity panel)
	#minimum base quality of 20
	#emit all sites as we want homozygous ref calls 
	for i in {1..22}
	do
		java -jar -Xmx28g $GATK -T UnifiedGenotyper -I ${file_we}_chr${i}.bam -L ${snp_list}/autosomes_with_Gronau_filters_chr${i}.bed -R ${reference_genome} -out_mode EMIT_ALL_SITES  -mbq 20 -inputPrior 0.0010 -inputPrior 0.4995 -o ${file_we}_Gronau_filters_chr${i}.vcf > ${file_we}_Gronau_filters_chr${i}_GATK.warn 2> ${file_we}_Gronau_filters_chr${i}_GATK.log&
	done

	wait

	#merge GATK log files together 

	for i in {1..22}
	do
		echo "${file_we}_Gronau_filters_chr${i}_GATK.log" >> ${file_we}_Gronau_filters_GATK.log
		cat ${file_we}_Gronau_filters_chr${i}_GATK.log >> ${file_we}_Gronau_filters_GATK.log

		echo "${file_we}_Gronau_filters_chr${i}_GATK.warn" >> ${file_we}_Gronau_filters_GATK.log
		cat ${file_we}_Gronau_filters_chr${i}_GATK.warn >> ${file_we}_Gronau_filters_GATK.log

	done
	
	rm ${file_we}_Gronau_filters_chr*_GATK.log
	rm ${file_we}_Gronau_filters_chr*_GATK.warn

	#clean up files if option specified
	files=${file_we}_chr*.bam*
	clean_up

	#check_DP_flag

	#if filtering using the "times average coverage" option, calculate the average depth of coverage across all Gronau sites
	if [[ ($min_depth = *"X"*) || ( $max_depth = *"X"* )]]; then
		average_genotype_coverage	
	fi

 	#filter out GC sites (and indels and missing sites)
 	echo "filtering GC sites..."
 	for i in {1..22}
 	do
		python ${filter_GC} ${file_we}_Gronau_filters_chr${i}.vcf ${file_we}_Gronau_GC_filters_chr${i}.vcf &
 	done 

 	wait

	files=${file_we}_Gronau_filters_chr*.vcf*
 	clean_up

	#filter by minimum and maximum depth thresholds
	specify_min_max_depth
	check_min_max_depth

	echo "filtering sites by minimum depth of coverage $min_depth and maximum depth of coverage $max_depth" 
	for i in {1..22}
	do
		python ${filter_vcf_minDepth_maxDepth} ${file_we}_Gronau_GC_filters_chr${i}.vcf ${file_we}_Gronau_GC_filters_chr${i}_minDP${min_depth}_maxDP${max_depth}.vcf ${min_depth} ${max_depth} &
	done

	wait


	files=${file_we}_Gronau_GC_filters_chr[1-9].vcf
 	clean_up
 	files=${file_we}_Gronau_GC_filters_chr[1-2][0-9].vcf
 	clean_up

 	#merge files if option specified, bgzip file/files and create tabix index
 	if [[ $merge =~ ^(Y|YES|Yes|yes|1|T|TRUE|True)$ ]] 
	then
		echo "merging files..."
		concat_list=$(ls -v1 ${file_we}_Gronau_GC_filters_chr*_minDP${min_depth}_maxDP${max_depth}.vcf | tr '\n' ' ')

		${vcf_concat} $concat_list | bgzip -c > ${file_we}_Gronau_GC_filters_minDP${min_depth}_maxDP${max_depth}_merge.vcf.gz
		tabix -h ${file_we}_Gronau_GC_filters_minDP${min_depth}_maxDP${max_depth}_merge.vcf.gz
		files=${file_we}_Gronau_GC_filters_chr*_minDP${min_depth}_maxDP${max_depth}.vcf
		clean_up
		echo "files merged"
		echo "all done!"
	else   
		for i in {1..22}
		do
			bgzip ${file_we}_Gronau_GC_filters_chr${i}_minDP${min_depth}_maxDP${max_depth}.vcf &
		done

		wait

		for i in {1..22}
		do
			tabix -h ${file_we}_Gronau_GC_filters_chr${i}_minDP${min_depth}_maxDP${max_depth}.vcf.gz &
		done

		echo "all done!"
	fi

done
