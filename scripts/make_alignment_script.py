"""
18th May 2018
Eppie Jones

Script to make bash alignment pipeline.
A file containing file information and another containing parameters and links is required.

Running the script:

python make_aligment_script.py <file_info.csv> <parameters_and_links.csv> <output_file.sh>

"""

import sys
import os

args=sys.argv
input_csv=open(args[1],'r')
par_csv=open(args[2],'r')
output_sh=open(args[3],'w')

sample_ids=[]
fastq_files=[]
barcodes=[]
library=[]
bam_files=[]

#parse parameter file
lines=par_csv.readlines()
lines = map(lambda s: s.strip(), lines) #get rid of new line character
reference,picard,gatk,softclip,cutadapt,bwa,mq=(lines[1].split(','))[1],(lines[2].split(','))[1],(lines[3].split(','))[1],(lines[4].split(','))[1],(lines[7].split(','))[1],(lines[8].split(','))[1],(lines[9].split(','))[1]

#parse sample info file
for each in input_csv:
	if each.startswith('sample_id'):
		pass
	else:
		each = each.strip()
		col = each.split(',')
		sample_ids.append(col[0])
		fastq_files.append(col[1])
		barcodes.append(col[2])
		library.append(col[3])

#get a list of unique sample ids
unique_samples = list(set(sample_ids))

#get list of unmerged output bam files
out_bams=[]
for fastq in fastq_files:
		if fastq.endswith(".fastq.gz"):
			fastq=fastq[:-9]
                elif fastq.endswith(".fq.gz"):
                        fastq=fastq[:-6]
		else:
			print "ERROR:fastq files do not have the expected file extenstion (fastq.gz ot fq.gz)"
			sys.exit()			

		out_bams.append(fastq + "_sort.bam")
#initialize bash script
output_sh.write('#!/bin/bash\n\n')

#add column names to output stat files
output_sh.write('echo -e \'\t\' "sample_id" "file" "nreads" "nalign" >  alignment_stats_fastq_to_bam.txt\n')
output_sh.write('echo -e \'\t\' "sample_id" "nmerge" "nMQ" "nrmdup"  >  alignment_stats_merged_bam.txt\n\n')

#write first bash loop to trim and align fastq files
output_sh.write('sample_id=(' + ' '.join(sample_ids) + ')\n')
output_sh.write('fastq_file=(' + ' '.join(fastq_files) + ')\n')
output_sh.write('barcode=(' + ' '.join(barcodes) + ')\n')
output_sh.write('library=(' + ' '.join(library) + ')\n\n')


output_sh.write('for file in "${fastq_file[@]}";do\n')
output_sh.write('file_we=${file%.fastq.gz}\n')
output_sh.write('file_we=${file_we%.fq.gz}\n')
output_sh.write('cutadapt ' + cutadapt + ' ${file} -o ${file_we}_cutadapt_trimmed.fq.gz > ${file_we}_cutadapt_trimmed.log &\n')
output_sh.write('done\n\n')
output_sh.write('wait\n\n')

output_sh.write('for file in "${fastq_file[@]}";do\n')
output_sh.write('file_we=${file%.fastq.gz}\n')
output_sh.write('file_we=${file_we%.fq.gz}\n')
output_sh.write('bwa aln ' + bwa + ' ' + reference + ' ${file_we}_cutadapt_trimmed.fq.gz > ${file_we}.sai 2> ${file_we}_bwa_aln.log\n')
output_sh.write('done\n\n')
output_sh.write('wait\n\n')

output_sh.write('COUNTER=0\n')
output_sh.write('for file in "${fastq_file[@]}";do\n')
output_sh.write('file_we=${file%.fastq.gz}\n')
output_sh.write('file_we=${file_we%.fq.gz}\n')
output_sh.write('bwa samse -r "@RG\\tID:${file_we}\\tSM:${sample_id[COUNTER]}\\tLB:${barcode[COUNTER]}_${library[COUNTER]}\\tPL:ILLUMINA" ' + reference + ' ${file_we}.sai ${file_we}_cutadapt_trimmed.fq.gz | samtools view -Sb -F 4 - > ${file_we}.bam &\n')
output_sh.write('let COUNTER=COUNTER+1\n')
output_sh.write('done\n\n')
output_sh.write('wait\n\n')

output_sh.write('for file in "${fastq_file[@]}";do\n')
output_sh.write('file_we=${file%.fastq.gz}\n')
output_sh.write('file_we=${file_we%.fq.gz}\n')
output_sh.write('samtools sort ${file_we}.bam ${file_we}_sort &\n')
output_sh.write('done\n\n')
output_sh.write('wait\n\n')

output_sh.write('for file in "${fastq_file[@]}";do\n')
output_sh.write('file_we=${file%.fastq.gz}\n')
output_sh.write('file_we=${file_we%.fq.gz}\n')
output_sh.write('samtools index ${file_we}_sort.bam &\n')
output_sh.write('done\n\n')

#output some stats and tidy up
output_sh.write('COUNTER=0\n')
output_sh.write('for file in "${fastq_file[@]}";do\n')
output_sh.write('file_we=${file%.fastq.gz}\n')
output_sh.write('file_we=${file_we%.fq.gz}\n')
output_sh.write('nreads=$(grep "Total reads processed" ${file_we}_cutadapt_trimmed.log | awk \'{print $4}\' -)\n')
output_sh.write('nalign=$(samtools flagstat ${file_we}_sort.bam | head -n1 | cut -f1 -d \' \')\n')
output_sh.write('echo -e \'\t\' ${sample_id[COUNTER]} $file $nreads $nalign >>  alignment_stats_fastq_to_bam.txt\n')
output_sh.write('rm ${file_we}.sai ${file_we}_cutadapt_trimmed.fq.gz ${file_we}.bam\n')
output_sh.write('let COUNTER=COUNTER+1\n') 
output_sh.write('done\n\n')

output_sh.write('wait\n\n')

#assign bam files to sample ids
for each in unique_samples:
	file_list=[]
	for i in out_bams:
		sample_id=(i.split("_"))[0]
    		if  each == sample_id:
			file_list.append("I=" + i)
		file_string = ' '.join(file_list)
	bam_files.append('"' + file_string + '"')

#write second bash loop to merge bam files from the same sample and then to filter by mapping quality, remove duplicates, realign around indels and softclip 2bp from the ends of reads
output_sh.write('COUNTER=0\n')
output_sh.write('unique_samples=(' + ' '.join(unique_samples) + ')\n')
output_sh.write('bams=(' + ' '.join(bam_files) + ')\n')

output_sh.write('for sample_id in "${unique_samples[@]}";do\n')
output_sh.write('n_bams=$(echo ${bams[COUNTER]} | tr " " "\\n" | grep -c "bam")\n')
output_sh.write('bam_remove=$(echo ${bams[COUNTER]} | sed s/I=// -)\n')

output_sh.write('if [[ ($n_bams > 1) ]];then\n')
output_sh.write('java -jar ' + picard + ' MergeSamFiles VALIDATION_STRINGENCY=SILENT ${bams[COUNTER]} O=${sample_id}_merge.bam\n')
output_sh.write('samtools view -bhq' + mq + ' ${sample_id}_merge.bam > ${sample_id}_merge_mq' + mq +'.bam\n')
output_sh.write('samtools rmdup -s ${sample_id}_merge_mq' + mq +'.bam ${sample_id}_merge_mq' + mq +'_rmdup.bam 2> ${sample_id}_merge_mq' + mq +'_rmdup.log\n')
output_sh.write('samtools index ${sample_id}_merge_mq' + mq +'_rmdup.bam\n')
output_sh.write('java -Xmx100g -jar ' + gatk + ' -T RealignerTargetCreator -R ' + reference + ' -I ${sample_id}_merge_mq' + mq +'_rmdup.bam -o ${sample_id}_merge_mq' + mq +'_rmdup_forIndelRealigner.intervals\n')
output_sh.write('java -Xmx100g -jar ' + gatk + ' -T IndelRealigner -R ' + reference + ' -I ${sample_id}_merge_mq' + mq +'_rmdup.bam -targetIntervals ${sample_id}_merge_mq' + mq +'_rmdup_forIndelRealigner.intervals -o ${sample_id}_merge_mq' + mq +'_rmdup_IR.bam\n')
output_sh.write('samtools index ${sample_id}_merge_mq' + mq +'_rmdup_IR.bam\n')
output_sh.write('samtools view -h ${sample_id}_merge_mq' + mq +'_rmdup_IR.bam | python ' + softclip + ' | samtools view -Sbh - > ${sample_id}_merge_mq' + mq +'_rmdup_IR_2bp_softclip.bam\n')
#output some stats
output_sh.write('nmerge=$(samtools flagstat ${sample_id}_merge.bam | head -n1 | cut -f1 -d \' \')\n')
output_sh.write('nMQ=$(samtools flagstat ${sample_id}_merge_mq' + mq +'.bam | head -n1 | cut -f1 -d \' \')\n')
output_sh.write('nrmdup=$(samtools flagstat ${sample_id}_merge_mq' + mq +'_rmdup.bam | head -n1 | cut -f1 -d \' \')\n')
output_sh.write('echo -e \'\t\' $sample_id $nmerge $nMQ $nrmdup  >>  alignment_stats_merged_bam.txt\n')
#output_sh.write('rm $bam_remove \n')
output_sh.write('fi\n')


output_sh.write('if [[ ($n_bams == 1) ]];then\n')
output_sh.write('bam=$(echo ${bams[COUNTER]} | sed s/I=// -)\n')
output_sh.write('samtools view -bhq' + mq + ' $bam > ${sample_id}_mq' + mq +'.bam\n')
output_sh.write('samtools rmdup -s ${sample_id}_mq' + mq +'.bam ${sample_id}_mq' + mq +'_rmdup.bam 2> ${sample_id}_mq' + mq +'_rmdup.log\n')
output_sh.write('samtools index ${sample_id}_mq' + mq +'_rmdup.bam\n')
output_sh.write('java -Xmx100g -jar ' + gatk + ' -T RealignerTargetCreator -R ' + reference + ' -I ${sample_id}_mq' + mq +'_rmdup.bam -o ${sample_id}_mq' + mq +'_rmdup_forIndelRealigner.intervals\n')
output_sh.write('java -Xmx100g -jar ' + gatk + ' -T IndelRealigner -R ' + reference + ' -I ${sample_id}_mq' + mq +'_rmdup.bam -targetIntervals ${sample_id}_mq' + mq +'_rmdup_forIndelRealigner.intervals -o ${sample_id}_mq' + mq +'_rmdup_IR.bam\n')
output_sh.write('samtools index ${sample_id}_mq' + mq +'_rmdup_IR.bam\n')
output_sh.write('samtools view -h ${sample_id}_mq' + mq +'_rmdup_IR.bam | python ' + softclip + ' | samtools view -Sbh - > ${sample_id}_mq' + mq +'_rmdup_IR_2bp_softclip.bam\n')

#output some stats
output_sh.write('nmerge=$(samtools flagstat ${bam} | head -n1 | cut -f1 -d \' \')\n')
output_sh.write('nMQ=$(samtools flagstat ${sample_id}_mq' + mq +'.bam | head -n1 | cut -f1 -d \' \')\n')
output_sh.write('nrmdup=$(samtools flagstat ${sample_id}_mq' + mq +'_rmdup.bam | head -n1 | cut -f1 -d \' \')\n')
output_sh.write('echo -e \'\t\' $sample_id $nmerge $nMQ $nrmdup >>  alignment_stats_merged_bam.txt\n')
#output_sh.write('rm $bam \n')
output_sh.write('fi\n')
output_sh.write('let COUNTER=COUNTER+1\n')
output_sh.write('done\n')

