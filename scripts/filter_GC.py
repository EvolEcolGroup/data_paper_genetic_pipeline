
"""
Eppie Jones
28th February 2018

Script to make a vcf file where GC/CG sites have been removed.
It will also filter out indels and sites where no alleles have been called. 
If there is missing information before or after a G or C, the G/C allele is removed (conservative approach). 

Running the script:
python filter_GC.py <input.vcf> <output.vcf>

"""
import sys
import argparse

parser = argparse.ArgumentParser(description='Filter GC sites from a vcf file. Indels and sites where no alleles have been called are also filtered out.')
parser.add_argument("input_vcf",help="input file in vcf format")
parser.add_argument("output_vcf",help="output file in vcf format")
parser.parse_args()

args=sys.argv
input_vcf=open(args[1],"r")
output_vcf=open(args[2],"w")

my_geno_list=[]
my_position_list=[]
my_line_list=[]

# start reading in input vcf file
for each in input_vcf:
	each=each.strip()

# write header lines to output file
	if each.startswith("#"):
		output_vcf.write(each + "\n")

# get genotype information from the vcf file (position, reference allele, alternate allele and 0/0, 0/1, 1/1 status) and save to lists
	if not each.startswith("#"):
		col=each.split()
		pos=col[1]
		ref_allele=col[3]
		alt_allele=col[4]
		geno=(col[9].split(":"))[0]

#appending position and line to lists after establishing genotype as a way to exclude missing genotypes
		if geno=="0/0":
			my_geno=ref_allele
			if len(my_geno)==1:
				my_geno_list.append(my_geno)
				my_position_list.append(pos)
				my_line_list.append(each)

		if geno=="1/1":
			my_geno=alt_allele
			if len(my_geno)==1:
				my_geno_list.append(my_geno)
				my_position_list.append(pos)
				my_line_list.append(each)

		if geno=="0/1":
			my_geno=ref_allele+alt_allele
			if len(my_geno)==2:
				my_geno_list.append(my_geno)
				my_position_list.append(pos)
				my_line_list.append(each)

# If the genotype contains a G and there is a C immediately upstream or downstream from it, ignore this datapoint. 
# If the genotype contains a C and there is a G immediately upstream or downstream from it, ignore this datapoint. 
for i in range(len(my_geno_list)):
	if i>0 and ("C" in my_geno_list[i]) and ("G" in my_geno_list[i-1]) and int(my_position_list[i])-int(my_position_list[i-1])==1:
		pass
	elif i<(len(my_geno_list)-1) and ("C" in my_geno_list[i]) and ("G" in my_geno_list[i+1]) and int(my_position_list[i+1])-int(my_position_list[i])==1:
		pass
	elif i>0 and ("G" in my_geno_list[i]) and ("C" in my_geno_list[i-1]) and int(my_position_list[i])-int(my_position_list[i-1])==1:
		pass
	elif i<(len(my_geno_list)-1) and ("G" in my_geno_list[i]) and ("C" in my_geno_list[i+1]) and int(my_position_list[i+1])-int(my_position_list[i])==1:
		pass

#remove C and G positions if they are preceded or proceeded by missing data
	elif i>0 and (("C" in my_geno_list[i]) or ("G" in my_geno_list[i])) and (int(my_position_list[i])-int(my_position_list[i-1])>1):
		pass
	elif i<(len(my_geno_list)-1) and (("C" in my_geno_list[i]) or ("G" in my_geno_list[i])) and int(my_position_list[i+1])-int(my_position_list[i])>1:
		pass

#write every other genotype to the output file 
	else:
		output_vcf.write(my_line_list[i] + "\n")

