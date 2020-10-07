"""
Eppie Jones
5th March 2018

filter_vcf_minDepth_maxDepth.py

Use program to set a minimum and maximum depth of coverage for genotype calls.

Running the script:

python filter_vcf_minDepth_maxDepth.py <input_file> <output_file> <minimum depth of coverage> <maximum depth of coverage>

"""

import sys
import os

args = sys.argv
input_file = open(args[1],'r')
output_file = open(args[2],'w')
mindepth = args[3]
maxdepth = args[4]

for each in input_file:
	each = each.strip()
	if each.startswith('#'):
		output_file.write(each + '\n')

	else:
		col = each.split()
		if col[5] != '.' and './.' not in each:


			format_field=col[8]
			format_field=format_field.split(":")
			
			info = col[9]
			info = info.split(':')

			if "AD" not in format_field:
				DP_location=format_field.index("DP")
				true_depth=info[DP_location]
				if int(true_depth) >= int(mindepth) and int(true_depth) <= int(maxdepth):
					output_file.write(each + '\n')		

			
			if "AD" in format_field:
				depth_list=[]

				DP_location=format_field.index("DP")
				DP = info[DP_location]
				depth_list.append(int(DP))
				
				AD_location=format_field.index("AD")
				AD = info[AD_location]
				AD = AD.split(',')
				AD = int(AD[0]) + int(AD[1])
				depth_list.append(AD)

				true_depth= min(depth_list)

				if int(true_depth) >= int(mindepth) and int(true_depth) <= int(maxdepth):
					output_file.write(each + '\n')				
input_file.close()
output_file.close()
