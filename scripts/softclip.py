import sys
import os
import fileinput

for each in fileinput.input():
	each = each.strip('\n')
	if each.startswith('@'):
		print each
	else:
		column = each.split()
		head_columns = '\t'.join(column[0:10])
		tail_columns = '\t'.join(column[11:])

		quality_string = column[10]
		length_quality_string = len(quality_string)
		soft_clipped_quality_string = '##' + quality_string[2:(length_quality_string - 2)] + '##'

		print head_columns + '\t' + soft_clipped_quality_string + '\t' + tail_columns 
