import os
import re
import sys
from print_table_functions import *
from hmgi_parser import *

if __name__ == '__main__':
	if not len(sys.argv) >= 4:
		print 'Usage %s <abundance file> <uniref annotation file> <outputfile> <optional: annotation folder> <if annot, then datafolder>' %sys.argv[0]
		raise IOError 
		
	short_udict = get_uniref_data(sys.argv[2]) #uniref file
	abund_dict = get_abund_data(sys.argv[1]) #abundance file
	print 'short and abund complete'
	check = True
	try:
		annot_folder = sys.argv[4]
	except:
		check = False

	if check:
		data_folder = sys.argv[5] #'/n/huttenhower_lab_nobackup/data/hmp/stool_annotation/data_files'
		gi_annot = get_gi_annot(annot_folder)
		print 'gi_annot complete'
		write_dict(gi_annot, 'folder','gi_annot', data_folder)
		annot = get_annotations(short_udict, gi_annot)
		print 'annot complete'
		write_annotated_table(abund_dict, short_udict, sys.argv[3], annot)
	else:
		write_annotated_table(abund_dict, short_udict, sys.argv[3], [])
