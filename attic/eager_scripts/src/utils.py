# Shafquat, Afrah
# shafquat@hsph.havard.edu
# March, 2014

import os
import sys
import re
import pdb

def write_dict(foo_dict, typ, tmp, data_folder):
	''' Utility function to write dicts to output folder
	    Input:
	    * foo_dict: Dict that needs to be written; can be (A) {key:[],...} or (B) {key:{nest_key:[]}}
	    * tmp: if 'hmgc_gi','hmgc_files', 'prev_hmgc' then treats as A, else as B
	    * data_folder: output folder path where files are written'''
	if typ == 'file':
        	foo = open(data_folder + '/' + tmp + '.txt', 'w')
		for i in foo_dict:
			foo.writelines([i + '\t' + str.join('\t', foo_dict[i]) + '\n'])
		foo.close()
	elif typ == 'folder':
		print tmp
		try:
			os.mkdir(data_folder+'/'+tmp)
		except:
			pass
		for fname in foo_dict:
			foo = open(data_folder+'/'+tmp+'/'+fname+'.txt', 'w')
			for key in foo_dict[fname]:
				foo.writelines([key + '\t' + str(foo_dict[fname][key]) + '\n'])
			foo.close()
	print 'Files written\t' + data_folder + '/' + tmp
