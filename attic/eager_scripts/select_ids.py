import os
import sys
import re
import pdb

if __name__ == '__main__':
	#python select_ids.py prevalent_abundant_ids.txt indices_headers.txt output_file_name.txt
	foo2 = open(sys.argv[1])#'pabs.txt')
	foo2 = foo2.readlines()

	foo = open(sys.argv[2]) #'tmp_inds.txt')
	foo = foo.readlines()

	ind_dict = {}
	for i in foo:
		ind_dict[re.sub('[\n\t\r]','',i.split('>')[1])] = int(i.split(':')[0])
		#ind_dict[re.sub('[\n\t\r]','',i.split('|')[-1])] = int(i.split(':')[0]) #anares
	foo3 = open(sys.argv[3],'w')
	for i in foo2:
		try:
			foo3.writelines([str(ind_dict[re.sub('[\r\t\n]','',i)])+'\n'])
			#foo3.writelines([str(ind_dict[re.sub('[\r\t\n]','',i.split('|')[-1])])+'\n'])
		except:
			pdb.set_trace()
	foo3.close()
