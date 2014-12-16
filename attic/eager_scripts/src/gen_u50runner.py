import os
import pdb
import re
import sys


def check_equals(fasta_inds, u90_inds):
	if not len(u90_inds) == len(fasta_inds):
		return False
	else:
		return True

def write_u50(fasta_file, fasta_inds, u90_inds, u50_out):
	
	fasta_inds = open(fasta_inds)
	fasta_inds = fasta_inds.readlines()

	u90_inds = open(u90_inds)
	u90_inds = [re.sub('[\r\t\n]','', i) for i in u90_inds.readlines() if '#' not in i]
	print len(fasta_inds)
	print len(u90_inds)
	#pdb.set_trace()	
	check = len(fasta_inds) == len(u90_inds)

	if not check:
		fasta_dict = {}
		inds_all = []
		for i in fasta_inds:
			split_i = i.split(':')
			fasta_dict[re.sub('[\r\n\t]','',split_i[1][1:]).strip()] = int(split_i[0])-1
			inds_all += [int(split_i[0])-1]
		abs_fasta = []
		for i in fasta_dict:
			if not i in u90_inds:
				abs_fasta += [fasta_dict[i]]
		abs_fasta.sort()
		fasta_file = open(fasta_file)
		fasta_file = fasta_file.readlines()
		u50_out = open(u50_out, 'w')
		for i in abs_fasta:
			start = i
			try:	
				ind = inds_all.index(i)+1
				end = inds_all[ind]
			except:
				end = -1
			u50_out.writelines(fasta_file[i:end])
		u50_out.close()
	else:
		print 'No file written'
		return None
		
if __name__ == '__main__':
	fasta_file = sys.argv[1]
	fasta_inds = sys.argv[2]
	u90_inds = sys.argv[3]
	u50_out = sys.argv[4]

	write_u50(fasta_file, fasta_inds, u90_inds, u50_out) 
