import os
import sys
import pdb
import re
import argparse
import numpy

'''Normalizes table given in first argument. Prints table as std.out'''
numpy.seterr(divide='ignore', invalid='ignore')
def main():
	help = ['-h', '--h', '--help']
	if sys.argv[1] in help:
		print 'Usage: python '+sys.argv[0]+' <input_table> > <normalized_table>'
		sys.exit()
		
	foo = open(sys.argv[1])
	multiplier = float(sys.argv[2])
	metadata = []
	genes_order = []
	dm = []
	for line in foo:
		if line.startswith('#'):
			print line.strip()
		else:
			split_i = line.split('\t')
			dm += [[float(i) for i in split_i[1:]]]
			genes_order += [split_i[0].strip()]
	norm_dm = numpy.array(dm)
	norm_dm = norm_dm*multiplier/sum(norm_dm)
	# pdb.set_trace()
	for i in range(len(norm_dm)):
		for j in range(len(norm_dm[i])):
			norm_dm[i][j] = numpy.nan_to_num(norm_dm[i][j])
	for i, gene in enumerate(genes_order):
		print '\t'.join([gene]+[str(j) for j in list(norm_dm[i])])

if __name__ == '__main__':
	main()
