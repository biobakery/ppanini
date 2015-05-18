import os
import sys
import pdb
import re
import argparse
import numpy
import logging


if __name__ == '__main__':
	foo = open(sys.argv[1])
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
	norm_dm = norm_dm/sum(norm_dm)
	for i, gene in enumerate(genes_order):
		#pdb.set_trace()
		print '\t'.join([gene]+[str(j) for j in list(norm_dm[i])])

