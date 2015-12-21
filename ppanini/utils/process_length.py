import os
import sys
import re
import pdb
import matplotlib

from matplotlib import pyplot
from src import utilities

'''Tmp file to parse results'''

filename = sys.argv[1]
fasta_foo = utilities.read_fasta(filename)

fasta_len = {}
for gene in fasta_foo:
	# if len(fasta_foo[gene]) <100:
	fasta_len[gene] = len(fasta_foo[gene])

x = fasta_len.values()
x = sorted(x)

pyplot.plot(x,'r.')
pyplot.savefig('tmp_length.pdf')

