
import os
import sys
import pdb
import re

'''Tmp file to parse results'''

foo = open(sys.argv[1])

genes_genomes = {}
for line in foo:
	split_i = line.split('\t')
	for i in split_i[1:]:
		sp = re.sub('[\r\t\n]','',i).strip()
		if '|' in i:
			sp = [j for j in i.split('|') if 'g__' in j and '.s__' in j][0]
		try:
			genes_genomes[split_i[0]] += [sp]
		except:
			genes_genomes[split_i[0]] = [sp]
# pdb.set_trace()
for gene in genes_genomes:
	genes_genomes[gene] = list(set(genes_genomes[gene]))
	print '\t'.join([gene]+genes_genomes[gene])
