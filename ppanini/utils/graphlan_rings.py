import os
import sys


'''Tmp file to parse results'''

foo = open(sys.argv[1])
genomes = {}
for line in foo:
	split_line = line.split('\t')
	genomes[split_line[0]] = float(split_line[1].strip())
max_val = float(max(genomes.values()))

for genome in genomes:
	genomes[genome] = genomes[genome]*5.0/max_val

with open(sys.argv[1]+'_rings.txt','w') as foo_mod:
	foo_mod.writelines(['ring_internal_separator_thickness\t1\t1.0\nring_width\t5\t0.5\n'])
	for genome in genomes:
		foo_mod.writelines([genome+'\tring_height\t1\t'+str(genomes[genome])+'\n'])

