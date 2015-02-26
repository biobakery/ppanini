import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--assemblies_folder', help='FOLDER containing ASSEMBLIES')
	parser.add_argument('-r', '--reads_folder', help='Folder containing READS files', default=4)
	parser.add_argument('-g', '--gff3s_folder', help='GFF3 Folder')
	parser.add_argument('-n', '--niche', help='GFF3 Folder')
	parser.add_argument('-o', '--output_table', help='Gene Table to write', default=sys.stdout)

	args = parser.parse_args()
	assemblies_folder = args.assemblies_folder
	gff3s_folder = args.gff3s_folder
	reads_folder = args.reads_folder
	assemblies = os.listdir(assemblies_folder)
	 #.scaffolds.fa
	reads = os.listdir(reads_folder)
	gff3s = os.listdir(gff3s_folder)
	samples = [i.split('.')[0] for i in gff3s] #.tar.bz2 or .tar.gz
	#.with_fasta.gff3
	niche = args.niche
	with open(args.output_table, 'w') as foo:
		foo.writelines([str.join('\t', ['#SAMPLES','CONTIG_ASSEMBLIES','READS', 'GFF3S', 'NICHE'])+'\n'])
		for sample in samples:
			a = assemblies_folder+'/'+sample+'.scaffolds.fa'
			print sample
			print a
			if not os.path.exists(a):
				a = assemblies_folder+'/'+sample+'.scaffold.fa'
			if not os.path.exists(a):
				print a
				raise IOError('CANT FIND FILE')

			if sample+'.tar.gz' in reads:
				r = reads_folder+'/'+ sample+'.tar.gz'
			else:
				r = reads_folder+'/'+ sample+'.tar.bz2'

			g = gff3s_folder+'/'+sample+'.with_fasta.gff3'
			foo.writelines([str.join('\t', [sample, a, r, g, niche])+'\n'])

