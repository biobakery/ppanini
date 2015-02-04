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
	samples = [i.split('.')[0] for i in assemblies] #.tar.bz2 or .tar.gz
	#.with_fasta.gff3

	with open(args.output_table, 'w') as foo:
		foo.writlines([str.join('\t', ['SAMPLES','CONTIG_ASSEMBLIES','READS', 'GFF3S', args.niche])+'\n'])
		for sample in samples:
			a = assemblies_folder+'/'+sample+'.scaffolds.fa'
			if sample+'.tar.gz' in reads:
				r = reads_folder+'/'+ sample+'.tar.gz'
			else:
				r = reads_folder+'/'+ sample+'.tar.bz2'
			g = gff3s_folder+'/'+sample+'.with_fasta.gff3'
			foo.writlines([str.join('\t', [sample, a, r, g])+'\n'])

