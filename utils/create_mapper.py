import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess


if __name__ == '__main__':
	labels = ['SAMS', 'BAMS', 'GFF3S', 'ABUNDANCE_TABLES', 'ANNOTATION', 'READS', 'CONTIG_ASSEMBLIES', 'FAAS', 'FNAS']
	parser = argparse.ArgumentParser()
	parser.add_argument('--assemblies', default=False, help='FOLDER containing ASSEMBLIES')
	parser.add_argument('--reads', default=False, help='Folder containing READS files', default=4)
	parser.add_argument('--gff3s', default=False, help='GFF3 Folder')
	parser.add_argument('--sams', default=False, help='SAMS Folder')
	parser.add_argument('--bams', default=False, help='SAMS Folder')
	parser.add_argument('--abund', default=False, help='SAMS Folder')
	parser.add_argument('--annot', default=False, help='SAMS Folder')
	parser.add_argument('--faas', default=False, help='SAMS Folder')
	parser.add_argument('--fnas', default=False, help='SAMS Folder')
	parser.add_argument('--niche',  default=False, help='GFF3 Folder')
	parser.add_argument('--samples',  default=False, help='GFF3 Folder', required=True)
	parser.add_argument('-o', '--output_table', help='Gene Table to write', default=sys.stdout)

	args = parser.parse_args()
	labels = {'SAMS':args.sams, \
			  'BAMS':args.bams, \
			  'NICHE':args.niche, \
			  'GFF3S':args.gff3s, \
			  'ABUNDANCE_TABLES':args.abund, \
			  # 'ANNOTATION':args.annot, \
			  'READS':args.reads, \
			  'CONTIG_ASSEMBLIES':args.assemblies, \
			  'FAAS':args.faas,\
			  'FNAS': args.fnas}
	samples = []
	with open(args.samples) as foo:
		for line in foo:
			samples +=[line.strip()]
	
	files = {}
	all_files = {}

	for sample in samples:
		files[sample] = {}
	keys = []
	for i in labels:
		if labels[i]:
			keys += [i]
			all_files[i] = os.listdir(labels[i])
			for fname in all_files[i]:
				sample = fname.rpartition('/')[0].split('.')[0].strip()
				files[sample][i] = labels[i]+'/'+fname
	# assemblies_folder = args.assemblies_folder
	# gff3s_folder = args.gff3s_folder
	# reads_folder = args.reads_folder
	# assemblies = os.listdir(assemblies_folder)
	#  #.scaffolds.fa
	# reads = os.listdir(reads_folder)
	# gff3s = os.listdir(gff3s_folder)
	# samples = [i.split('.')[0] for i in gff3s] #.tar.bz2 or .tar.gz
	# #.with_fasta.gff3
	annot_files = []
	if args.annot:
		annot_files = os.listdir(args.annot)
	line_end = '\n'
	meta_line = '\n'
	if args.niche:
		meta_line ='\tNICHE\n' 
		line_end = args.niche+'\n'
	with open(args.output_table, 'w') as foo:
		if 
		foo.writelines([str.join('\t', ['#SAMPLES']+keys)+meta_line])
		for sample in files:
			foo.writelines(['\t'.join([sample]+[files[sample][i] for i in keys])+line_end])
			# a = assemblies_folder+'/'+sample+'.scaffolds.fa'
			# print sample
			# print a
			# if not os.path.exists(a):
			# 	a = assemblies_folder+'/'+sample+'.scaffold.fa'
			# if not os.path.exists(a):
			# 	print a
			# 	raise IOError('CANT FIND FILE')

			# if sample+'.tar.gz' in reads:
			# 	r = reads_folder+'/'+ sample+'.tar.gz'
			# else:
			# 	r = reads_folder+'/'+ sample+'.tar.bz2'

			# g = gff3s_folder+'/'+sample+'.with_fasta.gff3'
			# foo.writelines([str.join('\t', [sample, a, r, g, niche])+'\n'])

