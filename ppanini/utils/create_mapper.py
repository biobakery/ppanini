import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess

'''Creates a mapper table for PREPPANINI input'''
def main():
	labels = ['SAMS', 'BAMS', 'GFF3S', 'ABUNDANCE_TABLES', 'ANNOTATION', 'READS', 'CONTIG_ASSEMBLIES', 'FAAS', 'FNAS']
	parser = argparse.ArgumentParser()
	parser.add_argument('--assemblies', default=False, help='FOLDER containing ASSEMBLIES')
	parser.add_argument('--reads', default=False, help='Folder containing READS files')
	parser.add_argument('--gff3s', default=False, help='GFF3 Folder')
	parser.add_argument('--sams', default=False, help='SAMS Folder')
	parser.add_argument('--bams', default=False, help='SAMS Folder')
	parser.add_argument('--abund', default=False, help='SAMS Folder')
	parser.add_argument('--annot', default=False, help='SAMS Folder')
	parser.add_argument('--faas', default=False, help='SAMS Folder')
	parser.add_argument('--fnas', default=False, help='SAMS Folder')
	parser.add_argument('--niche',  default=False, help='GFF3 Folder')
	parser.add_argument('--samples',  default=False, help='GFF3 Folder', required=True)
	parser.add_argument('-o', '--output-table', dest= 'output_table', help='Gene Table to write', default=sys.stdout)

	args = parser.parse_args()
	labels = {'SAMS':args.sams, \
			  'BAMS':args.bams, \
			  'GFF3S':args.gff3s, \
			  'ABUNDANCE_TABLES':args.abund, \
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
				sample = fname.rpartition('/')[-1].split('.')[0].strip()
				try:
					files[sample][i] = labels[i]+'/'+fname
				except:
					print sample
					pass

	annot_files = []
	if args.annot:
		x = os.listdir(args.annot)
		annot_files = [args.annot+'/'+i for i in x]
	line_end = '\n'
	meta_line = '\n'
	if args.niche:
		meta_line ='\tNICHE\n' 
		line_end = '\t'+args.niche+'\n'
	with open(args.output_table, 'w') as foo:
		if annot_files:
			print 'comes here'
			meta_line = '\tANNOTATION'+meta_line
		count = 0
		foo.writelines([str.join('\t', ['#SAMPLES']+keys)+meta_line])
		for sample in files:
			if annot_files:
				if count < len(annot_files):
					foo.writelines(['\t'.join([sample]+[files[sample][i] for i in keys])+'\t'+annot_files[count]+line_end])
				else:
					foo.writelines(['\t'.join([sample]+[files[sample][i] for i in keys])+'\t'+line_end])
				count +=1
			else:
				foo.writelines(['\t'.join([sample]+[files[sample][i] for i in keys])+line_end])


if __name__ == '__main__':
	main()