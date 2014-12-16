import os
import re
import sys
import pdb
import argparse
import numpy

def read_contig_gid(gff3_folder):
	foo_list = os.listdir(gff3_folder)
	gid_contig = {}
	for fname in foo_list:
		sample = fname.split('.')[0].strip()
		gid_contig[sample]={}
		foo = open(gff3_folder+'/'+fname)
		lines = foo.readlines()
		for line in lines:
			if not line[0] == '#':
				if line.split('\t')[2] == 'CDS':
					split_i = line.split('\t')
					if not split_i[8].split('=')[1].split(';')[0].strip() in gid_contig[sample]:
						gid_contig[sample][split_i[8].split('=')[1].split(';')[0].strip()] = [split_i[0].strip(), float(split_i[4])-float(split_i[3])+1]
	return gid_contig

def read_contig_abund(abund_folder):
	foo_list = os.listdir(abund_folder)
	contig_abund = {}
	for fname in foo_list:
		sample = fname.split('.')[0].strip()
		contig_abund[sample]={}
		foo = open(abund_folder+'/'+fname)
		lines = foo.readlines()
		for line in lines:
			split_i = line.split('\t')
			contig_abund[sample][split_i[0].strip()]=float(re.sub('[\t\r\n]', '', split_i[2]))/(float(split_i[1])-100)
	
	return contig_abund

def read_gid_uniref(annot_folder, gid_contig):
	foo_list = os.listdir(annot_folder)
	gid_annot = {}
	for fname in foo_list:
		sample = fname.split('.')[0].strip()
		gid_annot[sample]={}
		foo = open(abund_folder+'/'+fname)
		lines = foo.readlines()
		for line in lines:
			if not '#' == line[0]:
				split_i = line.split('\t')
				id = float(split_i[2]) * float(split_i[3])/(gid_contig[sample][split_i[0]][1]/3.0)
				if not split_i[0] in gid_annot[sample]:
					if id >= 90.0:
						gid_annot[sample][split_i[0]] = [split_i[1], id]
				elif id >= gid_annot[sample][split_i[0]][1]:
						gid_annot[sample][split_i[0]] = [split_i[1], id]
	return gid_annot

def read_gid_uniref_list(annot_file, gid_contig):
	foo = open(annot_file)
	lines = foo.readlines()
	
	gid_annot = {}
	for line in lines:
		if not '#' == line[0]:
			sample = line.split('.')[0]
			if sample in gid_contig:
				if not sample in gid_annot:
					gid_annot[sample] = {}
				split_i = line.split('\t')
				id = float(split_i[2]) * float(split_i[3])/(gid_contig[sample][split_i[0]][1]/3.0)
				if not split_i[0] in gid_annot[sample]:
					if id >= 90.0:
						gid_annot[sample][split_i[0]] = [split_i[1], id]
				elif id >= gid_annot[sample][split_i[0]][1]:
						gid_annot[sample][split_i[0]] = [split_i[1], id]
	return gid_annot


def put_together(gid_contig, contig_abund, gid_annot, out_file):
	gid_contig = gid_contig # Samples: genes ID: [contig ID, length of gene]
	contig_abund = contig_abund # Samples: contig: FPKM
	gid_annot = gid_annot #Samples: genes ID: [annotation, %id]
	foo = open(out_file,'w')
	samples_list = gid_contig.keys()
	foo.writelines([str.join('\t', ['Gene ID', 'Annotation']+samples_list)+'\n'])
	for sample in gid_contig:
		for gid in gid_contig[sample]:
			gid_abund = list(numpy.zeros(len(samples_list))) # make a zeros row
			gid_abund[samples_list.index(sample)] = contig_abund[sample][gid_contig[sample][gid][0]] # add abundance of gene in sample
			str_gid_abund = [str(i) for i in gid_abund] #convert the gid_abund to list of strings
			if sample in gid_annot and gid in gid_annot[sample]:
					annot = gid_annot[sample][gid][0]
			else:
					annot = 'UniRef_unknown'
			foo.writelines([str.join('\t', [gid]+[annot]+str_gid_abund)+'\n'])
	foo.close()



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-r','--reads_abundance', help='Folder containing samtools results')
	parser.add_argument('-m','--contig_gene_mapping', help='Folder containing GFF3 files')
	parser.add_argument('-a','--gene_annotation', nargs = '?' , help='List of UniRef annotation results')
	parser.add_argument('-o','--output_file', help='Output file')
	args = parser.parse_args()

	reads_abundance = args.reads_abundance
	contig_gene_mapping = args.contig_gene_mapping
	gene_annotation = args.gene_annotation

	gid_contig = read_contig_gid(contig_gene_mapping)
	contig_abund = read_contig_abund(reads_abundance)
	gid_annot = read_gid_uniref_list(gene_annotation, gid_contig)

	put_together(gid_contig, contig_abund, gid_annot, args.output_file)


	