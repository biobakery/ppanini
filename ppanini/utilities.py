import os
import re
import sys
import pdb
import Bio
import numpy
import logging
import argparse
import subprocess
import multiprocessing

from Bio import Seq

logger = logging.getLogger(__name__)
def read_ppanini_imp_genes_table(filename):
	keys = { 'abundance_rank':0, 'prevalence_rank':0}
	abund = []
	prev = []
	genes = []
	ppanini_score = []
	
	with open(filename) as foo:
		for line in foo:
			split_line = [re.sub('[\r\t\n]','', i) for i in line.split('\t')]
			if line.startswith('#'):
				for i, val in enumerate(split_line):
					if 'abundance_rank' in val:
						keys['abundance_rank'] = i
					elif 'prevalence_rank' in val:
						keys['prevalence_rank'] = i
					#elif 'beta' in val:
					#	keys['beta'] = i
					elif 'ppanini_score' in val:
						keys['ppanini_score'] = i	
			else:
				# pdb.set_trace()
				#print line
				genes +=[split_line[0]]
				ppanini_score +=[None]
				abund +=[float(split_line[keys['abundance_rank']])]
				prev +=[float(split_line[keys['prevalence_rank']])]
	#print prev
	ppanini_table = {'genes': genes, 'abundance_rank': abund, 'prevalence_rank': prev, 'ppanini_score':ppanini_score }
	return ppanini_table
def read_fasta(fasta_filename):
	'''Reads a fasta_file and returns a fasta dict
	Input: fasta_filename = path_to_fasta_file
	Output: fasta_seq = {sequence_header: sequence, ...}'''
	logger.debug('read_fasta '+fasta_filename)

	fasta_file = open(fasta_filename)
	fasta_seq = {}
	name = ''

	for line in fasta_file:
		if not line.startswith('#'):
			if line.startswith('>'):
				name = re.sub('>','', line.split(' ')[0].strip())
			else:
				if name not in fasta_seq:
					fasta_seq[name] =  re.sub('[\r\t\n]','', line)
				else:
					fasta_seq[name] +=  re.sub('[\r\t\n]','', line)
	return fasta_seq
def parse_table(m8_filename, fasta_filename):
	'''Parse the BLAST results to give gene hits to genomes
	Input: 
	m8_filename = filename of blast results
	fasta_filename = filename of corresponding fasta file
	
	Output: 
	table = {gene: [List of genomes]}'''

	fasta_dict = utilities.read_fasta(fasta_filename)

	for seq in fasta_dict:
		fasta_dict[seq] = float(len(fasta_dict[seq]))
	table = {}
	foo = open(m8_filename)
	for line in foo:
		split_i = line.split('\t')
		try:
			threshold = float(split_i[2])*float(split_i[3])/fasta_dict[split_i[0]]
		except:
			raise Exception('Gene '+split_i[0]+' not found in fasta: '+fasta_filename)
		if threshold > 90.0:
			split_sp = split_i[1].split('|')
			sp = [i for i in split_sp if 'g__' in i and '.s__' in i]
			if split_i[0] not in table:
				table[split_i[0]] = sp
			elif sp not in table[split_i[0]]:
				table[split_i[0]] += sp
	return table

def read_parsed(m8_filename):
	'''Read parsed table for {gene: genomes}
	Input: 
	m8_filename = filename of blast results

	Output: 
	table = {gene: [List of genomes]}'''

	table = {}
	foo = open(m8_filename)

	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		try:
			table[split_i[0]] += [split_i[1]]
		except:
			table[split_i[0]] = [split_i[1]]
	print"Total No. of genes:", len(table) 
	return table
def read_data(mg_file, ppanini_output_file):
	metagenomic_table  = read_parsed(mg_file)
	uniq_genomes = []
	for gene in metagenomic_table:
		for genome in metagenomic_table[gene]:
			if genome not in uniq_genomes:
				uniq_genomes +=[genome]
	no_uniq_genomes = len(uniq_genomes)
	print 'No. of unique genomes: '+str(no_uniq_genomes)
	ppanini_output = read_ppanini_imp_genes_table(ppanini_output_file)
	return metagenomic_table, ppanini_output, no_uniq_genomes 
def pullgenes_fromcontigs(contig_file, gff3_file, fna_file, faa_file):
	'''Pulls genes from contigs using the coordinates from GFF3 file provided

	Input: contig_file: FASTA file containing contigs in FNA format.
		   gff3_file: File containing coordinates of genes in GFF3 format
		   fna_file: filepath for genes written in nucleotides sequences
		   faa_file: filepath for genes written in amino-acid sequences'''
	logger.debug('pullgenes_fromcontigs '+contig_file+' '+gff3_file)

	gene_contig_mapper = {}
	contig_gene_mapper = {}
	gene_start_stop = {}

	[gene_contig_mapper, gene_start_stop, contig_gene_mapper] = read_gff3(gff3_file)

	genes_fasta = {}
	contigs_fasta_dict = read_fasta(contig_file)
	for contig in contigs_fasta_dict:
		if contig in contig_gene_mapper:
			for gene in contig_gene_mapper[contig]:
				[start_x, stop_x, strand] = gene_start_stop[gene]
				try:
					if strand == '+':
						genes_fasta[gene] = contigs_fasta_dict[contig][start_x-1:stop_x+1]
					else:
						contig_len = len(contigs_fasta_dict[contig])
						start_minus = -1*stop_x
						stop_minus = -1*(start_x-1)
						genes_fasta[gene] = Bio.Seq.reverse_complement(contigs_fasta_dict[contig])[start_minus:stop_minus]
				except:
					raise Exception('Circular DNA Detected')
					new_stop_x = stop_x - len(contigs_fasta_dict[contig])+1
					genes_fasta[gene] = contigs_fasta_dict[contig] + contigs_fasta_dict[contig][:new_stop_x]
	
	write_fasta(genes_fasta, fna_file, False) #FNA
	write_fasta(genes_fasta, faa_file, True) #FAA


def read_gff3(filename):
	'''Reads GFF3 files and returns the relevant information
	Input: filename = path_to_gff3_file
	Output: [gene_contig_mapper, gene_start_stop, contig_gene_mapper]
	gene_contig_mapper = {gene: contig, ...}
	gene_start_stop = {gene: [start, stop, strand], ...}
	contig_gene_mapper = {contig: [List of genes], ...}'''
	
	logger.debug('read_gff3 '+filename)

	gene_contig_mapper = {}
	contig_gene_mapper = {}
	gene_start_stop = {}

	with open(filename,'r') as foo:
		for line in foo:
			if re.match('(\w+)\t(\w+)\tgene\t(\w+)', line):
				split_line = [re.sub('[\r\t\n]', '', i).strip() for i in line.split('\t')]
				gid = split_line[-1].split('=')[-1]
				gene_contig_mapper[gid] = split_line[0]
				gene_start_stop[gid] = [int(split_line[3]), int(split_line[4]), split_line[6]]
				if split_line[0] in contig_gene_mapper:
					contig_gene_mapper[split_line[0]] += [gid]
				else:
					contig_gene_mapper[split_line[0]] = [gid]
	return [gene_contig_mapper, gene_start_stop, contig_gene_mapper]


def write_fasta(seqs_dict, filename, to_translate):
	'''Writes dictionary of fasta sequences into text file
	Input: seqs_dict = {geneID: sequence}
		   filename = path_to_output_genes_fastafile'''
	logger.debug('write_fasta to '+filename)

	with open(filename,'w') as foo:

		test = ''
		for i in seqs_dict:
			test = seqs_dict[i]
			break
		format = is_protein(test)
		to_translate =  to_translate and not format
		
		if to_translate: # if not FAA already and to be translated
			for seq in seqs_dict:
				t_seq = Bio.Seq.translate(seqs_dict[seq], to_stop=True)
				len_seq = str(len(t_seq)*3)
				foo.writelines(['>' + seq + '|' + len_seq + '\n'])
				foo.writelines([t_seq+'\n'])
		else:
			ind = 1
			if format:
				ind = 3 #amino acids * 3 nucleotides
			for seq in seqs_dict:
				len_seq = str(len(seqs_dict[seq]*ind))
				foo.writelines(['>' + seq + '|' + len_seq + '\n'])
				foo.writelines([seqs_dict[seq] + '\n'])

def is_protein(sequence):
	'''Returns True if the sequence is protein.
	Input: (str) format sequence
	Output: (bool) True if amino acids; False if nucleotides sequence'''

	logger.debug('is_protein')

	format = False  #'FNA'
	try:
		Bio.Seq.translate(sequence)
	except:
		format = True #'FAA'
	return format

def create_folders(list_folders):
	'''Creates the list of folders if they dont exist already
	Input: list_folders = [List of folders to be created]

	Output: Folders created'''

	logger.debug('create_folders '+'\t'.join(list_folders))
	
	for fname in list_folders:
		try:
			os.stat(fname)
		except:
		    os.mkdir(fname)
		

def read_dict(gene_annotations_file):
	'''Reads tabulated file into a dictionary

	Input: gene_annotations_file = path_to_output_gene_annotations_table

	Output: dictX = {geneID: annotation}'''

	logger.debug('read_dict '+gene_annotations_file)

	dictX = {}

	with open(gene_annotations_file) as foo:
		for line in foo:
			if not line.startswith('#'):
				split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split('\t')]
				dictX[split_line[0]] = split_line[1]
	return dictX

def read_dict_num(gene_annotations_file):
	'''Reads tabulated file into a dictionary

	Input: gene_annotations_file = path_to_output_gene_annotations_table

	Output: dictX = {geneID: annotation}'''

	logger.debug('read_dict)num '+gene_annotations_file)

	dictX = {}

	with open(gene_annotations_file) as foo:
		for line in foo:
			if not line.startswith('#'):
				split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split('\t')]
				if len(split_line) <2:
					split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split(' ')]
				try:
					dictX[split_line[0]] = float(split_line[1])
				except:
					pdb.set_trace()
	return dictX

def write_dict(dictX, gene_annotations_file):
	'''Writes dictionary of genes and their annotations into text file
	Input: dictX = {geneID: annotation}
		   gene_annotations_file = path_to_output_gene_annotations_table'''

	logger.debug('write_dict to '+gene_annotations_file)

	with open(gene_annotations_file, 'w') as foo:
		foo.writelines('#GENEID\tANNOTATION\n')
		for i in dictX:
			foo.writelines(['\t'.join([i, dictX[i] ])+'\n'])

def is_present(metadata, meta_type):
	'''Returns True if meta_type is present in metadata extracted from mappert_file

	Input: metadata = [metadata strings]; Rows with # as first character in table
		   meta_type = Type of metadata that you are querying e.g. FASTAS, NICHE etc.

	Output: [line, ind]
			line = The corresponding line from metadata, [] if not present
			ind = The index of the line in the metadata sequence, [] if not present'''

	logger.debug('is_present '+meta_type)

	line = []
	ind = []
	for i, val in enumerate(metadata):
		if val.upper().startswith(meta_type):
			line = val
			ind = i
			break
	return [line, ind]

if __name__ == '__main__':
	pass
