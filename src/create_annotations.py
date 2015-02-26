import os
import sys
import numpy
import re
import Bio
from Bio import Seq
import argparse


def read_fasta(fasta_filename):
	'''Reads a fasta_file and returns a fasta dict
	Input: fasta_filename = path_to_fasta_file
	Output: fasta_seq = {sequence_header: sequence, ...}'''
	
	fasta_file = open(fasta_filename)
	fasta_seq = {}

	name = ''
	
	for line in fasta_file.readlines():
		if not line.startswith('#'):
			if line.startswith('>'):
				name = line.split(' ')[0][1:].strip()
			else:
				if name not in fasta_seq:
					fasta_seq[name] =  re.sub('[\r\t\n]','', line)
				else:
					fasta_seq[name] +=  re.sub('[\r\t\n]','', line)
	return fasta_seq

def parse_annotation_table(annotations_file, fasta_sequences, thld_ref):
	'''Parses annotations result from RAPSEARCH to ensure they pass threshold
	Input: annotations_file = path_to_rapsearch_results_file
		   fasta_sequences = {sequence_header: sequence}
		   thld_ref = int, threshold for sequence acceptance
	
	Output: [sample_annotations, rapsearch50_seqs]
	where, sample_annotations = {sequence: UniRef annotation}
		   rapsearch50_seqs = {sequence_header:  sequence} for all sequences that need to be redone'''

	foo = open(annotations_file)
	sample_annotations = {}
	for line in foo.readlines():
		if not line.startswith('#'):
			#FILE_FORMAT: Seq_Name<\t>Seq_Length<\t>Mapped_Reads<\t>Unmapped_Reads
			split_line = [re.sub('[\t\r\n]', '', i) for i in line.split('\t')]

			query_length = float(len(fasta_sequences[split_line[0]]))
			aln_len = float(split_line[3])
			identity = float(split_line[2])
			thld = aln_len*identity/query_length
			
			if thld >= thld_ref:
				sample_annotations[split_line[0].strip()] = split_line[1].strip()
	
	rapsearch50_seqs = {}
	for seq in fasta_sequences:
		if not seq in sample_annotations:
			rapsearch50_seqs[seq] = fasta_sequences[seq]

	return [sample_annotations, rapsearch50_seqs]

def write_dict(dictX, gene_annotations_file):
	'''Writes dictionary of genes and their annotations into text file
	Input: dictX = {geneID: annotation}
		   gene_annotations_file = path_to_output_gene_annotations_table
	'''

	with open(gene_annotations_file, 'w') as foo:
		foo.writelines('#GENEID\tANNOTATION')
		for i in dictX:
			foo.writelines([str.join('\t', [i, dictX[i]])+'\n'])

def write_fasta(seqs_dict, filename):
	'''Writes dictionary of fasta sequences into text file
	Input: seqs_dict = {geneID: sequence}
		   filename = path_to_output_genes_fastafile
	'''
	with open(filename,'w') as foo:
		test = seqs_dict.values()[0]
		format = True  #'FNA'
		try:
			Bio.Seq.translate(test) 
		except:
			format = False #'FAA'

		if format:
			for seq in seqs_dict:
				foo.writelines(['>'+seq+'\n'])
				foo.writelines([Bio.Seq.translate(seqs_dict[seq], to_stop=True)+'\n'])
		else:
			for seq in seqs_dict:
				foo.writelines(['>'+seq+'\n'])
				foo.writelines([seqs_dict[seq]+'\n'])

def generate_annotation(gene_x_file, all_paths, nprocesses):
	'''Runs rapsearch2 for genes_fasta file against UniRef90 and UniRef50
	Input: gene_x_file: path_to_genes_fasta_file
		   all_paths = {'uniref90': path_to_uniref90_index, 
		   				'uniref50': path_to_uniref50_index, 
		   				'umap90_50': path_to_uniref90_uniref50_mapping}
		   nprocesses = int, Number of processes
		   sample = sample_name

	Output: all_annotations = {gene: UniRef annotation}'''

	gene_x_fname = gene_x_file.rpartition('/')[-1] 
	try:
		os.mkdir('tmp/u90')
	except:
		pass
	
	try:
		os.mkdir('tmp/u50')
	except:
		pass

	try:
		os.mkdir('tmp/annot')
	except:
		pass
	
	fasta_sequences = read_fasta(gene_x_file)
	
	all_annotations = {}

	u90_out_fname =  'tmp/u90/' + gene_x_fname
	u50_out_fname =  'tmp/u50/' + gene_x_fname 
	u50_gene_input = 'tmp/u50input_' + gene_x_fname
	gene_annotations_file = 'tmp/annot/allannot_parsed_' + gene_x_fname + '.m8'

	os.system(all_paths['rapsearch']+'/rapsearch -q ' + gene_x_file + ' \
						 -d ' + all_paths['uniref90'] + ' \
						 -o ' + u90_out_fname + ' \
						 -u 2 \
						 -b 0 \
						 -v 1 \
						 -z' + str(nprocesses))
	
	[sample_annotations90, rapsearch50_seqs] = parse_annotation_table(u90_out_fname+'.m8', fasta_sequences, 90.0)
	
	all_annotations = sample_annotations90

	if not rapsearch50_seqs == {}:
		write_fasta(rapsearch50_seqs, u50_gene_input)
		os.system(all_paths['rapsearch']+'/rapsearch -q ' + u50_gene_input + ' \
							 -d ' + all_paths['uniref50'] + '\
							 -o ' + u50_out_fname + ' \
							 -u 2 \
							 -b 0 \
							 -v 1 \
							 -z ' + str(nprocesses))
		
		[sample_annotations50, rapsearchukn_seqs] = parse_annotation_table(u50_out_fname+'.m8', rapsearch50_seqs, 50.0)
		
		for gid in sample_annotations50:
			if gid not in all_annotations:
				all_annotations[gid] = sample_annotations50[gid]
			else:
				raise Exception('GeneID for UniRef50 exists in Uniref90; something wrong with annotations!!!')
	
	write_dict(all_annotations, gene_annotations_file)
	return all_annotations


if __name__ == '__main__':
	#python create_annotations.py fasta uniref90 uniref50 rapsearch_path nprocesses >> gives all_annotations files
	gene_x_file = sys.argv[1]
	nprocesses = int(sys.argv[5])
	all_paths = {'uniref90': sys.argv[2], 'uniref50': sys.argv[3], 'rapsearch': sys.argv[4]}
	
	generate_annotation(gene_x_file, all_paths, nprocesses)
