import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing


def read_fasta(foo):
	'''Returns fasta dict with Seqname:sequence'''
	foo = open(foo)
	fasta_seq = {}

	name = ''
	
	for line in foo.readlines():
		if not line.startswith('#'):
			if line.startswith('>'):
				name = line.split(' ')[0][1:].strip()
			else:
				if name in fasta_seq:
					fasta_seq[name] =  re.sub('[\r\t\n]','', line)
				else:
					fasta_seq[name] +=  re.sub('[\r\t\n]','', line)
	return fasta_seq

def parse_annotation_table(annotations_file, fasta_sequences, thld_ref):
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

	with open(gene_annotations_file, 'w') as foo:
		foo.writelines('#GENEID\tANNOTATION')
		for i in dictX:
			foo.writelines([str.join('\t', [i, dictX[i]])+'\n'])


def generate_annotation(gene_x_file, all_paths, nprocesses, sample):
	#assembly_x_withpath = mapper[sample]['#ASSEMBLIES'] #Full path inlcuded to Assembly X
	genes_x_fname = gene_x_file.rpartition('/')[-1] #Name of Assembly X
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
	u50_gene_input = 'tmp/u50input_' + genes_x_fname
	gene_annotations_file = 'tmp/annot/allannot_parsed_' + genes_x_fname + '.m8'

	os.system('rapsearch -q ' + gene_x_file + ' \
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
		os.system('rapsearch -q ' + u50_gene_input + ' \
							 -d ' + all_paths['uniref50'] + '\
							 -o ' + u50_out_fname + ' \
							 -u 2 \
							 -b 0 \
							 -v 1')
		
		[sample_annotations50, rapsearchukn_seqs] = parse_annotation_table(u50_out_fname+'.m8', rapsearch50_seqs, 50.0)
		
		for gid in sample_annotations50:
			if gid not in all_annotations:
				all_annotations[gid] = sample_annotations50[gid]
			else:
				raise Exception('GeneID for UniRef50 exists in Uniref90; something wrong with annotations!!!')
	all_annotations['key'] = sample
	
	write_dict(all_annotations, gene_annotations_file)

	return all_annotations

def write_fasta(seqs_dict, filename):
	with open(filename,'w') as foo:
		for seq in seqs_dict:
			foo.writelines(['>'+seq+'\n'])
			foo.writelines([seqs_dict[seq]+'\n'])

def get_annotation(mapper, all_paths, nprocesses):

	pool = multiprocessing.Pool(processes=nprocesses)
	results = [pool.apply_async(generate_annotation, args=(mapper[sample]['FASTA'], all_paths, nprocesses, sample)) for sample in mapper]
	
	out = sum([1 for p in results])
	
	if out/len(results) < 1:
		raise Exception('The Annotation calculations encountered an issue; please check the code')

	annotations_dict = {}
	for p in results:
		tmp_dict = p.get()
		annotations_dict[tmp_dict['key']] = tmp_dict

	#mapper['annotations_dict'] = annotations_dict
	
	return annotations_dict

