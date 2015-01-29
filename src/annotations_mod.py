import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing
import Bio
from Bio import Seq

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


def generate_annotation(gene_x_file, all_paths, nprocesses, sample):
	'''Runs rapsearch2 for genes_fasta file against UniRef90 and UniRef50
	Input: gene_x_file: path_to_genes_fasta_file
		   all_paths = {'uniref90': path_to_uniref90_index, 
		   				'uniref50': path_to_uniref50_index, 
		   				'umap90_50': path_to_uniref90_uniref50_mapping}
		   nprocesses = int, Number of processes
		   sample = sample_name

	Output: all_annotations = {gene: UniRef annotation}'''

	#assembly_x_withpath = mapper[sample]['#ASSEMBLIES'] #Full path inlcuded to Assembly X
	
	gene_x_fname = gene_x_file.rpartition('/')[-1] #Name of Assembly X
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

def get_annotations_fromgenes(mapper, all_paths, nprocesses):
	'''Return dictionary of annotations for all samples' genes fasta files
	Input: mapper = {sample:{'READS': path_to_reads_file, 
		   					 'CONTIG_ASSEMBLIES': path_to_assemblies, 
		   					 'FASTAS': path_to_fasta_file,
		   					 'SAMS': path_to_sam_file,
		   					 'BAMS': path_to_bam_file,
		   					 'GFF3S': path_to_gff3_file,
		   					 'NICHE': niche_metadata}, ...}
		   all_paths = {'uniref90': path_to_uniref90_index, 
		   				'uniref50': path_to_uniref50_index, 
		   				'umap90_50': path_to_uniref90_uniref50_mapping}
		   nprocesses = int, Number of processes

	Output: annotations_dict = {sample :{gene: UniRef annotation}}'''
	pool = multiprocessing.Pool(processes=nprocesses)
	results = [pool.apply_async(generate_annotation, args=(mapper[sample]['FASTAS'], all_paths, nprocesses, sample)) for sample in mapper]

	out = sum([1 for p in results])
	
	if out/len(results) < 1:
		raise Exception('The Annotation calculations encountered an issue; please check the code')

	annotations_dict = {}
	for p in results:
		tmp_dict = p.get()
		annotations_dict[tmp_dict['key']] = tmp_dict
	
	return annotations_dict

def get_annotations_fromcontigs(mapper, all_paths, nprocesses):
	'''Return dictionary of annotations for genes from each sample's contig_assembly

	Input: mapper = {sample:{'READS': path_to_reads_file, 
		   					 'CONTIG_ASSEMBLIES': path_to_assemblies, 
		   					 'FASTAS': path_to_fasta_file,
		   					 'SAMS': path_to_sam_file,
		   					 'BAMS': path_to_bam_file,
		   					 'GFF3S': path_to_gff3_file,
		   					 'NICHE': niche_metadata}, ...}
		   all_paths = {'uniref90': path_to_uniref90_index, 
		   				'uniref50': path_to_uniref50_index, 
		   				'umap90_50': path_to_uniref90_uniref50_mapping}
		   nprocesses = int, Number of processes

	Output: [annotations_dict, gene_contig_mapper]
	annotations_dict = {sample :{gene: UniRef annotation}}
	gene_contig_mapper = {sample: {gene: contig}}'''	

	gene_contig_mapper = {}
	contig_gene_mapper = {}
	gene_start_stop = {}
	
	for sample in mapper:
		[gene_contig_mapper_persample, gene_start_stop_persample, contig_gene_mapper_persample] = read_gff3(mapper[sample]['GFF3S'])
		gene_start_stop[sample] = gene_start_stop_persample
		contig_gene_mapper[sample] = contig_gene_mapper_persample
		gene_contig_mapper[sample] = gene_contig_mapper_persample

	mapper_withgenes = pull_genes_fromcontigs(mapper, gene_start_stop, contig_gene_mapper)

	annotations_dict = get_annotations_fromgenes(mapper_withgenes, all_paths, nprocesses)

	return [annotations_dict, gene_contig_mapper]

def pull_genes_fromcontigs(mapper, gene_start_stop, contig_gene_mapper):
	''' Pull genes from contig_assemblies and returns mapper with created genes_fasta_files
	Input: mapper = {sample:{'READS': path_to_reads_file, 
		   					 'CONTIG_ASSEMBLIES': path_to_assemblies, 
		   					 'GFF3S': path_to_gff3_file,
		   					 'NICHE': niche_metadata}, ...}
		   gene_start_stop = {sample: {gene: [start, stop, strand]}, ...}
		   contig_gene_mapper = {sample: {contig: [List of genes]}}
	
	Output: mapper = {sample:{'READS': path_to_reads_file, 
		   					 'CONTIG_ASSEMBLIES': path_to_assemblies, 
		   					 'FASTAS': path_to_fasta_file,
		   					 'GFF3S': path_to_gff3_file,
		   					 'NICHE': niche_metadata}, ...}'''
	try:
		os.mkdir('tmp/fasta_files')
	except:
		pass

	for sample in mapper:
		genes_fasta = {}
		contigs_fasta_dict = read_fasta(mapper[sample]['CONTIG_ASSEMBLIES'])
		for contig in contigs_fasta_dict:
			if contig in contig_gene_mapper[sample]:
				for gene in contig_gene_mapper[sample][contig]:
					[start_x, stop_x, strand] = gene_start_stop[sample][gene]
					try:
						if strand == '+':
							genes_fasta[gene] = contigs_fasta_dict[contig][start_x-1:stop_x+1]
						else:
							contig_len = len(contigs_fasta_dict[contig])
							start_minus = contig_len - stop_x
							stop_minus = contig_len - start_x+1
							genes_fasta[gene] = Bio.Seq.reverse_complement(contigs_fasta_dict[contig])[start_minus:stop_minus]
					except:
						raise Exception('Circular DNA Detected')
						new_stop_x = stop_x - len(contigs_fasta_dict[contig])+1
						genes_fasta[gene] = contigs_fasta_dict[contig] + contigs_fasta_dict[contig][:new_stop_x]
	
		sample_fasta_filename = 'tmp/fasta_files/'+sample+'.fasta'
		write_fasta(genes_fasta, sample_fasta_filename)
		mapper[sample]['FASTAS'] = sample_fasta_filename
	return mapper

def read_gff3(filename):
	'''Reads GFF3 files and returns the relevant information
	Input: filename = path_to_gff3_file
	Output: [gene_contig_mapper, gene_start_stop, contig_gene_mapper]
	gene_contig_mapper = {gene: contig, ...}
	gene_start_stop = {gene: [start, stop, strand], ...}
	contig_gene_mapper = {contig: [List of genes], ...}'''

	gene_contig_mapper = {}
	contig_gene_mapper = {}
	gene_start_stop = {}

	with open(filename,'r') as foo:
		foo_lines = foo.readlines()
		for line in foo_lines:
			if re.match('(\w+)\t(\w+)\tgene\t(\w+)', line):
				split_line = [re.sub('[\r\t\n]','', i).strip() for i in line.split('\t')]
				gid = split_line[-1].split('=')[-1]
				gene_contig_mapper[gid] = split_line[0]
				gene_start_stop[gid] = [int(split_line[3]), int(split_line[4]), split_line[6]]
				
				if split_line[0] in contig_gene_mapper:
					contig_gene_mapper[split_line[0]] += [gid]
				else:
					contig_gene_mapper[split_line[0]] = [gid]

	return [gene_contig_mapper, gene_start_stop, contig_gene_mapper]
