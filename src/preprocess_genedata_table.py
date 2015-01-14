import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing

def parse_mapper(mapper_file):
# #SAMPLE<\t>#READS<\t>#ASSEMBLIES<\t>#NICHE?<\t>#GENES<\t>CONTIG_GENE
	mapper = {}
	mapper_foo = open(mapper_file)
	mapper_foo = mapper_foo.readlines()

	header = [re.sub('[\t\n\r]', '', i) for i in mapper_foo[0].split('\t')]

	niche_flag = [True for i in header if i=='#NICHE']
	
	for i, val_i in enumerate(mapper_foo):
		split_val_i = [re.sub('[\t\n\r]', '', i) for i in val_i.split('\t')]
		mapper[i] = {}
		for j, val_j in enumerate(split_val_i):
			mapper[i][header[j]] = val_j

	return [mapper, niche_flag]


def generate_abundance(assembly_x_withpath, reads_x):
	#assembly_x_withpath = mapper[sample]['#ASSEMBLIES'] #Full path inlcuded to Assembly X
	assembly_x = assembly_x_withpath.rpartition('/')[-1] #Name of Assembly X
	#reads_x = mapper[sample]['#READS'] #Path included to Reads of X
	assembly_x_index = 'tmp/' + assembly_x + '_index' # 
	assembly_x_bam = 'tmp/' + assembly_x + '.bam'
	assembly_x_bam_presort = 'tmp/' + assembly_x + '.sorted'
	assembly_x_bam_sorted = 'tmp/' + assembly_x + '.sorted.bam'
	assembly_x_stats = 'tmp/' + assembly_x + '.txt'
	
	os.system('bowtie2-build --quiet ' + assembly_x_withpath + ' ' + assembly_x_index)
	os.system('tar -xOvf ' + reads_x + ' | \
			   bowtie2 -x ' + assembly_x_index + ' -U - --no-unal --very_sensitive | \
			   samtools view -bS - > ' + assembly_x_bam)
	os.system('samtools sort ' + assembly_x_bam + ' ' + assembly_x_bam_presort)
	os.system('samtools index ' + assembly_x_bam_sorted)
	os.system('samtools idxstats ' + assembly_x_bam_sorted + ' > ' + assembly_x_stats)
	
	return True

def get_abundance(mapper, all_paths):
	#tar -xOvf reads__file | bowtie2 -x filename_index -U - --no-unal --very_sensitive | \
	#samtools view -bS - > newfile.bam
	#samtools sort filename newfilename.sorted
	#samtools index newfilename.sorted.bam
	#samtools idxstats newfilename > statsfilename

	pool = multiprocessing.Pool(processes=4)
	results = [pool.apply_async(generate_abundance, args=(mapper[sample]['#ASSEMBLIES'], mapper[sample]['#READS']) for sample in mapper]
	
	out = sum([p.get() for p in results])
	
	if out/len(results) < 1:
		raise Exception('The Abundance calculations encountered an issue; please check the code')

	for sample in mapper:
		assembly_x_stats = 'tmp/' + mapper[sample]['#ASSEMBLIES'].rpartition('/')[-1] + '.txt'
		mapper[sample]['abund_file'] = assembly_x_stats

	return read_abundance_tables(mapper, all_paths)

def generate_annotation(gene_x_file, all_paths):
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
	
	u90_out_fname =  'tmp/u90/' + gene_x_fname
	u50_out_fname =  'tmp/u50/' + gene_x_fname
	genes_header_indices = 'tmp/' + genes_x_fname + '_header_indices.txt'
	parsed_u90_out = 'tmp/' + genes_x_fname + '.m8'  
	u50_gene_input = 'tmp/u50input_' + genes_x_fname
	allannot_toparse = 'tmp/allannot_toparse_' + genes_x_fname+'.m8'
	gene_annotations_file = 'tmp/annot/allannot_toparse_' + genes_x_fname + '.m8'

	os.system('rapsearch -q ' + gene_x_file + ' -d ' + all_paths['uniref90'] + ' -o ' + u90_out_fname + ' -u 2 -b 0 -v 1')
	os.system('python utils/parse_annot.py ' + u90_out_fname + '.m8 tmp ' + genes_x_file) #produces parsed_u90_out
	os.system('cat ' + genes_x_file + ' | grep -n \'>\' > '+ genes_header_indices)
	os.system('python utils/gen_u50runner.py ' + genes_x_file + ' ' + genes_header_indices + ' ' + parsed_u90_out + ' ' + u50_gene_input)
	if os.path.isfile(u50_gene_input):
		os.system('rapsearch -q ' + u50_gene_input + ' -d ' + all_paths['uniref50'] + ' -o ' + u50_out_fname + ' -u 2 -b 0 -v 1')
		os.system('cat ' + parsed_u90_out +' '+ u50_out_fname +'.m8 > ' + allannot_toparse)
		os.system('python utils/parse_annot.py ' + allannot_toparse + ' tmp/annot ' + genes_x_file) #produces allannot_toparse under annot
	else:
		os.system('cp ' + parsed_u90_out + ' '+ gene_annotations_file)

	return True

def get_annotation(mapper, all_paths):
	#tar -xOvf reads__file | bowtie2 -x filename_index -U - --no-unal --very_sensitive | \
	#samtools view -bS - > newfile.bam
	#samtools sort filename newfilename.sorted
	#samtools index newfilename.sorted.bam
	#samtools idxstats newfilename > statsfilename

	pool = multiprocessing.Pool(processes=4)
	results = [pool.apply_async(generate_annotation, args=(mapper[sample]['#GENES'], all_paths) for sample in mapper]
	
	out = sum([p.get() for p in results])
	
	if out/len(results) < 1:
		raise Exception('The Annotation calculations encountered an issue; please check the code')

	for sample in mapper:
		genes_x_annotations_file = 'tmp/annot/allannot_toparse_' + gene_x_file.rpartition('/')[-1] + '.m8'
		mapper[sample]['abund_file'] = assembly_x_stats

	return read_abundance_tables(mapper, all_paths)

		
def read_abundance_tables(mapper, all_paths):
	abundance_dict = {}
	for sample in mapper:
		abundance_dict[mapper[sample]['#SAMPLE']] = {}
		foo = open(mapper[sample]['abund_file'])
		for line in foo.readlines():
			if not (line.startswith('#') and line.startswith('*')):
				#FILE_FORMAT: Seq_Name<\t>Seq_Length<\t>Mapped_Reads<\t>Unmapped_Reads
				split_line = [re.sub('[\t\r\n]', '', i) for i in line.split('\t')]
				#abundance_dict[sample][contig_name] = [number_mapped_reads, sequence_length]
				abundance_dict[mapper[sample]['#SAMPLE']][split_line[0].strip()] = float(split_line[2])/(float(split_line[1])-100.0)
	return abundance_dict


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	#parser.add_argument('-a', '--assemblies_folder', help='Folder containing ORF-called contig assemblies of samples')
	#parser.add_argument('-r', '--reads_folder', help='Folder containing reads from samples')
	#parser.add_argument('-n', '--niche_present', choices=[True, False], help='Is Niche-classification present? [True or False]]')
	parser.add_argument('-m', '--mapper_file', help='Mapper file associating read files with their corresponding assembly files and, if available, GFF3 files and NICHE informations')
	#parser.add_argument('-g', '--gff3_folder', help='Folder containing GFF3 files')
	parser.add_argument('-umap', '--uniref90_50', help='UniRef90 XML file')
	parser.add_argument('-u90', '--uniref90_fasta', help='UniRef90 fasta file')
	parser.add_argument('-u50', '--uniref50_fasta', help='UniRef50 fasta file')
	#parser.add_argument('-rs2', '--rapsearch2_path', help='Path to run RAPSEARCH2') #
	#parser.add_argument('-st', '--samtools_path', help='Path to run SAMTOOLS')
	#parser.add_argument('-bwt2', '--bowtie2_path', help='Path to run bowtie2_path')

	args = parser.parse_args()

	try:
		os.mkdir('tmp')
	except:
		pass

	#'assemblies': args.assemblies_folder, \
				 #'reads': args.reads_folder, \
	all_paths = {'uniref_map': args.uniref90_50, \
				 'uniref90': args.uniref90_fasta, \
				 'uniref50': args.uniref50_fasta, \
				 'rapsearch2': args.rapsearch2_path, \
				 'samtools': args.samtools_path, \
				 'bowtie2': args.bowtie2_path}

	#if args.gff3_folder is not None:
	#	all_paths['gff3s'] = args.gff3_folder

	[mapper, niche_flag] = parse_mapper(args.mapper_file)
	abundance_dict = get_abundance(mapper, all_paths) #sample:contig: mapped reads/seq_length
	get_annotations(mapper, all_paths)

