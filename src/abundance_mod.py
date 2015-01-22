import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing

def generate_abundance_viabwt2(assembly_x_withpath, reads_x):
	'''To go from ASSEMBLIES, READS to IDXSTATS'''
	#assembly_x_withpath = mapper[sample]['#ASSEMBLIES'] #Full path inlcuded to Assembly X
	assembly_x = assembly_x_withpath.rpartition('/')[-1] #Name of Assembly X
	#reads_x = mapper[sample]['#READS'] #Path included to Reads of X
	assembly_x_index = 'tmp/' + assembly_x + '_index' # 
	assembly_x_bam = 'tmp/' + assembly_x + '.bam'
	
	os.system('bowtie2-build --quiet ' + assembly_x_withpath + ' ' + assembly_x_index)
	os.system('tar -xOvf ' + reads_x + ' | \
			   bowtie2 -x ' + assembly_x_index + ' -U - --no-unal --very_sensitive | \
			   samtools view -bS - > ' + assembly_x_bam)
	generate_abundance_bam(assembly_x_bam)
	
	return True

def generate_abundance_viasam(assembly_x_sam_withpath):
	'''To go from SAM to IDXSTATS'''
	
	assembly_x_sam = assembly_x_sam_withpath.rpartition('/')[-1] #Name of Assembly X
	assembly_x_bam = 'tmp/' + assembly_x_sam + '.bam'
	
	os.system('samtools view -bS ' + assembly_x_sam_withpath + ' > ' + assembly_x_bam)
	generate_abundance_bam(assembly_x_bam)
	
	return True

def generate_abundance_viabam(assembly_x_bam_withpath):
	'''To go from BAM to IDXSTATS'''
	assembly_x_bam = assembly_x_withpath.rpartition('/')[-1]
	assembly_x_bam_presort = 'tmp/' + assembly_x_bam + '.sorted'
	assembly_x_bam_sorted = 'tmp/' + assembly_x_bam + '.sorted.bam'
	assembly_x_stats = 'tmp/' + assembly_x_bam + '.txt'
	
	os.system('samtools sort ' + assembly_x_bam_withpath + ' ' + assembly_x_bam_presort)
	os.system('samtools index ' + assembly_x_bam_sorted)
	os.system('samtools idxstats ' + assembly_x_bam_sorted + ' > ' + assembly_x_stats)
	return True


def get_abundance(mapper, nprocesses, w_i):
	
	workflows = {1: [generate_abundance_viabam,'BAMS'], 2: [generate_abundance_viasam, 'SAMS'], 3: [generate_abundance_viabwt2, 'ASSEMBLIES']}
	
	pool = multiprocessing.Pool(processes=nprocesses)
	
	if w_i == 1:
		results = [pool.apply_async(workflows[w_i][0], args=(mapper[sample]['BAMS'])) for sample in mapper]
	elif w_i == 2:
		results = [pool.apply_async(workflows[w_i][0], args=(mapper[sample]['SAMS'])) for sample in mapper]
	elif w_i == 3:
		results = [pool.apply_async(workflows[w_i][0], args=(mapper[sample]['ASSEMBLIES'], mapper[sample]['READS'])) for sample in mapper]
	else:
		raise Exception('Invalid workflow inserted! Should be either 1, 2 or 3')
	
	out = sum([p.get() for p in results])
	
	if out/len(results) < 1:
		raise Exception('The Abundance calculations encountered an issue; please check the code')

	for sample in mapper:
		assembly_x_stats = 'tmp/' + mapper[sample][workflows[w_i][1]].rpartition('/')[-1] + '.txt'
		mapper[sample]['abundance_file'] = assembly_x_stats


	return read_abundance_tables(mapper, all_paths)

		
def read_abundance_tables(mapper):
	abundance_dict = {}
	for sample in mapper:
		abundance_dict[sample] = {}
		foo = open(mapper[sample]['abundance_file'])
		for line in foo.readlines():
			if not (line.startswith('#') and line.startswith('*')):
				#FILE_FORMAT: Seq_Name<\t>Seq_Length<\t>Mapped_Reads<\t>Unmapped_Reads
				split_line = [re.sub('[\t\r\n]', '', i) for i in line.split('\t')]
				#abundance_dict[sample][contig_name] = [number_mapped_reads, sequence_length]
				abundance_dict[mapper[sample]][split_line[0].strip()] = float(split_line[2])/(float(split_line[1])-100.0) #avg. read length?
	
	return abundance_dict