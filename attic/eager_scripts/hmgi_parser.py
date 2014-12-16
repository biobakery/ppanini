# Shafquat, Afrah
# shafquat@hsph.havard.edu
# March, 2014

import os
import sys
import re
import pdb


def cluster_reads(foo, data_folder):
	'''Returns a dict with HMGC centroids mapped to list of HMGI IDs Input 
    * foo = filepath to the HMP uclust output where the 10th column denotes the GC ID and 9th column denotes the GI ID
    * data_folder= path to output folder, where hmgc_gi.txt file will be created

    Output
    * hmgc_gi = {'SRS0XXX.XXX': ['SRS0XXXX:XXX', ...]}'''
	print 'Processing HMGC: HMGI mapping'
	foo = open(foo)
	foo = foo.readlines()
	foo = foo[1:]  # First line is assumend to be header
	hmgc_gi = {}
	for i in foo:  # Reads the file, 10th column is GC ID, 9th is GI ID
		split_i = i.split('\t')
		if not re.sub('[\t\n\r]', '', split_i[9]).strip().split('-')[0] in hmgc_gi:
			hmgc_gi[re.sub('[\t\n\r]', '', split_i[9]).strip().split('-')[0]] = [re.sub('[\t\n\r]', '', split_i[8]).strip().split('-')[0]]
		else:
			hmgc_gi[re.sub('[\t\r\n]', '', split_i[9]).strip().split('-')[0]] += [re.sub('[\t\n\r]', '', split_i[8]).strip().split('-')[0]] # writes the hmgc_gi.txt to the output folder
	for gc in hmgc_gi:
		if not gc in hmgc_gi[gc]:
			hmgc_gi[gc] += [gc]
	write_dict(hmgc_gi, 'file','hmgc_gi', data_folder)
	return hmgc_gi


def gi_reads(folder, data_folder):
	'''Returns a dict recording the contig ID each GI ID maps to.
    Input
    * folder: path to folder containing gff3 files with only the *gene* lines
    * data_folder: path to the output folder, where gi_reads folder will be created

    Output
    * gi_reads: {filename:{SRSXXX.XXXX:scaffoldID, ...}, ...}    '''
	print 'Processing GI to Contig mapping'
	files = os.listdir(folder)
	gi_rid = {}
	for i in files:
		name = re.sub('[\t\r\n]', '', i.split('.')[0]).strip()
		gi_rid[name] = {}
		foo = open(folder + '/' + i)
		foo = foo.readlines()
		for j in foo:
			j_split = j.split('\t')
			if len(j_split) == 9:
				id = re.sub('[\n\t\r]', '', j_split[8].split(';')[0].strip()).split('ID=')[1].split('-')[0]
				if not id in gi_rid[name]:
					gi_rid[name][id] = re.sub('[\n\t\r]', '', j_split[0].strip()) # writes the gi_reads folder to data_folder containing all the mappings of# gi IDs :contig IDs
	write_dict(gi_rid, 'folder','gi_rids', data_folder)
	return gi_rid


def id_abundance(folder, data_folder):
	'''Returns a dict containing contig IDs and their abundances (normalized FPKM) in the HMP samples
	    Input:
	    * folder: folderpath containing the *samtools idxstats* results for the HMP samples
	    * data_folder: folderpath for the output folder where reads_abund_norm will be created
	
	    Normalization: 200 bases are subtracted from the length as an approximation to counter the edge-read effect 

	    Output:
	    * reads: {filename:{ContigID: normalized_abundance, ...}, ...}'''
	print 'Processing Contigs: Abundance'
	files = os.listdir(folder)
	reads = {}
	for i in files:
		name = re.sub('[\n\r\t]', '', i.split('_')[1].split('.')[0])
		foo = open(folder + '/' + i)
		foo = foo.readlines()
		reads[name] = {}
		for j in range(len(foo)-1):
			j_split = foo[j].split('\t')
			try:
				reads[name][re.sub('[\n\t\r]', '', j_split[0].strip())] = float(re.sub('[\n\t\r]', '', j_split[2].strip())) / (float(re.sub('[\t\n\r]', '', j_split[1].strip()))- 200.00)
			except:
				pdb.set_trace()			

# normalization happens here # creates a folder under data_folder that contains all the normalized# abundances for the HMP samples contigs (NOT GENES)
	write_dict(reads, 'folder','reads_abund_norm', data_folder)
	return reads


def write_dict(foo_dict, typ, tmp, data_folder):
	''' Utility function to write dicts to output folder
	    Input:
	    * foo_dict: Dict that needs to be written; can be (A) {key:[],...} or (B) {key:{nest_key:[]}}
	    * tmp: if 'hmgc_gi','hmgc_files', 'prev_hmgc' then treats as A, else as B
	    * data_folder: output folder path where files are written'''
	if typ == 'file':#hmgc_gi' or tmp == 'hmgc_files' or tmp == 'prev_hmgc' or tmp == 'prev_gc':
        	foo = open(data_folder + '/' + tmp + '.txt', 'w')
		for i in foo_dict:
			foo.writelines([i + '\t' + str.join('\t', foo_dict[i]) + '\n'])
		foo.close()
	elif typ == 'folder':
		print tmp
		try:
			os.mkdir(data_folder+'/'+tmp)
		except:
			pass
		for fname in foo_dict:
			foo = open(data_folder+'/'+tmp+'/'+fname+'.txt', 'w')
			for key in foo_dict[fname]:
				foo.writelines([key + '\t' + str(foo_dict[fname][key]) + '\n'])
			foo.close()
	print 'Files written\t' + data_folder + '/' + tmp
