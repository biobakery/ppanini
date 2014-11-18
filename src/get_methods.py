import os
import sys
import re
import pdb
import utils
import numpy

def get_id_abund(folder, data_folder):
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
				reads[name][re.sub('[\n\t\r]', '', j_split[0].strip())] = float(re.sub('[\n\t\r]', '', j_split[2].strip())) / (float(re.sub('[\t\n\r]', '', j_split[1].strip()))- 100.00)
			except:
				pdb.set_trace()			
	utils.write_dict(reads, 'folder','id_abund', data_folder)


def get_prev_gc(hmgc_gi, data_folder):
	hmgc_files = {}
	for gc in hmgc_gi:
		no_files = list(set([gi.split('.')[0].strip() for gi in hmgc_gi[gc]]))
		if not gc.split('.')[0] in no_files:
			no_files += [gc.split('.')[0].strip()]
		hmgc_files[gc] = no_files
	len_hmgc = [len(hmgc_files[i]) for i in hmgc_files]
	thold = int(0.1*max(len_hmgc))
	prev_gc = {}
	for gc in hmgc_files:
		if len(hmgc_files[gc]) >= thold:
			prev_gc[gc] = [str(len(hmgc_files[gc]))]
	utils.write_dict(hmgc_files, 'file','gc_files', data_folder)
	utils.write_dict(prev_gc, 'file','prev_gc', data_folder)


def get_gc_gi(foo, data_folder):
	foo = open(foo)
	foo = foo.readlines()
	hmgc_gi = {}
	for i in foo:  
		if not i[0] == '#':
			split_i = i.split('\t')
			if not re.sub('[\t\n\r]', '', split_i[9]).strip() in hmgc_gi:
				hmgc_gi[re.sub('[\t\n\r]', '', split_i[9]).strip()] = [re.sub('[\t\n\r]', '', split_i[8]).strip()]
			else:
				hmgc_gi[re.sub('[\t\r\n]', '', split_i[9]).strip()] += [re.sub('[\t\n\r]', '', split_i[8]).strip()] 
	for gc in hmgc_gi:
		if not gc in hmgc_gi[gc]:
			hmgc_gi[gc] += [gc]
	utils.write_dict(hmgc_gi, 'file','gc_gi', data_folder)

def get_gi_id(folder, data_folder):
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
				id = re.sub('[\n\t\r]', '', j_split[8].split(';')[0].strip()).split('ID=')[1].split('-')[0]+'-T1-C'
				if not id in gi_rid[name]:
					gi_rid[name][id] = re.sub('[\n\t\r]', '', j_split[0].strip())
	utils.write_dict(gi_rid, 'folder','gi_id', data_folder)


def get_abund_gc(hmgc_gi, gi_contig, contig_abund, prev_gc, hmgc_files, data_folder):
	abund_gc = {}
	for gc in prev_gc:
		abund_gc[gc] = {}
		gc_file = hmgc_files[gc]
		for f in gc_file:
			abund_gc[gc][f] = []
		for gi in hmgc_gi[gc]:
			abund_gc[gc][gi.split('.')[0].strip()] += [contig_abund[gi.split('.')[0].strip()][gi_contig[gi.split('.')[0].strip()][gi]]]  # gc:fname:[00,00,00]
	mean_abund_gc = {}
	median_abund_gc = {}
	
	for gc in abund_gc:
		mean_abund_gc[gc] = numpy.mean([sum(i) for i in abund_gc[gc].values()])
		median_abund_gc[gc] = numpy.median([sum(i) for i in abund_gc[gc].values()])
	thold = numpy.percentile(mean_abund_gc.values(),10)
	prev_abund_gc = {}
	for gc in mean_abund_gc:
		if mean_abund_gc[gc] >= thold:
			prev_abund_gc[gc] = [len(hmgc_files[gc]), mean_abund_gc[gc], median_abund_gc[gc]]
	
	foo = open(data_folder+'/gc_prevalence_abundance.txt','w')
	for gc in prev_abund_gc:
		tmp_gis = []
		tmp_abd = []
		for gi in hmgc_gi[gc]:
			tmp_gis += [gi]
			tmp_abd += [contig_abund[gi.split('.')[0].strip()][gi_contig[gi.split('.')[0]][gi]]]
		foo.writelines([gc+'\t'+str.join('\t',tmp_gis)+'\n'])
		foo.writelines([gc+'\t'+str.join('\t', [str(i) for i in tmp_abd])+'\n'])
		foo.writelines([gc+'\t'+str.join('\t', [str(i) for i in prev_abund_gc[gc]])+'\n'])
	foo.close()
	return prev_abund_gc

