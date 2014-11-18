import os
import re
import sys
import numpy
import time
import pdb
from hmgi_parser import *


def get_prev_hmgc(hmgc_gi):
    '''
    Returns prevalent hmgc centroids
    Input:
    something here
    '''
    # Returns prevalent hmgc
    # hmgc_files: gc: number of files
    # prev_gc: [gc: gis]
    hmgc_files = {}
    for gc in hmgc_gi:
        no_files = list(set([gi.split('.')[0].strip() for gi in hmgc_gi[gc]]))
	if not gc.split('.')[0] in no_files:
		no_files += [gc.split('.')[0].strip()]
        hmgc_files[gc] = no_files
    len_hmgc = [len(hmgc_files[i]) for i in hmgc_files]
    thold = int(0.1*max(len_hmgc))
    print str(thold) + '\t THRESHOLD'
    prev_gc = {}
    for gc in hmgc_files:
        if len(hmgc_files[gc]) >= thold:
            prev_gc[gc] = hmgc_gi[gc]
    print 'Found ' + str(len(prev_gc)) + ' prevalent centroids'
    return [prev_gc, hmgc_files]


def get_imp_hmgc(hmgc_gi, gi_contig, contig_abund, data_folder):
	# get important centroids
	# Returns a dict of only the prevalent hmgc
	[prev_gc, hmgc_files] = get_prev_hmgc(hmgc_gi)
	write_dict(hmgc_files, 'file','hmgc_files', data_folder)
	write_dict(prev_gc, 'file','prev_gc', data_folder)
	abund_gc = {}
	for gc in prev_gc:
		abund_gc[gc] = {}
		gc_file = hmgc_files[gc]
		for f in gc_file:
			abund_gc[gc][f] = []
		for gi in prev_gc[gc]:
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
		for gi in prev_gc[gc]:
			tmp_gis += [gi]
			tmp_abd += [contig_abund[gi.split('.')[0].strip()][gi_contig[gi.split('.')[0]][gi]]]
		foo.writelines([gc+'\t'+str.join('\t',tmp_gis)+'\n'])
		foo.writelines([gc+'\t'+str.join('\t', [str(i) for i in tmp_abd])+'\n'])
		foo.writelines([gc+'\t'+str.join('\t', [str(i) for i in prev_abund_gc[gc]])+'\n'])
	foo.close()
	return prev_abund_gc

if __name__ == '__main__':
	t = time.time()
	print 'Time started:\t' + str(t)
	gff3_folder = sys.argv[1]
	cluster_file = sys.argv[2]
	bam_folder = sys.argv[3]
	data_folder = sys.argv[4]
	print 'Reading hmgc_gi:\t' + cluster_file
	hmgc_gi = cluster_reads(cluster_file, data_folder)  # gc:gis
	contig_abund = id_abundance(bam_folder, data_folder)  # contig:abun
	gi_contig = gi_reads(gff3_folder, data_folder)  # f:contig:abun
	imp_hmgc = get_imp_hmgc(hmgc_gi, gi_contig, contig_abund, data_folder)
	body_site = cluster_file.split('/')[-1].split('_')[-1].split('.')[0]
	foo = open(data_folder + '/prevalent_abundant_centroids_' + body_site + '.txt', 'w')
	for i in imp_hmgc:
		foo.writelines([i + '\t' + str.join('\t', [str(j) for j in imp_hmgc[i]]) + '\n'])
	foo.close()
	print 'Time elapsed:\t' + str((time.time() - t)/60) + ' mins'
