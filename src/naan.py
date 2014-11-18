import os
import re
import sys
import numpy
import time
import pdb
import utils
import argparse
import read_methods
import get_methods

if __name__ == '__main__':
	t = time.time()
	parser = argparse.ArgumentParser()
	parser.add_argument('-r','--reads_abund', help='Folder for abundance data for genes [Can be of genes or contigs]')
	parser.add_argument('-c','--cluster_file', help='UCLUST .uc file that describes the clustering of centroids to gene IDs')
	parser.add_argument('-m','--gi_id', nargs = '?' , help='Folder for gene ID to contig ID mapping [Required if contig abundance specified in --reads-abund]')
	parser.add_argument('-f','--fasta', help='Fasta file containing the centroid sequences')
	parser.add_argument('-o','--output_folder', help='Fasta file containing the centroid sequences')
	args = parser.parse_args()
		
	print 'Build Tree'
	try:
		os.mkdir(args.output_folder)
	except:
		pass
	try:
		os.mkdir(args.output_folder+'/prevabund_id')
	except:
		pass
	try:
		os.mkdir(args.output_folder+'/u90_annot')
	except:
		pass
	try:
		os.mkdir(args.output_folder+'/u50_annot')
	except:
		pass

	get_methods.get_gc_gi(args.cluster_file, args.output_folder+'/prevabund_id')
	get_methods.get_id_abund(args.reads_abund, args.output_folder+'/prevabund_id')
	get_methods.get_gi_id(args.gi_id, args.output_folder+'/prevabund_id')
	
	print 'gc_gi, id_abund, gi_id map completed ' + str((time.time()-t)/60) + ' mins elapsed'	

	gc_gi = read_methods.read_gc_file(args.output_folder+'/prevabund_id/gc_gi.txt')
	get_methods.get_prev_gc(gc_gi, args.output_folder+'/prevabund_id')

	gi_id = read_methods.read_gi_folder(args.output_folder+'/prevabund_id/gi_id')
	id_abund = read_methods.read_abund_folder(args.output_folder+'/prevabund_id/id_abund')
	prev_gc = read_methods.read_prev_gc_file(args.output_folder+'/prevabund_id/prev_gc.txt')
	gc_files = read_methods.read_gc_file(args.output_folder+'/prevabund_id/gc_files.txt')
		
	print 'Prevalent centroids determined: ' + str((time.time()-t)/60) + ' mins elapsed'	

	abund_hmgc = get_methods.get_abund_gc(gc_gi, gi_id, id_abund, prev_gc, gc_files, args.output_folder+'/prevabund_id/')

	foo = open(args.output_folder+'/prevabund_id/prevalent_abundant_centroids.txt', 'w')
	for i in abund_hmgc:
		foo.writelines([i + '\t' + str.join('\t', [str(j) for j in abund_hmgc[i]]) + '\n'])
	foo.close()
	
	print 'Abundant centroids determined: '+str((time.time()-t)/60)+' mins elapsed'
