import unittest
import importlib
import re
import tempfile
import os
import logging
import ppanini

from ppanini.tests import test_config
from ppanini.utils import preppanini


demo_faa = test_config.demo_faa
demo_fna = test_config.demo_fna
demo_gff3 = test_config.demo_gff3
#demo1_test = test_config.demo_test


class TestAnnotateGenesBasicFunctions(unittest.TestCase):
	"""Test the functions found in annotate_genes.py"""
	#get_annotations_dict(centroid_annotations, centroid_gis):
	#get_clusters_dict(gene_centroid_clusters_file_path):
	#run_vclust(vsearch_folder, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, perc_id, nprocesses):
	#run_uclust(usearch_folder, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, perc_id, nprocesses):
	#run_rapsearch(query_file, db, out_fname, all_paths, nprocesses):
	#run_diamond(query_file, db, out_fname, all_paths, nprocesses):
	#parse_annotation_table(annotations_file, fasta_sequences, thld_ref):

	def test_annotate_genes(self):
		pass