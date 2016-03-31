import unittest
import importlib
import re
import tempfile
import os
import logging
import ppanini

from ppanini.tests import test_config

demo_faa = test_config.demo_faa
demo_fna = test_config.demo_fna
demo_gff3 = test_config.demo_gff3
#demo1_test = test_config.demo_test


class TestPPANINIBasicFunctions(unittest.TestCase):
	"""Test the functions found in preppanini.py"""
	def test_XYZ(self):
		self.assertEqual(1, 1)
	def test_MNO(self):
		self.assertTrue(1)
		#self.assertDictEqual


##Test for Contigs to FAAS
##Test for FNAS to FAAS (if only FNA provided)
##Test for if centroid not in UC file created as self containing cluster?
##Test for join tables
## Test for all genes being present
##Test for all unannotated genes being present
##Test for clustering containing all genes
##Test for read_fasta and write_fasta consistency
#T#Test for gene length consistency when writing fasta




##Test for Contigs to FAAS
##Test for FNAS to FAAS (if only FNA provided)
##Test for if centroid not in UC file created as self containing cluster?
##Test for join tables
## Test for all genes being present
##Test for all unannotated genes being present
##Test for clustering containing all genes
##Test for read_fasta and write_fasta consistency
#T#Test for gene length consistency when writing fasta

# utilities.create_folder
# utilities.is_present
# annotate_genes.uclust, vclust, get_clusters_dict