import unittest
import importlib
import re
#import tempfile
import os
import logging
import ppanini
import subprocess

from ppanini import config
#from ppanini.tests import var
from ppanini.demo import var

demo1_faa = var.demo1_faa
demo1_fna = var.demo1_fna
demo1_gff3 = var.demo1_gff3
demo1_test = var.demo1_test

class TestPPANINIBasicFunctions(unittest.TestCase):
	"""Test the functions found in ppanini.py"""
	
	def test_searchPrograms(self):
		try:
			print config.usearch
			subprocess.call([config.usearch, "--version"])
		except:
			try:
				subprocess.call([config.vsearch, "--version"])
			except:
				self.fail('No Clustering program found\nLooked for USEARCH and VSEARCH.')

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
