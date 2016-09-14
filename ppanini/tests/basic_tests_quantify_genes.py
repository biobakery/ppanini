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


class TestQuanitfyGenesBasicFunctions(unittest.TestCase):
	"""Test the functions found in quantify_genes.py"""

	#generate_abundance_viabwt2(assembly_x_withpath, reads_x, sample, out):
	#generate_abundance_viasam(assembly_x_sam_withpath, sample, out):
	#generate_abundance_viabam(assembly_x_bam_withpath, sample, out):
	#read_abundance_tables(mapper, norm_flag):

	def test_quantify_genes(self):
		pass