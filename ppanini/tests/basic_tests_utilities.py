import unittest
import importlib
import re
import tempfile
import os
import logging
import ppanini
import time
import random
from ppanini.tests import test_config
from ppanini import utilities


demo_faa = test_config.demo_faa
demo_fna = test_config.demo_fna
demo_gff3 = test_config.demo_gff3

#demo1_test = test_config.demo_test


class TestUtilitiesBasicFunctions(unittest.TestCase):
	"""Test the functions found in utilities.py"""

	def test_is_present(self):
		self.assertEqual(utilities.is_present(['Soil\tWater\tSoil','DEF'],'#NICHE'),[[],[]])
		self.assertIsNotNone(utilities.is_present(['#NICHE\tSoil\tWater\tSoil','DEF'],'#NICHE'))
		self.assertIsNotNone(utilities.is_present(['#Niche\tSoil\tWater\tSoil','DEF'],'#NICHE'))
		self.assertIsNotNone(utilities.is_present(['#niche\tSoil\tWater\tSoil','DEF'],'#NICHE'))

	def test_write_dict(self):
		pass
	
	def test_create_folders(self):
		pass
	def test_is_protein(self):
		fna_dict=utilities.read_fasta(demo_fna)
		ind = random.choice(range(0, len(fna_dict)-1))
		self.assertFalse(utilities.is_protein(fna_dict.values()[ind]))

		faa_dict=utilities.read_fasta(demo_faa)
		ind2 = random.choice(range(0, len(faa_dict)-1))
		self.assertTrue(utilities.is_protein(faa_dict.values()[ind]))

	def test_write_fasta(self):
		pass

	def test_read_gff3(self):
		pass

	def test_pullgenes_fromcontigs(self):
		test_fna =os.path.join(test_config.data_folder,"demo_test.fna") #translated by function
		test_faa =os.path.join(test_config.data_folder,"demo_test.faa")
		
		faa_dict=utilities.read_fasta(demo_faa) #original

		utilities.pullgenes_fromcontigs(demo_fna, demo_gff3, test_fna, test_faa)
		# print "pulling done"+str(time.time())
		test_fna_dict = utilities.read_fasta(test_fna)
		# print "reading fna done"+str(time.time())
		test_faa_dict = utilities.read_fasta(test_faa)
		# print "reading faa done"+str(time.time())
		ind = random.choice(test_faa_dict.keys())
		ind2 = ind.split('|')[0]+'-T1-C'
		self.assertEqual(test_faa_dict[ind], faa_dict[ind2])
		os.remove(test_fna)
		os.remove(test_faa)

	def test_read_fasta(self):
		faa_dict=utilities.read_fasta(demo_faa)
		self.assertEqual(len(faa_dict), 5509) #number of genes read

	def test_read_ppanini_imp_genes_table(self):
		pass

## Test for all genes being present
##Test for all unannotated genes being present
##Test for clustering containing all genes
##Test for read_fasta and write_fasta consistency
#T#Test for gene length consistency when writing fasta
