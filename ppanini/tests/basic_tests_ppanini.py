import unittest
import importlib
import re
#import tempfile
import os
import logging
import ppanini
import subprocess
import shutil
import numpy
from ppanini import ppanini as pp
from ppanini import config
from ppanini.tests import test_config
from ppanini import utilities

class TestPPANINIBasicFunctions(unittest.TestCase):
	"""Test the functions found in ppanini.py"""
	
	config.basename='ppanini_test'
	config.temp_folder =os.path.join(test_config.data_folder,"temp")
	utilities.create_folders([config.temp_folder])

	def test_read_gene_table(self):
		"""Tests the function read_gene_table"""

		config.input_table= test_config.demo_ppanini_input
		[uniref_dm, gis_dm, metadata] = pp.read_gene_table(config)
		
		#Tests the type of the returned objects
		self.assertTrue(type(uniref_dm)==dict)
		self.assertTrue(type(gis_dm)==dict)
		self.assertTrue(type(metadata)==list)
		# Number of genes is >0
		self.assertTrue(len(uniref_dm)>0 | len(gis_dm)>0) 

		#Tests the contents of the returned data structures
		
		#Number of UniRef50 clusters
		self.assertEqual(len([i for i in uniref_dm if 'UniRef50' in i]), 86) 
		#Shape of the uniref_dm matrix
		self.assertEquals(numpy.array(uniref_dm.values()).shape, (879, 38)) 
		#Shape of the gis_dm matrix
		self.assertEquals(numpy.array(gis_dm.values()).shape, (91, 38)) 


		shutil.rmtree(config.temp_folder)
	