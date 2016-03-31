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


		def test_get_centroids_fromUCLUST(self, gis_dm):
			"""Tests the function get_centroids_fromUCLUST"""

			config.uclust_file = test_config.demo_ppanini_clusters
			cluster_dict = pp.get_centroids_fromUCLUST(gis_dm.keys(), config)
			
			#Type of output
			self.assertTrue(type(cluster_dict)==dict) 
			#For this dataset, 91 centroids
			self.assertEqual(len(cluster_dict), 91) 
			#All centroids should contain themselves
			self.assertEqual(sum(1 for i in cluster_dict if i in cluster_dict[i]), 91) 
		
		def test_get_centroids(self, uniref_dm, gi_dm):
			"""Tests the function get_centroids"""
			
			def test_get_centroids_table(self,gc_dm):
				"""Tests the function get_centroids_table"""
				
				[norm_data_matrix, centroids_list]= pp.get_centroids_table(gc_dm, metadata, config)
				#Number of centroids and rows in dm
				self.assertEqual(norm_data_matrix.shape[0], len(centroids_list)) 

			gc_dm = pp.get_centroids(uniref_dm, gi_dm, config)
			#gc_dm contains all the UniRef90 clusters
			self.assertEqual(sum([1 for i in uniref_dm if i in gc_dm]), len(uniref_dm))
			#gc_dm is <= uniref90 + gi centroids
			self.assertLessEqual(len(gc_dm), len(uniref_dm)+len(gi_dm)) 

			test_get_centroids_table(self, gc_dm)


		test_get_centroids_fromUCLUST(self, gis_dm)
		test_get_centroids(self, uniref_dm, gis_dm)
		shutil.rmtree(config.temp_folder)
	
	def test_get_clusters(self):
		"""Tests the function get_clusters"""

		usearch_path = config.usearch
		vsearch_path = config.vsearch
		config.usearch =''
		config.vsearch =''
		with self.assertRaises(Exception): ##Test if usearch/vsearch not present; it throws an exception
			pp.get_clusters(config)
