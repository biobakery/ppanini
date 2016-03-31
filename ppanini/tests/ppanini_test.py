import os
import sys
import unittest
import importlib
import ppanini

py_version = [2 , 7]

if (sys.version_info[0] != py_version[0] or sys.version_info[1] != py_version[1]):

	sys.exit("CRITICAL ERROR: Python version doesn't match"+\
			 "\nRequired version: Python: "+'.'.join([str(i) for i in py_version])+\
			 "\nCurrent version: Python "+str(sys.version_info))	

try:
	from ppanini import tests
except:
	sys.exit("CRITICAL ERROR: Unable to find the PPANINI python package.\nPlease check your install")

def get_unittests():
	directory_of_tests=os.path.dirname(os.path.abspath(__file__))
	
	print directory_of_tests
	
	basic_suite = unittest.TestLoader().discover(directory_of_tests, pattern='basic_tests_*.py')
	advanced_suite = unittest.TestLoader().discover(directory_of_tests, pattern='advanced_tests_*.py')

	return unittest.TestSuite([basic_suite,advanced_suite])

def main():
	full_suite = get_unittests()
	unittest.TextTestRunner(verbosity=2).run(full_suite)