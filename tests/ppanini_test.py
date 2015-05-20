import os
import sys
import unittest
import importlib


py_version = [2 , 7]

if (sys.version[0] != py_version[0] or sys.version[1] != py_version[1]):
	sys.exit("CRITICAL ERROR: Python version doesn't match"+\
			 "\nRequired version: Python: "+'.'.join([str(i) for i in py_version])+\
			 "\nCurrent version: Python "+'.'.join([str(i) for i in sys.version[:2]]))

try:
	from ppanini import tests
except:
	sys.exit("CRITICAL ERROR: Unable to find the PPANINI python package.\nPlease check your install")

python_version()


def main():
	directory_of_tests=os.path.dirname(os.path.abspath(__file__))

	basic_suite = unittest.TestLoader().discover(directory_of_tests, pattern='basic_tests_*.py')
    advanced_suite = unittest.TestLoader().discover(directory_of_tests, pattern='advanced_tests_*.py')
    full_suite = unittest.TestSuite([basic_suite,advanced_suite])
   
    return full_suite 