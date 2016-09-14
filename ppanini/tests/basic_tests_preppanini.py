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


class TestPrePPANINIBasicFunctions(unittest.TestCase):
	"""Test the functions found in preppanini.py"""
	def test_preppanini(self):
		pass
