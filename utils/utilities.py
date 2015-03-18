import os
import sys
import pdb
import re
import numpy


def create_folders(list_folders):
	'''Creates the list of folders if they dont exist already
	Input: list_folders = [List of folders to be created]

	Output: Folders created'''

	for fname in list_folders:
		try:
			os.mkdir(fname)
		except:
			pass