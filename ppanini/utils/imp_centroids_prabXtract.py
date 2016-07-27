import os
import sys
import re
import pdb

'''Takes the first col of first file, and extracts gene abundance rows from the second file'''
def main():
	help = ['-h', '--h', '--help']
	if sys.argv[1] in help:
		print 'Usage: python '+sys.argv[0]+' <imp_centroids_list> <centroids_abundance_matrix_file> > <imp_centroids_abundance_matrix_file>'
		sys.exit()
		
	centroids_list = []
	with open(sys.argv[1]) as centroids:
		for line in centroids:
			if not line.startswith('#'):
				centroids_list += [line.split('\t')[0].strip()]
	
	with open(sys.argv[2]) as centroids_fasta:
		for line in centroids_fasta:
			if line.startswith('#'):
				print line.strip()
			elif line.split('\t')[0].strip() in centroids_list:
					print line.strip()
if __name__ == '__main__':
	main()