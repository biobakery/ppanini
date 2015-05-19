import os
import sys
import re
import pdb

'''Takes the first col of first file, and extracts FASTA sequences from the second file'''

if __name__ == '__main__':
	help = ['-h', '--h', '--help']
	if sys.argv[1] in help:
		print 'Usage: python '+sys.argv[0]+' <imp_centroids_list> <fasta_file> > <imp_centroids_fasta_file>'
		sys.exit()
		
	centroids_list = []
	with open(sys.argv[1]) as centroids:
		for line in centroids:
			if not line.startswith('#'):
				centroids_list += [line.split('\t')[0].strip()]
	
	with open(sys.argv[2]) as centroids_fasta:
		check = False
		genes_written = 0
		for line in centroids_fasta:
			if line.startswith('>'):
				if re.sub('[\r\t\n>]','', line).split(' ')[0].strip() in centroids_list:
					genes_written +=1
					check = True
					print line.strip()
				else:
					check = False
					if genes_written == len(centroids_list):
						break
			elif check:
				print line.strip()