import os
import sys
import pdb
import re
import argparse

def clean_make(files, tmp, n, index_file, output_folder, body_site):
	indz = open(index_file)
	indz = indz.readlines()
	tags_list = [int(re.sub('[\n\r\t]','',i))-1 for i in indz] #indices of > in main file
	tags_list.sort()
	foo = open(files,'r') #main file
	foo_lines = foo.readlines()
	if tmp == 'n': # divide main file into several sub files using fixed fasta header no.
		n = int(n)
		ii =[tags_list[i] for i in range(0,len(tags_list), n)]
		for i in range(len(ii)):
			if not ii[i] == ii[-1]:
				name = files.split('/')[-1]
				foo_new = open(output_folder+'/extracted_'+str(ii[i])+'_'+str(ii[i+1])+'_'+name,'w')
				foo_new.writelines(foo_lines[ii[i]:ii[i+1]])
				foo_new.close()
			else:
				name =files.split('/')[-1]
				foo_new = open(output_folder+'/extracted_'+str(ii[i])+'_end_'+name, 'w')
				foo_new.writelines(foo_lines[ii[i]:])
				foo_new.close()
	elif tmp == 'ind': #extract specific fasta sequences
		foo2= open(n)
		foo2 = foo2.readlines()
		inds_foo2 = [int(re.sub('[\r\n\t]','',i))-1 for i in foo2]
		name = files.split('/')[-1]
		foo_new = open(output_folder+'/extracted_'+body_site+'.txt','w')
		for i in range(len(inds_foo2)):
			end = tags_list.index(inds_foo2[i])+1
			if end >= len(tags_list):
				end = len(foo_lines)
			else:
				end = tags_list[end]
			#print str(i)+'\t'+str(end)
			foo_new.writelines(foo_lines[inds_foo2[i]:end])
		foo_new.close()
	print 'Process completed'

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f','--file',help='File that needs to be split')
	parser.add_argument('-t','--type',help='Type of split\nOptions: tmp or n or ind\nn: Split file with the number of fasta headers specified using flag -h 1000\nind: Only extract the specified indicies\n')
	parser.add_argument('-n','--num-fasta',help='If n then the number of fasta headers per file, (ii)  file of selected integers')
	parser.add_argument('-a','--all',help='File containing indices of all headers')
	parser.add_argument('-o','--output',help='Folder to put the files in')
	parser.add_argument('-b','--body-site',help='Body site e.g. stool, anteriornares')

	if len(sys.argv) == 7:
		#n = int(sys.argv[2])
		clean_make(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
	else:
		print('Usage: python %s <filepath> <tmp n or ind> <number_of_fasta_headers> <index_all> <output_folder> <body_site: stool>' %sys.argv[0])
		print('Usage: python %s <filepath> <tmp n or ind> <index_select> <index_all> <output_folder> <body_site>' %sys.argv[0])
