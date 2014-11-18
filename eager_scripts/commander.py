import sys
import os


def print_files(input_folder, output_folder, db_folder, analysis):
	#print 'ENTERS COMMANDER'
	indexer_files = os.listdir(input_folder)
	#print indexer_files
	if analysis == 'u90':
		id = 0.9
	else:
		id =0.5
	head = ['#! /bin/sh\n', '#BSUB -u shafquat@hsph.harvard.edu\n','#BSUB -q normal_serial\n','#BSUB -g /usearch\n']
	if not 'bsub_'+analysis in os.listdir('.'):
		os.mkdir('./bsub_'+analysis)
	for i in indexer_files:
		#print 'comes'
		name = i
		foo = open('./bsub_'+analysis+'/'+name+'_bsub.bsub','w')
		foo.writelines(head)
		head_ii = ['#BSUB -J '+analysis+'_'+name+'\n', '#BSUB -o out/'+analysis+'_'+name+'.out\n', '#BSUB -e err/'+analysis+'_'+name+'.err\n']
		foo.writelines(head_ii)
		foo.writelines(['\n','module load bio/usearch7.0.1001\n'])
		#foo.writelines(['\nfiles=$(ls '+db_folder+')\nread -a array_files <<< $files\nfor element in "${array_files[@]}"\ndo\nusearch -usearch_local '+input_folder+'/'+name+' -db '+db_folder+'/$element -target_cov 0.8 -query_cov 0.9 -id 0.9 -blast6out '+output_folder+'/'+name+'_$element.m8\ndone\n'])
		foo.writelines(['\nfiles=$(ls '+db_folder+')\nread -a array_files <<< $files\nfor element in "${array_files[@]}"\ndo\nusearch -usearch_global '+input_folder+'/'+name+' -db '+db_folder+'/$element -id '+str(id)+' -blast6out '+output_folder+'/'+name+'_$element.m8\ndone\n'])
		foo.writelines(['\ncat '+output_folder+'/'+name+'_*.m8 > '+output_folder+'/'+name+'_'+analysis+'.m8\n'])
		foo.close()

if __name__ == '__main__':
        if len(sys.argv) == 5:
                input_folder = sys.argv[1]
                output_folder = sys.argv[2]
                db_folder = sys.argv[3]
                analysis = sys.argv[4]
                print_files(input_folder, output_folder, db_folder, analysis)
        else:
                print 'Usage: python %s <input_folder> <output_folder> <db_folder> <analysis>' %sys.argv[0]
