import os
import sys
import re
import argparse

import logging

from .. import quantify_genes
from .. import utilities
from .. import config

logger = logging.getLogger(__name__)
def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser( 
        description="PPANINI gene caller", 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        "-i", "--contig", 
        required = True,
        help="contigs file (fna)",
        )
    parser.add_argument( 
        "-f", "--fastq", 
        required = True,
        help="reads file (fastq)",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for outputs",
        )
    parser.add_argument(
        "--output-basename",
        help="the basename for the output files\n[DEFAULT: " +
        "input file basename]",
        default=config.file_basename,
        metavar="<sample_name>")
    parser.add_argument(
        "-r","--resume", 
        help="bypass commands if the output files exist\n", 
        action="store_true",
        default=config.resume)
    
    args = parser.parse_args()
    return args


def main():
	args = get_args()
	# Set the basename of the output files if specified as an option
	if args.output_basename:
	    config.file_basename=args.output_basename
	else:
	    # Determine the basename of the input file to use as output file basename
	    input_file_basename=os.path.basename(args.contig)
	    # Remove gzip extension if present
	    if re.search('.gz$',input_file_basename):
	        input_file_basename='.'.join(input_file_basename.split('.')[:-1])
	    # Remove input file extension if present
	    if '.' in input_file_basename:
	        input_file_basename='.'.join(input_file_basename.split('.')[:-1])
	
	    config.file_basename=input_file_basename
        
	#Steps
	#Commands to generate gene families abundance from assemblies and sequences files
	
	# append sample file name to the contig names
	
	# contig as input
	new_contig_file = utilities.append_filename2cotignames(args.contig)
	
	# gene call using prodigal
	genes_file_gff, genes_file_fna, genes_file_faa = utilities.genecall(new_contig_file)
	
	# make index database using bowti2-build
	index_name = index(args.fastq)
	
	# gene abundance using Bowtie2
	alignment_file = utilities.alignment(args.fastq, index_name)
	
	# Run featureCounts
	
	
	# find gene clusters mapped to the same uniref90 using diamond
	#uniref = /n/huttenhower_lab/data/humann2_databases/uniref_annotated/uniref90/v1.1_uniref90/uniref90_annotated.1.1.dmnd
	#diamond_alignment(genes_file_faa,uniref, unaligned_reads_file_fasta)
	
	'''
	
	
	* run for all samples (on-by-one):  sh ~/ppanini_stuff/scripts/run_diamond.sh 
	
	$ python /n/huttenhower_lab/tools/ppanini/ppanini/utils/infer_abundance.py  SRSXXXXX.uniref90hits --min-percid .9 --min-qcover .8 --min-scover .8 —output SRSXXXXX  --json
	* for many: sh ppanini_stuff/scripts/run_infer_abund.sh  
	$ python /n/huttenhower_lab/tools/ppanini/ppanini/utils/fasta_select.py -i hmp_sub_nares.faa -f infer_output/no_hits.txt —output no_hits_reads.faa
	$ python /n/huttenhower_lab/tools/ppanini/ppanini/utils/fasta_select.py -i hmp_sub_nares.faa -f infer_output/hits.txt  —output hits_reads.faa
	
	$ source new-modules.sh
	$ module load cd-hit/4.6.4-fasrc02
	$ cd-hit -c .9 -aL .8 -G 0 -T 8 -i no_hits_reads.faa -o no_hits_reads.clust90
	$ python /n/huttenhower_lab/tools/ppanini/ppanini/utils/ppanini_cluster2genes.py -i no_hits_reads.clust90.clstr --output cd_hit_clust_temp —json
	
	$ module load bowtie2/2.3.1-fasrc01
	$ bowtie2-build -f renamed_contigs_SRS015051.fna  renamed_contigs_SRS015051_bowtie2_index_db 
	* for all samples $ sh ~/ppanini_stuff/scripts/mkbowtie2_dbs.sh 
	$ bowtie2 -q -p 8 -x SRS015051_bowtie2_index_db -U SRS015051.fastq.gz -S SRS015051.sam
	* for all samples sh ~/ppanini_stuff/scripts/align_bowtie2.sh 
	
	$ source new-modules.sh 
	$ module load subread/1.5.1-fasrc01 
	$ featureCounts -T 8 -g ID -t CDS  -a ../../prodigal_output/hmp_sub_nares/renamed_contigs_SRS015051.gff -o counts.txt SRS015051.sam 
	* many samples: $ sh ~/ppanini_stuff/scripts/run_featureCounts.sh 
	$ python ~/Documents/Hutlab/ppanini/ppanini/utils/ppanini_join_tables.py -i hmp_sub_nares/ -o hmp_sub_nares.tsv
'''

if __name__ == '__main__':
	main()
