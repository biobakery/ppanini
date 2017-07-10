import os
import sys
import re
import argparse
import logging

#from .. import quantify_genes
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
    
	utilities.make_directory(config.file_basename)    
	#Steps
	#Commands to generate gene families abundance from assemblies and sequences files
	
	# append sample file name to the contig names
	
	# contig as input
	new_contig_file = utilities.append_filename2cotignames(args.contig)
	print new_contig_file
	# gene call using prodigal
	genes_file_gff, genes_file_fna, genes_file_faa = utilities.genecall(new_contig_file)
	
	# make index database using bowti2-build
	index_name = index(args.fastq)
	
	# reads alignment using Bowtie2
	alignment_file = utilities.alignment(args.fastq, index_name)
	
	# gene abundance using featureCounts
	abundance_file = utilities.abundance(genes_file_gff, alignment_file)
	
	# join gene abundance tables
	

if __name__ == '__main__':
	main()