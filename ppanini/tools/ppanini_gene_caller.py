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
        required = True,
        help="Path for outputs",
        )
    parser.add_argument(
        "--output-basename",
        dest = 'output_basename',
        help="the basename for the output files\n[DEFAULT: " +
        "input file basename]",
        default=config.file_basename,
        metavar="<sample_name>")
    parser.add_argument(
        "-r","--resume", 
        help="bypass commands if the output files exist\n", 
        action="store_true",
        default=config.resume)
    parser.add_argument(
        "--threads", 
        help="number of threads/processes\n[DEFAULT: " + str(config.threads) + "]", 
        metavar="<" + str(config.threads) + ">", 
        type=int,
        default=config.threads)
    
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    config.threads = args.threads
    config.resume = args.resume

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
    config.output_folder = args.output
    config.temp_dir= config.output_folder+'/'+config.file_basename+'/'+os.path.basename(os.path.normpath(config.output_folder))+'_temp'
    
    #Steps:   
    
    #make a directory or outputs
    utilities.make_directory(config.output_folder)
    utilities.make_directory(config.temp_dir)
    #Steps
    #Commands to generate gene families abundance from assemblies and sequences files
    
    # append sample file name to the contig names
    
    
    #make a directory or outputs
    utilities.make_directory(config.temp_dir) 
    
    # contig as input
    new_contig_file = utilities.append_filename2cotignames(args.contig)
    
    # make directory for prodigal output
    config.temp_dir +='/prodigal_output/'
    utilities.make_directory(config.temp_dir) 
    
    #config.file_basename = ''
    # gene call using prodigal
    genes_file_gff, genes_file_fna, genes_file_faa = utilities.genecall(new_contig_file)
    
    # make directory for bowtie2 output
    config.temp_dir = args.output+'/bowtie2_output/'
    utilities.make_directory(config.temp_dir)
    
    # make index database using bowtie2-build
    index_name = utilities.index(new_contig_file)
    
    # reads alignment using Bowtie2
    alignment_file = utilities.alignment(args.fastq, index_name)
    
    # make directory for featureCounts abundance output
    config.temp_dir = args.output+'/featureCounts_output/'
    utilities.make_directory(config.temp_dir)
    
    # gene abundance using featureCounts
    abundance_file = utilities.abundance(genes_file_gff, alignment_file)
    
    # move the three main output under main output folder from temp files
    utilities.execute_command("mv", abundance_file, config.output_folder, config.output_folder+'/'+os.path.basename(os.path.normpath(abundance_file)))


if __name__ == '__main__':
	main()
