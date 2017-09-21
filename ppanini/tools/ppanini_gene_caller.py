import os
import sys
import re
import argparse
import logging
import shutil

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
        "-u", "--uniref", 
        #dest = 'uniref',
        required = True,
        help="UniRe90 database",
        )
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
    parser.add_argument(
        "--one-contig",
        dest ="one_contig", 
        help="If there is only contig file for all samples, then this option should eb provided", 
        action="store_true")
    
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    config.threads = args.threads
    config.resume = args.resume
    config.one_contig = args.one_contig

    # Set the basename of the output files if specified as an option
    if args.output_basename:
        config.file_basename=args.output_basename
    else:
        # Determine the basename of the input file to use as output file basename
        input_file_basename=os.path.basename(args.fastq)
        # Remove gzip extension if present
        if re.search('.gz$',input_file_basename):
            input_file_basename='.'.join(input_file_basename.split('.')[:-1])
        # Remove input file extension if present
        if '.' in input_file_basename:
            input_file_basename='.'.join(input_file_basename.split('.')[:-1])
    
        config.file_basename=input_file_basename
    config.output_folder = args.output
    config.temp_dir= config.output_folder+'/'+config.file_basename+'_temp'
    
    #Steps:   
    
    #make a directory or outputs
    utilities.make_directory(config.output_folder)
    utilities.make_directory(config.temp_dir)
    #Steps
    #Commands to generate gene families abundance from assemblies and sequences files
    
    # append sample file name to the contig names
    
    
    #make a directory or outputs
    utilities.make_directory(config.temp_dir) 
    
    # if each sample has its own contig
    if not args.one_contig:
        new_contig_file = utilities.append_filename2cotignames(args.contig)
        
        # make directory for prodigal output
        temp_dir = config.temp_dir
        config.temp_dir = temp_dir+'/prodigal_output/'
        utilities.make_directory(config.temp_dir) 
        genes_file_gff, genes_file_fna, genes_file_faa = utilities.genecall(new_contig_file)
    
    # if only one contig is used then no changing in names is needed     
    else:
        new_contig_file = args.contig 
        # check if there is already a prodigal output
        if not os.path.isfile(config.output_folder+"/prodigal.gff"): # if there is no prodigal output   
            # gene call using prodigal
            genes_file_gff, genes_file_fna, genes_file_faa = utilities.genecall(new_contig_file)
        else:
            genes_file_gff = config.output_folder+ '/prodigal.gff'
            #genes_file_fna = config.output_folder+ '/prodigal.fna'
            genes_file_faa = config.output_folder + '/prodigal.faa'
        
    
    # make directory for bowtie2 output
    temp_dir = config.temp_dir
    config.temp_dir = temp_dir+'/bowtie2_output/'
    utilities.make_directory(config.temp_dir)
    
    # make index database using bowtie2-build
    index_name = utilities.index(new_contig_file)
    
    # reads alignment using Bowtie2
    alignment_file = utilities.alignment(args.fastq, index_name)
    
    # make directory for featureCounts abundance output
    config.temp_dir = temp_dir+'/featureCounts_output/'
    utilities.make_directory(config.temp_dir)
    
    # gene abundance using featureCounts
    abundance_file = utilities.abundance(genes_file_gff, alignment_file)
    
    # Run diamond
    alignment_file = utilities.diamond_alignment(genes_file_faa, args.uniref )
    
    
    # Infer abundance for sufficient hits to  uniref90 and no_hits
    hits_map, no_hits_map = utilities.Infer_aligmnets(alignment_file, config.temp_dir)
    
    # wirite mape file to the ouput
    with open(config.output_folder+'/'+config.file_basename+'_hits_map.txt' , 'w') as foo:
        foo.writelines([str.join('\t', hits_map)+'\n']) #header

    
    # select sequence for insufficient hits
    no_hits_genes_faa = utilities.select_sequnces(genes_file, no_hits_map, output_name = config.file_basename+'_no_hits.faa')
    
    # select sequence for sufficient hits 
    hits_genes_faa = utilities.select_sequnces(genes_file, hits_map, output_name = config.file_basename+'_hits.faa')
    
    
    # move the three main output under main output folder from temp files
    # if there is more than one contig and the prodigal outputs ahvn't been produces (first sample run)
    shutil.move(abundance_file, config.output_folder+'/'+os.path.basename(os.path.normpath(abundance_file)))
    if not os.path.isfile(config.output_folder+"/prodigal.gff"):
        if config.one_contig:
            
            shutil.move(genes_file_gff, config.output_folder+'/prodigal.gff')
            shutil.move(genes_file_faa, config.output_folder+'/prodigal.faa')
            print ("Three main output files for ppanini_press are written in: \n%s\n%s\n%s")% (config.output_folder+'/'+os.path.basename(os.path.normpath(abundance_file)),
                                                config.output_folder+'/prodigal.gff', 
                                                config.output_folder+'/prodigal.faa')
        else:
            shutil.move(genes_file_gff, config.output_folder+'/'+os.path.basename(os.path.normpath(genes_file_gff)))
            shutil.move(genes_file_faa, config.output_folder+'/'+os.path.basename(os.path.normpath(genes_file_faa)))
            print ("Three main output files for ppanini_press are written in: \n%s\n%s\n%s")% (config.output_folder+'/'+os.path.basename(os.path.normpath(abundance_file)),
                                                config.output_folder+'/'+os.path.basename(os.path.normpath(genes_file_gff)), 
                                                config.output_folder+'/'+os.path.basename(os.path.normpath(genes_file_faa)))
        
        
if __name__ == '__main__':
	main()
