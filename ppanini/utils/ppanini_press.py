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
        "-i", "--gene-path",
        dest = 'gene', 
        required = True,
        help="Prodigal output",
        )
    parser.add_argument( 
        "-u", "--uniref90-db", 
        dest = 'uniref90',
        required = True,
        help="UniRe90 database",
        )
    parser.add_argument( 
        "-u", "--uniref-db", 
        dest = 'uniref',
        required = True,
        help="UniRe90 database",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for outputs",
        )
    parser.add_argument(
        "-r","--resume", 
        help="bypass commands if the output files exist\n", 
        action="store_true",
        default=config.resume)
    
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    
    config.temp_dir= args.output
    
    #Steps:   
    
    #make a directory or outputs
    config.temp_dir +='/diamond_output/'
    utilities.make_directory(config.temp_dir)
    
    # Concatenate all FAA files from prodigal outputs
    for gene_file in os.listdir(config.gene):          
        # only use FAA files
        if gene.endswith('.faa'):
            temp_out_files.append(temp_out_file)
 
    utilities.execute_command("cat",temp_out_files,temp_out_files,[alignment_file],
        alignment_file)

    # Run diamond
    utilities.diamond_alignment(alignment_file,args.uniref, unaligned_reads_file_fasta)
    #diamond blastp --quiet --query $FILE --db /n/huttenhower_lab/data/humann2_databases/uniref_annotated/uniref90/v1.1_uniref90/uniref90_annotated.1.1.dmnd --threads 2 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --out ${diamond_output}/${sample}.uniref90hits &
    
    # Infer abundance for sufficient hits to  uniref90 and no_hits
    
    # Cluster no sufficient hits using CD-Hit
    
    # Generate mapping file for clusters to genes (with no sufficient hit to UniRef90)
    
    # Join gene families
    
    new_contig_file = utilities.append_filename2cotignames(args.contig)
    
    # make directory for prodigal output
    config.temp_dir +='/prodigal_output/'
    utilities.make_directory(config.temp_dir) 
    


if __name__ == '__main__':
	main()
