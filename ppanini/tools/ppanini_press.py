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
        description="PPANINI Press: clusters genes to gene families including annotated genes to UniRef90 and homology-based clustered genes.", 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        "-i", "--gene-path",
        dest = 'gene_path', 
        required = True,
        help="a directory path to ppanini_gene_caller outputs which includes txt, gff, and faa files for each sample.",
        )
    parser.add_argument( 
        "-u", "--uniref", 
        #dest = 'uniref',
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
    parser.add_argument(
        "--threads", 
        help="number of threads/processes\n[DEFAULT: " + str(config.threads) + "]", 
        metavar="<" + str(config.threads) + ">", 
        type=int,
        default=config.threads)
    parser.add_argument(
        "--scale",
        dest= 'scale', 
        help="scale the abundance table\n",
        choices=["rpk","count"], 
        default='rpk')
    
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    
    config.temp_dir= args.output+'/'+os.path.basename(os.path.normpath(args.output))+'_temp'
    config.output_folder= args.output
    config.threads = args.threads
    config.resume = args.resume
    #Steps:   
    
    #make a directory or outputs
    utilities.make_directory(config.output_folder)
    utilities.make_directory(config.temp_dir)
    
    # Concatenate all FAA files from prodigal outputs
    temp_out_files=[]
    for gene_file in os.listdir(args.gene_path):          
        # only use FAA files
        if gene_file.endswith('.faa'):
            temp_out_files.append(args.gene_path+'/'+gene_file)
    genes_file = utilities.name_temp_file('genes.faa')
    utilities.execute_command("cat", temp_out_files, temp_out_files, genes_file)
    
    # Run diamond
    print (genes_file, args.uniref )
    alignment_file = utilities.diamond_alignment(genes_file, args.uniref )
    
    
    # Infer abundance for sufficient hits to  uniref90 and no_hits
    hits_map, no_hits_map = utilities.Infer_aligmnets(alignment_file, config.temp_dir)
    
    
    # select sequence for insufficient hits
    no_hits_genes_faa = utilities.select_sequnces(genes_file, no_hits_map, output_name = 'no_hits.faa')
    
    # select sequence for sufficient hits 
    hits_genes_faa = utilities.select_sequnces(genes_file, hits_map, output_name = 'hits.faa')

    # Cluster no sufficient hits using CD-Hit
    cluster_gene_file, cluster_alignments = utilities.cluster_genes(no_hits_genes_faa)
    
    
    # Generate mapping file for clusters to genes (with no sufficient hit to UniRef90)
    #ppanini_cluster2genes -i ${infer_output}/no_hits_reads.clust90.clstr --output ${infer_output}/cd_hit_clust_temp
    mapping_cluster = utilities.mapping_clusters_genes(cluster_gene_file)
    
        
    # Join gene families and write them to the output directory as gene_families_table.txt
    gene_families_table = utilities.gene2genefamilies(args.gene_path, mapping_cluster, hits_map, args.scale)


if __name__ == '__main__':
	main()
