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
        "-i", "--genes-path",
        dest = 'gene', 
        required = True,
        help="Prodigal outputs for all samples which includes gff and faa files",
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
    
    config.temp_dir= args.output
    config.threads = args.threads
    #Steps:   
    
    #make a directory or outputs
    utilities.make_directory(config.temp_dir)
    
    # Concatenate all FAA files from prodigal outputs
    temp_out_files=[]
    for gene_file in os.listdir(args.gene):          
        # only use FAA files
        if gene_file.endswith('.faa'):
            temp_out_files.append(args.gene+'/'+gene_file)
    genes_file = utilities.name_temp_file('genes.faa')
    print temp_out_files, genes_file
    utilities.execute_command("cat",temp_out_files,temp_out_files,[genes_file],
        genes_file)
    print genes_file
    #alignment_file = config.temp_dir + '/genes.uniref90hits'

    # Run diamond
    alignment_file = utilities.diamond_alignment(genes_file, args.uniref )
    
    
    # Infer abundance for sufficient hits to  uniref90 and no_hits
    hits_map, no_hits_map = utilities.Infer_aligmnets(alignment_file, config.temp_dir)
    
    
    # select sequence for sufficient hits and insufficient hits
    no_hits_genes_faa = utilities.seletc_sequnces(alignment_file, no_hits_map, output_name = 'no_hits.faa')

    # Cluster no sufficient hits using CD-Hit
    cluster_gene_file, cluster_alignments = utilities.cluster_genes(no_hits_genes_faa)
    
    
    # Generate mapping file for clusters to genes (with no sufficient hit to UniRef90)
    #ppanini_cluster2genes -i ${infer_output}/no_hits_reads.clust90.clstr --output ${infer_output}/cd_hit_clust_temp
    mapping_cluster = utilities.mapping_clusters_genes(cluster_gene_file)
    
    # Join gene families
    gene_families_table = gene2genefamilies(tables, mapping_cluster, mapping_uniref)


if __name__ == '__main__':
	main()
