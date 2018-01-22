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
        "-o", "--output", 
        default=None,
        required = True,
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
        #metavar="<" + str(config.threads) + ">", 
        type=int,
        default=config.threads)
    parser.add_argument(
        "--scale",
        dest= 'scale', 
        help="scale the abundance table\n",
        choices=["rpk","count"], 
        default='rpk')
    parser.add_argument(
        "--memory",
        dest= 'cd_hit_memory', 
        help="memory for -M option in CD-Hit \n",
        default=config.cd_hit_memory)
    
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    
    config.temp_dir= args.output+'/temp'
    config.output_folder= args.output
    config.threads = args.threads
    config.resume = args.resume
    config.cd_hit_memory = args.cd_hit_memory
    #Steps:   
    
    #make a directory or outputs
    utilities.make_directory(config.output_folder)
    utilities.make_directory(config.temp_dir)
    
    # Concatenate all FAA files from prodigal outputs
    temp_no_hits_name = []
    temp_no_hits_faa = []
    temp_hits_name = []
    temp_hits_faa = []
    temp_uniref_gene_file = []
    for gene_file in os.listdir(args.gene_path+'/no_hits/'):          
        if gene_file.endswith('no_hits.faa'):
            temp_no_hits_faa.append(args.gene_path+'/no_hits/'+gene_file)
        elif gene_file.endswith('no_hits.txt'):
            temp_no_hits_name.append(args.gene_path+'/no_hits/'+gene_file)
    for gene_file in os.listdir(args.gene_path+'/hits/'):  
        if gene_file.endswith('hits.faa'):
            temp_hits_faa.append(args.gene_path+'/hits/'+gene_file)
        elif gene_file.endswith('hits.txt'):
            temp_hits_name.append(args.gene_path+'/hits/'+gene_file)
    for map_file in os.listdir(args.gene_path+'/hits/'):  
        if map_file.endswith('uniref_gene_map.txt'):
            temp_uniref_gene_file.append(args.gene_path+'/hits/'+map_file)
    no_hits_faa = utilities.name_temp_file('no_hits_genes.faa')
    no_hits_name = utilities.name_temp_file('no_hits_name.txt')
    hits_faa = utilities.name_temp_file('hits_genes.faa')
    hits_name = utilities.name_temp_file('hits_name.txt')
    uniref_gene_file = utilities.name_temp_file('uniref_gene_map.txt')
    
    utilities.execute_command("cat", temp_no_hits_faa, temp_no_hits_faa, [no_hits_faa], no_hits_faa)    
    utilities.execute_command("cat", temp_no_hits_name, temp_no_hits_name, [no_hits_name], no_hits_name)
    utilities.execute_command("cat", temp_hits_faa, temp_hits_faa, [hits_faa], hits_faa)    
    utilities.execute_command("cat", temp_hits_name, temp_hits_name, [hits_name], hits_name)
    utilities.execute_command("cat", temp_hits_name, temp_hits_name, [hits_name], hits_name)
    utilities.execute_command("cat", temp_uniref_gene_file, temp_uniref_gene_file, [uniref_gene_file], uniref_gene_file)

    # Cluster no sufficient hits using CD-Hit
    cluster_gene_file, cluster_alignments = utilities.cluster_genes(no_hits_faa)
    
    
    # Generate mapping file for clusters to genes (with no sufficient hit to UniRef90)
    #ppanini_cluster2genes -i ${infer_output}/no_hits_reads.clust90.clstr --output ${infer_output}/cd_hit_clust_temp
    mapping_cluster = utilities.mapping_clusters_genes(cluster_gene_file)
    
        
    # Join gene families and write them to the output directory as gene_families_table.txt
    gene_families_table = utilities.gene2genefamilies(args.gene_path, mapping_cluster, uniref_gene_file, args.scale)
    
    # move the the gene families output table under main output folder from temp files
    
    shutil.move(gene_families_table, config.output_folder+'/'+os.path.basename(os.path.normpath(gene_families_table)))

    print ("The gene families output table for ppanin is written in: \n%s")% (config.output_folder+'/'+os.path.basename(os.path.normpath(gene_families_table)))


if __name__ == '__main__':
	main()
