
#default values for user options
version = '0.7.0'
input_table = ''
output_folder = ''
gene_catalog = '' 
uclust_file = ''
usearch = False
vsearch = False
basename = 'ppanini'
log_level = 'DEBUG'
verbose = 'DEBUG'
nprocesses = 1 
tshld_abund = 75
tshld_prev = 75
beta  = 0.5
bypass_prev_abund = False
tshld = None
centroid_prev_abund = None
niches = []
#all_prevalence = None
#all_mean_abund = None
niche_flag = None
ppanini_niche_score_labels = []




temp_folder = ''
centroids_list = ''
essantial_genes_uniref90_id = None
bypass_clustering = False
genomic_score = False
uniref2go = ''

# uniref formatting
uniref_delimiter="|"
uniref_gene_index=-2
uniref_length_index=-1

# bowtie2 options and threshold
bowtie2_large_index_threshold=4000000000
bowtie2_index_ext_list=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",
    ".rev.1.bt2",".rev.2.bt2"]
bowtie2_large_index_ext=".1.bt2l"
bowtie2_version={
    "flag" : "--version",
    "major" : 2,
    "minor" : 2,
    "line" : 0,
    "column" : 2}

bowtie2_build_opts=["-q"]
bowtie2_align_opts=["--very-sensitive"]
bowtie2_index_name="_bowtie2_index"

featureCounts_opts=["-T 8 -g ID -t CDS"]


# translated alignment options
translated_alignment_choices = ["usearch","rapsearch","diamond", "vsearch"]
translated_alignment_selected = translated_alignment_choices[2]

# prodigal options
prodigal_opts=["-q"]
threads = 8

# file naming
temp_dir=""
unnamed_temp_dir=""
file_basename=""

# diamond options
diamond_database_extension=".dmnd"
diamond_opts_uniref50=["--max-target-seqs",20,"--sensitive","--outfmt",6]
diamond_opts_uniref90=["--max-target-seqs",20,"--outfmt",6]
diamond_cmmd_protein_search="blastp"
diamond_cmmd_nucleotide_search="blastx"
diamond_version={
    "flag" : "--version",
    "major" : 0,
    "minor" : 8,
    "second minor" : 22,
    "line" : 0,
    "column" : 2}

resume =  False
if __name__=='__main__':
	pass
