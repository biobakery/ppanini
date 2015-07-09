Shafquat, Afrah

[shafquat@hsph.harvard.edu](mailto:shafquat@hsph.harvard.edu)

July 10, 2015

#**PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration**

PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g. prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools.

## **REQUIREMENTS***

* PPANINI
* * matplotlib
* *  python 2.7
* * Biopython
* * Numpy 1.6.1??

* PREPPANINI
* bowtie2
* samtools
* usearch/vsearch
* diamond/usearch/rapsearch2

## **INSTALLATION**

To install, execute the following command in your Terminal/Commmand prompt:

```
#!cmd
hg clone http://bitbucket.org/biobakery/ppanini
```

The prerequisites for executing this command are: 

* [Mercurial](https://mercurial.selenic.com/wiki/Download)

Once cloned, run the following command:

```
#!cmd
export PYTHONPATH=$PYTHONPATH:<INSERT PATH to PPANINI HERE>
```


## Running PPANINI

```
#!python

usage: ppanini.py [-h] -i INPUT_TABLE [-o OUTPUT_FOLDER]
                  [--gene_catalog GENE_CATALOG] [--uc UC] [--usearch USEARCH]
                  [--vsearch VSEARCH] [--basename BASENAME]
                  [--log_level LOG_LEVEL] [--threads THREADS]
                  [--tshld_abund TSHLD_ABUND] [--tshld_prev TSHLD_PREV]
                  [--quad QUAD] [--bypass_prev_abund]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TABLE, --input_table INPUT_TABLE
                        REQUIRED: Gene abundance table with metadata
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Folder containing results
  --gene_catalog GENE_CATALOG
                        GENE CATALOG
  --uc UC               UCLUST file containg centroids and clustered genes
  --usearch USEARCH     Path to USEARCH
  --vsearch VSEARCH     Path to VSEARCH
  --basename BASENAME   BASENAME for all the output files
  --log_level LOG_LEVEL
                        Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
  --threads THREADS     Number of threads
  --tshld_abund TSHLD_ABUND
                        [X] Percentile Cutoff for Abundance; Default=75th
  --tshld_prev TSHLD_PREV
                        Threshold: val-2*SE > tshld_prev
  --quad QUAD           Quadrant analysis
  --bypass_prev_abund   Bypass quantifying abundance and prevalence
```


## **INPUTS**

* ``-i or --input_table`` : Gene Abundance Table containing annotated gene abundance values in CPM or counts per million
* * Such tables can be obtained using (i) HUMAnN2, (ii) preppanini.py or (iii) manually creating the table using samtools (idxstats) etc.
* * See the mock gene table for an example. ppanini/input/mock_gene_table.tsv

```
#!text
#SAMPLES SAMPLE_X  SAMPLE_Y
geneID_XYZ|UniRef90_XYZ  0.09 0.00
geneID_MNO|UniRef90_unknown  0.00 0.09
```

* * **Metadata**: *Optional*

```
#!text
#NICHES	N/A	NICHE_SAMPLE_X	NICHE_SAMPLE_Y
#SAMPLES SAMPLE_X SAMPLE_Y
```

* * * Each NICHE corresponds to the type of sample i.e. Human Stool, Skin, Soil, Rainforest etc. 
* * * This data is used to calculate the alpha- and beta- prevalence of the gene centroids i.e. prevalence within a specific niche and/or prevalence across different niches
* * * In absence of niche data, only alpha-prevalence is calculated.

* ``--output_folder``: folder containing all the output files
* ``--gene_catalog``: File containing the entire genes catalog for the metagenomic niche
* ``--uc``:  [Optional] File containing the clustering information for all the genes in input file
* ``--basename``: name prefix for all intermediate output files produced
* ``--log_level``: level of debugging information to be provided; Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
* ``--threads``: Number of threads to be used while clustering
* ``--usearch``: Runs USEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/usearch]
* ``--vsearch``: Runs VSEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/vsearch]
* ``--tshld_abund``: Percentile threshold used to prioritize genes. Default value 75th percentile of the gene abundance observed.
* ``--tshld_prev``: Prevalence cut-off used to prioritize genes. Default value 1/10 samples i.e. val - 2*Standard Error(distribution) > 0.1
* ``--quad``: The quadrant of genes to prioritize {1: High Abundance, Low Prevalence; 2: High Abundance, High Prevalence; 3: Low Abundance, High Prevalence, 4: Low Abundance, Low Prevalence}
* ``--bypass_abund_prev``: To bypass the calculation of important genes


## **OUTPUT**

Returns a list of "important" genes.

```
#!text
#GENEID  PREVALENCE MEAN_ABUNDANCE
geneID_XYZ  0.25   0.05
```

**For NICHE-SPECIFIC ANALYSIS**:

```
#!text
#GENEID   MEAN_ABUNDANCE  ALPHA_PREVALENCE_NICHEX ALPHA_PREVALENCE_NICHEY BETA_PREVALENCE
geneID_XYZ  0.05  0.35  0.50   0.42 
```

=========================================================================

#PREPROCESSING SCRIPTS

```
#!cmd

usage: preppanini.py [-h] -m MAPPER_FILE [--basename BASENAME]
                     [--bypass_abundance] [--bypass_annotation]
                     [--bypass_clust] [--bypass_write_table]
                     [--usearch USEARCH] [--vsearch VSEARCH]
                     [--diamond DIAMOND] [--rapsearch RAPSEARCH]
                     [--threads THREADS] [--uniref90 UNIREF90]
                     [--uniref50 UNIREF50] [--to_normalize]
                     [--log_level LOG_LEVEL]

optional arguments:
  -h, --help            show this help message and exit
  -m MAPPER_FILE, --mapper_file MAPPER_FILE
                        Mapper file containing paths to data
  --basename BASENAME   BASENAME for all the output files
  --bypass_abundance    Bypass quantifying abundance
  --bypass_annotation   Bypass annotating genes
  --bypass_clust        Bypass annotating genes
  --bypass_write_table  Bypass writing table
  --usearch USEARCH     Path to USEARCH
  --vsearch VSEARCH     Path to VSEARCH
  --diamond DIAMOND     Path to DIAMOND
  --rapsearch RAPSEARCH
                        Path to RAPSEARCH
  --threads THREADS     Number of threads
  --uniref90 UNIREF90   UniRef90 INDEX file
  --uniref50 UNIREF50   UniRef50 INDEX file
  --to_normalize        Default HUMAnN2 table; if sam-idxstats table; enable
  --log_level LOG_LEVEL
                        Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
```

* ``--mapper_file``: REQUIRED. Contains all information about which FAAS (Gene FAA format) corresponds to which SAMPLE and to which ABUNDANCE FILE; etc. Please see the demo input for more examples. All the possible options are: ['SAMS', 'BAMS', 'NICHE', 'GFF3S', 'ABUNDANCE_TABLES', 'ANNOTATION', 'READS', 'CONTIG_ASSEMBLIES', 'FAAS', 'FNAS', 'SAMPLES']
* ``--basename``: name prefix for all intermediate output files produced
* ``--bypass_abundance``: bool flag to bypass running samtools/bowtie2 to calculate the gene abundances
* ``--bypass_annotation``: bool flag to bypass running sequence homology search for the gene sequences against UniRef90/50.
* ``--bypass_write_table``: bool flag to bypass writing ppanini input table
* ``--usearch``: Runs USEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/usearch]
* ``--vsearch``: Runs VSEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/vsearch]
* ``--diamond``: Runs DIAMOND for sequence-homology search against UniRef90 using the path provided, including the name. E.g. [/n/usr/bin/diamond]
* ``--rapsearch``: Runs RAPSEARCH for sequence-homology search against UniRef90 using the path provided, including the name. E.g. [/n/usr/bin/rapsearch2]
* ``--threads``: Number of threads to be used while clustering or sequence-homology search tool.
* ``--uniref90``: UniRef90 database index filename for appropriate for the search tool being used. For information on how to create these indices, please refer to DIAMOND/RAPSEARCH user manual.
* ``--uniref50``: UniRef50 database index filename for appropriate for the search tool being used. For information on how to create these indices, please refer to DIAMOND/RAPSEARCH user manual.
* ``--to-normalize``: bool flag [True when sequence length normalization is required. E.g. when samtools is being used to translate Hits to FPKM or CPM]; [False when normalization has already been performed and an abundance table is being read. E.g. HUMAnN2 gene abundance table]
* ``--log_level``: level of debugging information to be provided; Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]


========================================================================================

Alternate names

* Annotation of Novel Genes via Integrated Re-assignment Approach (ANGIRA)