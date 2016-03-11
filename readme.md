#**PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration**
PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g. prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools.
----

 * Download the PPANINI software ([ppanini.tar.gz](https://bitbucket.org/biobakery/ppanini/downloads/biobakery-ppanini-0.6.0.tar.gz)) then follow the [steps to install and run](#markdown-header-getting-started-with-ppanini).

 * For additional information, please see the [PPANINI User Manual](http://huttenhower.sph.harvard.edu/ppanini/manual).

 * Please direct questions to the [PPANINI google group](https://groups.google.com/forum/#!forum/ppanini-users) (subscribe to receive PPANINI news).

 * If you use the PPANINI software, please cite our manuscript: **Systematic approach to prioritization of 'important' microbial genes**, Gholamali Rahnavard, Afrah Shafquat, Bahar Sayoldin, Eric A. Franzosa, Curtis Huttenhower (under preparation)
  	

----


## Contents ##
* [PPANINI](#markdown-header-ppanini)
    * [PPANINI Workflow](#markdown-header-ppanini-workflow)
    * [Requirements](#markdown-header-requirements)
    * [Installation](#markdown-header-installation)
* [Getting Started with PPANINI](#markdown-header-getting-started-with-ppanini) 
    * [Options](#markdown-header-options) 
    * [Input](#markdown-header-input)
    * [Output](#markdown-header-output)  
* [PREPPANINI: Creating a PPANINI table](#markdown-header-preppanini-creating-a-ppanini-table)
    * [Requirements](#markdown-header-requirements)
    * [How to run](#markdown-header-how-to-run)

------------------------------------------------------------------------------------------------------------------------------
# PPANINI #
## PPANINI Workflow ##
![ppanini_workflow.png](https://bitbucket.org/repo/rnag7z/images/3722133193-ppanini_workflow.png)

## REQUIREMENTS ##
* [matplotlib](http://matplotlib.org/)
* [Python 2.7.*](https://www.python.org/download/releases/2.7/)
* [Biopython](http://biopython.org/wiki/Download)
* [Numpy 1.9.*](http://www.numpy.org/)

## INSTALLATION ##
Download and unpack the PPANINI software: 
 ([ppanini.tar.gz](https://bitbucket.org/biobakery/ppanini/downloads/biobakery-ppanini-0.6.0.tar.gz)) 
```
$ tar -zxvf ppanini.tar.gz
$ cd ppanini
```
Install the PPANINI software
```
$ python setup.py install
```

or download the PPANINI software from the repository:

```
#!cmd
$ hg clone http://bitbucket.org/biobakery/ppanini
$ cd ppanini
$ python setup.py install
```

The prerequisites for executing this command are: 

* [Mercurial](https://mercurial.selenic.com/wiki/Download)


------------------------------------------------------------------------------------------------------------------------------

# Getting Started with PPANINI #
## Options##

```
#!python
usage: ppanini.py [-h] -i INPUT_TABLE [-o OUTPUT_FOLDER]
                  [--gene-catalog GENE_CATALOG] [--uc UC] [--usearch USEARCH]
                  [--vsearch VSEARCH] [--basename BASENAME]
                  [--log-level LOG_LEVEL] [--threads THREADS]
                  [--tshld-abund TSHLD_ABUND] [--tshld-prev TSHLD_PREV]
                  [--beta BETA] [--bypass-clustering]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TABLE, --input_table INPUT_TABLE
                        REQUIRED: Gene abundance table with metadata
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Folder containing results
  --gene-catalog GENE_CATALOG
                        GENE CATALOG
  --uc UC               UCLUST file containg centroids and clustered genes
  --usearch USEARCH     Path to USEARCH
  --vsearch VSEARCH     Path to VSEARCH
  --basename BASENAME   BASENAME for all the output files
  --log-level LOG_LEVEL
                        Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
  --threads THREADS     Number of threads
  --tshld-abund TSHLD_ABUND
                        [X] Percentile Cutoff for Abundance; Default=75th
  --tshld-prev TSHLD_PREV
                        Percentile cutoff for Prevalence
  --beta BETA           Beta parameter for weights on percentiles
  --bypass-clustering   Bypass clustering
```


## Input ##

* ``-i or --input-table`` : Gene Abundance Table containing annotated gene abundance values in CPM or counts per million
Such tables can be obtained using (i) ``HUMAnN2``, (ii) ``preppanini`` or (iii) manually creating the table using ``samtools`` (idxstats) etc.
See the mock gene table for an example. ppanini/input/mock_gene_table.tsv

```
#!text
#SAMPLES SAMPLE_X  SAMPLE_Y
geneID_XYZ|UniRef90_XYZ  0.09 0.00
geneID_MNO|UniRef90_unknown  0.00 0.09
```

**Metadata**: *Optional*

```
#!text
#NICHES	N/A	NICHE_SAMPLE_X	NICHE_SAMPLE_Y
#SAMPLES SAMPLE_X SAMPLE_Y
```

Each NICHE corresponds to the type of sample e.g. Human Stool, Skin, Leaf, Rainforest etc. 

This data is used to calculate the alpha- and beta- prevalence of the gene centroids i.e. prevalence within a specific niche and/or prevalence across different niches. In absence of niche data, only alpha-prevalence is calculated.

* ``--output-folder``: folder containing all the output files
* ``--gene-catalog``: File containing the entire genes catalog for the metagenomic niche (**REQUIRED** if uc file not provided)
* ``--uc``:  [Optional] File containing the clustering information for all the genes in input file (**REQUIRED** if gene_catalog not provided)
* ``--basename``: name prefix for all intermediate output files produced
* ``--log-level``: level of debugging information to be provided; Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
* ``--threads``: Number of threads to be used while clustering
* ``--usearch``: Runs USEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/usearch]
* ``--vsearch``: Runs VSEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/vsearch]
* ``--bypass-clustering``: To bypass clustering of unannotated genes based on homology of sequences (this option is best used when the clustering information[given via --uc] or the gene catalog, both dont exist, **OR** if you prefer to rank each unannotated gene individually)
* ``--beta``: The parameter used in the methodology to prioritize genes via their metagenomic properties. 
* ``--tshld_prev``: Percentile cut-off for prevalence
* ``--tshld_abund``: Percentile cut-off for abundance



## Output ##

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
------------------------------------------------------------------------------------------------------------------------

## PREPPANINI: Creating a PPANINI table ##
#### Requirements ###
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://samtools.sourceforge.net/)
* [usearch](http://www.drive5.com/usearch/) **or** [vsearch](https://github.com/torognes/vsearch)
* [diamond](http://omictools.com/diamond-s8955.html) **or** [usearch](http://www.drive5.com/usearch/)  **or** [rapsearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2/)
### How to run
```
#!cmd

usage: preppanini [-h] -m MAPPER_FILE [--basename BASENAME]
                     [--bypass-abundance] [--bypass-annotation]
                     [--bypass-clust] [--bypass-write-table]
                     [--usearch USEARCH] [--vsearch VSEARCH]
                     [--diamond DIAMOND] [--rapsearch RAPSEARCH]
                     [--threads THREADS] [--uniref90 UNIREF90]
                     [--to-normalize] [--log-level LOG_LEVEL]

optional arguments:
  -h, --help            show this help message and exit
  -m MAPPER_FILE, --mapper-file MAPPER_FILE
                        Mapper file containing paths to data
  --basename BASENAME   BASENAME for all the output files
  --bypass-abundance    Bypass quantifying abundance
  --bypass-annotation   Bypass annotating genes
  --bypass-clust        Bypass annotating genes
  --bypass-write-table  Bypass writing table
  --usearch USEARCH     Path to USEARCH
  --vsearch VSEARCH     Path to VSEARCH
  --diamond DIAMOND     Path to DIAMOND
  --rapsearch RAPSEARCH
                        Path to RAPSEARCH
  --threads THREADS     Number of threads
  --uniref90 UNIREF90   UniRef90 INDEX file
  --to-normalize        Default HUMAnN2 table; if sam-idxstats table; enable
  --log-level LOG_LEVEL
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

------------------------------------------------------------------------------------------------------------------------


# **Maintained by**
	
[Gholamali Rahnavard](mailto:rahnavar@hsph.harvard.edu) and 
[Afrah Shafquat](mailto:shafquat@hsph.harvard.edu)
