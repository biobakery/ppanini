#**PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration**
----

 * Download the PPANINI software ([ppanini.tar.gz](https://bitbucket.org/biobakery/ppanini/downloads/biobakery-ppanini-0.6.0.tar.gz)) then follow the [steps to install and run](#markdown-header-getting-started-with-ppanini).

 * For additional information, please see the [PPANINI User Manual](http://bitbucket.org/biobakery/ppanini).

 * Please direct questions to the [PPANINI google group](https://groups.google.com/forum/#!forum/ppanini-users) (subscribe to receive PPANINI news).

 * If you use the PPANINI software, please cite our manuscript: **Systematic approach to prioritization of 'important' microbial genes**, Gholamali Rahnavard, Afrah Shafquat, Bahar Sayoldin, Eric A. Franzosa, Curtis Huttenhower (under preparation)
  	

----

PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g. prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools.

# **PPANINI Workflow**
![ppanini_workflow.png](https://bitbucket.org/repo/rnag7z/images/994033213-ppanini_workflow.png)

------------------------------------------------------------------------------------------------------------------------------

# **SETUP**

## **REQUIREMENTS**

### **PPANINI**
* [matplotlib](http://matplotlib.org/)
* [Python 2.7.*](https://www.python.org/download/releases/2.7/)
* [Biopython](http://biopython.org/wiki/Download)
* [Numpy 1.9.*](http://www.numpy.org/)

------------------------------------------------------------------------------------------------------------------------------
## **INSTALLATION**
Download and unpack the PPANINI software: [ppanini.tar.gz](https://bitbucket.org/biobakery/ppanini/downloads/biobakery-ppanini-0.6.1.tar)
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

The prerequisites for executing this command are: 

[Mercurial](https://mercurial.selenic.com/wiki/Download)
```

------------------------------------------------------------------------------------------------------------------------------

# **Getting Started with PPANINI**

```
#!python

usage: ppanini [-h] -i INPUT_TABLE [-o OUTPUT_FOLDER]
                  [--gene-catalog GENE_CATALOG] [--uc UC] [--usearch USEARCH]
                  [--vsearch VSEARCH] [--basename BASENAME]
                  [--log-level LOG_LEVEL] [--threads THREADS]
                  [--tshld-abund TSHLD_ABUND] [--tshld-prev TSHLD_PREV]
                  [--beta BETA]

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
```


## **INPUTS**

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

Each NICHE corresponds to the type of sample i.e. Human Stool, Skin, Soil, Rainforest etc. 
This data is used to calculate the alpha- and beta- prevalence of the gene centroids i.e. prevalence within a specific niche and/or prevalence across different niches
In absence of niche data, only alpha-prevalence is calculated.

* ``--output-folder``: folder containing all the output files
* ``--gene-catalog``: File containing the entire genes catalog for the metagenomic niche (**REQUIRED** if uc file not provided)
* ``--uc``:  [Optional] File containing the clustering information for all the genes in input file (**REQUIRED** if gene_catalog not provided)
* ``--basename``: name prefix for all intermediate output files produced
* ``--log-level``: level of debugging information to be provided; Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
* ``--threads``: Number of threads to be used while clustering
* ``--usearch``: Runs USEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/usearch]
* ``--vsearch``: Runs VSEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/vsearch]
* ``--tshld-abund``: Percentile threshold used to prioritize genes. Default value 75th percentile of the gene abundance observed.
* ``--tshld-prev``: Prevalence cut-off used to prioritize genes. Default value 1/10 samples i.e. val - 2*Standard Error(distribution) > 0.1
* ``--quad``: The quadrant of genes to prioritize {1: High Abundance, Low Prevalence; 2: High Abundance, High Prevalence; 3: Low Abundance, High Prevalence, 4: Low Abundance, Low Prevalence}
* ``--bypass-abund-prev``: To bypass the calculation of important genes


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
------------------------------------------------------------------------------------------------------------------------

# **Maintained by**
	
[Gholamali Rahnavard](mailto:rahnavar@hsph.harvard.edu) and 
[Afrah Shafquat](mailto:shafquat@hsph.harvard.edu)