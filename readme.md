# **PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration**
----

PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g., prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools.


 * For additional information and a quick demo, please see the [PPANINI tutorial](https://github.com/biobakery/biobakery/wiki/ppanini).

 * Please direct questions to the [PPANINI support forum](https://forum.biobakery.org/c/Microbial-community-profiling/PPANINI/25).

 * If you use the PPANINI software, please cite our manuscript: Gholamali Rahnavard, Afrah Shafquat, Himel Mallick, Jason Lloyd-Price, Kevin Bonham, Bahar Sayoldin, Eric A. Franzosa, Curtis Huttenhower, **Identifying important uncharacterized genes using metagenomes and metatranscriptomes**. huttenhower.sph.harvard.edu/ppanini
  	

----


## Contents ##
- [PPANINI](#ppanini)
  - [PPANINI Workflow](#ppanini-workflow)
  - [REQUIREMENTS](#requirements)
  - [INSTALLATION](#installation)
- [Getting Started with PPANINI](#getting-started-with-ppanini)
  - [TEST PPANINI](#test-ppanini)
  - [Options](#options)
  - [Input](#input)
  - [Output](#output)
- [Guides to PPANINI utility scripts](#guides-to-ppanini-utility-scripts)
  - [ppanini_join_tables](#ppanini_join_tables)
  - [ppanini_barplot](#ppanini_barplot)
  - [Convert abundance units](#convert-abundance-units)
  - [Creating a gene table using PPANINI](#creating-a-gene-table-using-ppanini)
    - [Requirements](#requirements-1)
    - [Workflow](#workflow)
    - [How to run](#how-to-run)
      - [Step 1: gene abundance](#step-1-gene-abundance)
      - [Step 2: join tables](#step-2-join-tables)
      - [Step 3: gene families abundance](#step-3-gene-families-abundance)
  - [Creating a gene table using HUMAnN2](#creating-a-gene-table-using-humann2)
  - [Contributions ##](#contributions-)


------------------------------------------------------------------------------------------------------------------------------
# PPANINI #
## PPANINI Workflow ##
![overview_b.png](https://bitbucket.org/repo/rnag7z/images/3094733523-overview_b.png)
## REQUIREMENTS ##
* [Matplotlib](http://matplotlib.org/)
* [Python 2.7.*](https://www.python.org/download/releases/2.7/)
* [Biopython](http://biopython.org/wiki/Download)
* [Numpy 1.9.*](http://www.numpy.org/)
* [Pandas (version >= 0.18.1)](http://pandas.pydata.org/getpandas.html)

## INSTALLATION ##
```
$ pip install ppanini
```
------------------------------------------------------------------------------------------------------------------------------


# Getting Started with PPANINI #
## TEST PPANINI##

To test if PPANINI is running correctly, you may run the following command in the terminal:

```
#!cmd

ppanini_test

```

Which yields:

```
test_create_folders (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_is_present (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_is_protein (basic_tests_utilities.TestUtilitiesBasicFunctions) ... /n/huttenhower_lab/tools/ppanini/depends/lib/Bio/Seq.py:2309: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
  BiopythonWarning)
ok
test_pullgenes_fromcontigs (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_read_fasta (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_read_gff3 (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_read_ppanini_imp_genes_table (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_write_dict (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_write_fasta (basic_tests_utilities.TestUtilitiesBasicFunctions) ... ok
test_annotate_genes (basic_tests_annotate_genes.TestAnnotateGenesBasicFunctions) ... ok
test_read_gene_table (basic_tests_ppanini.TestPPANINIBasicFunctions)
Tests the function read_gene_table ... Gene Table contains 2 metadata lines .
Gene Table contains 998 gene or centroid lines.
ok
test_preppanini (basic_tests_preppanini.TestPrePPANINIBasicFunctions) ... ok
test_quantify_genes (basic_tests_quantify_genes.TestQuanitfyGenesBasicFunctions) ... ok

----------------------------------------------------------------------
Ran 13 tests in 4.101s

OK
```

## Options##

```
#!python
usage: ppanini [-h] -i INPUT_TABLE [-o OUTPUT_FOLDER] [--basename BASENAME]
               [--uniref2go UNIREF2GO] [--log-level LOG_LEVEL]
               [--tshld-abund TSHLD_ABUND] [--tshld-prev TSHLD_PREV]
               [--beta BETA] [--version]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TABLE, --input_table INPUT_TABLE
                        REQUIRED: Gene abundance table with metadata
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Folder containing results
  --basename BASENAME   BASENAME for all the output files
  --uniref2go UNIREF2GO
                        uniref to GO term mapping file
  --log-level LOG_LEVEL
                        Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
  --tshld-abund TSHLD_ABUND
                        [X] Percentile Cutoff for Abundance; Default=75th
  --tshld-prev TSHLD_PREV
                        Percentile cutoff for Prevalence
  --beta BETA           Beta parameter for weights on percentiles
  --version             prints the version
```


## Input ##

* ``-i or --input-table:`` Gene families abundance table containing annotated gene abundance values in RPK or relative abundance.
Such tables can be obtained using [ppanini_gene_caller](#markdown-header-creating-a-gene-table).
See the mock gene table as an example. ppanini/input/mock_gene_table.tsv
* HUMAnN2 genefamilies output also can be used as input for PPANINI

```
#!text
#SAMPLES SAMPLE_1  SAMPLE_2  SAMPLE_3
UniRef90_XYZ  0.09 0.00  0.03
UniRef90_MNP  0.00 0.09  0.0
```

* ``--output-folder``: folder containing all the output files
* ``--basename``: name prefix for all intermediate output files produced
* ``--log-level``: level of debugging information to be provided; Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
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


# Guides to PPANINI utility scripts #

## ppanini_join_tables ##

## ppanini_barplot ##

* **Basic usage:** `$ ppanini_barplot -i1 $PPANINI_GENE_INPUT.txt-i2 $PPANINI_OUTPUT_TABLE.txt`
* `$PPANINI_GENE_INPUT.txt` = a gene family abundance table used by PPANINI (it also in the temp directory in the output 
* `$PPANINI_OUTPUT_TABLE.txt` = the PPANINI output table 
* Run with `-h` to see additional command line options

Produces a stacked bar plot for gene families characterization.

```ppanini_barplot -i1 temp/ppanini_abundance_table.txt -i2 ppanini_table.txt```

```
usage: ppanini_barplot [-h] [-i1 <input table>] [-i2 <input table>]
                       [--summary-table SUMMARY_TABLE] [-o <feature id>]

PPANINI plotting tool

optional arguments:
  -h, --help            show this help message and exit
  -i1 <input table>, --ppanini-input <input table>
                        Gene abundance table
  -i2 <input table>, --ppanini-output <input table>
                        PPANINI output table
  --summary-table SUMMARY_TABLE
                        Summary table
  -o <feature id>, --output <feature id>
                        output plot)
```
![nares_ppanini_barplot.png](https://bitbucket.org/repo/rnag7z/images/559636210-nares_ppanini_barplot.png)

## Convert abundance units ##

Gene families abundance can be converted to relative abundance or “copies per million” with [humann2_renorm_table](https://github.com/biobakery/biobakery/wiki/humann2#32-normalizing-rpks-to-relative-abundance) script, while the unit of abundance doesn't effect the PPANINI results as PPANINI uses the rank of abundances internally. 

------------------------------------------------------------------------------------------------------------------------
## Creating a gene table using PPANINI ##
### Requirements ###
* [Prodigal V2.60](http://prodigal.ornl.gov/)
* [bowtie2 2.2.1](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [diamond 0.8.22](http://omictools.com/diamond-s8955.html) 
* [featureCounts/ subread 1.5.1](http://bioinf.wehi.edu.au/featureCounts/)
### Workflow ###
![ppanini_gene_caller.png](https://bitbucket.org/repo/rnag7z/images/4047935082-ppanini_gene_caller.png)

### How to run ###
There are two main steps to get gene families abundant: 1) Generating gene abundance for each sample, 2) clustering genes, including known genes that they map to UniRef90 and an uncharacterized gene with 90% homology-based similarity.

#### Step 1: gene abundance ####

To get the abundance of genes for a sample, short read file (FASTQ or FASTA) and contig file (FNA) is needed. Please use the same output name for all samples in this step to be used as input for Step 3. For example for all sample use ```-o ppanini_genecaller_output``` and use sample name as basename for the files using ```--output-basename $SAMPLE_ID``` where $SAMPLE_ID holds the name of sample or sample file.

```
#!cmd

usage: ppanini_gene_caller [-h] -i CONTIG -f FASTQ -o OUTPUT
                           [--output-basename <sample_name>] -u UNIREF [-r]
                           [--threads <1>] [--one-contig]

PPANINI gene caller

optional arguments:
  -h, --help            show this help message and exit
  -i CONTIG, --contig CONTIG
                        contigs file (fna)
  -f FASTQ, --fastq FASTQ
                        reads file (fastq)
  -o OUTPUT, --output OUTPUT
                        Path for outputs
  --output-basename <sample_name>
                        the basename for the output files
                        [DEFAULT: input file basename]
  -u UNIREF, --uniref UNIREF
                        UniRe90 database
  -r, --resume          bypass commands if the output files exist
  --threads <1>         number of threads/processes
                        [DEFAULT: 1]
  --one-contig          If there is only one contig file for all samples, then this option should be provided
```

#### Step 2: join tables ####

```
#!cmd
usage: ppanini_join_tables [-h] [-v] -i INPUT -o OUTPUT
                           [--file_name FILE_NAME] [-s]
                           [--mapping-uniref MAPPING_UNIREF]
                           [--mapping-cluster MAPPING_CLUSTER] [-r]
                           [--scale {rpk,count}]

Join gene, pathway, or taxonomy tables

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         additional output is printed
  -i INPUT, --input INPUT
                        the directory of tables
  -o OUTPUT, --output OUTPUT
                        the table to write
  --file_name FILE_NAME
                        only join tables with this string included in the file name
  -s, --search-subdirectories
                        search sub-directories of input folder for files
  --mapping-uniref MAPPING_UNIREF
                        Mapping file: gene to uniref90 file
  --mapping-cluster MAPPING_CLUSTER
                        Mapping file: cluster to genes file
  -r, --resume          bypass commands if the output files exist
  --scale {rpk,count}   scale the abundance table
```

After running panini_gene_caller for all samples, the results can be merged into a table with samples as columns and genes as rows:
1) make a directory (e.g. `tables`) and copy all gene abundance there. for example, you can copy all TXT file under the previouse step to a folder.
```mkdir tables```
```cp ppanini_gene_caller_output/*.txt tables``` 

```
#!cmd
 ppanini_join_tables -i tables/ -o gene_abundance.tsv
```
#### Step 3: gene families abundance ####
To generate gene families, panini_press puts genes that map to the same UniRef90 ID into the same gene family and for genes that don't map sufficiently to a UniRef90 ID, it clusters them with the same criteria that have been used for clustering and mapping into UniRef90 database. For the abundance of gene families, ppanini_press sums up the abundance of all genes belonging to that gene family. 
**Input:** 1) a path to the Step 1's output ```ppanini_gene_caller``` which includes TXT, GFF, and FAA files for all samples. 2) the path to UniRef90 database.

```
#!cmd
usage: ppanini_press [-h] -i GENE_PATH [-o OUTPUT] [-r] [--threads THREADS]
                     [--scale {rpk,count}] [--memory CD_HIT_MEMORY]

PPANINI Press: clusters genes to gene families including annotated genes to UniRef90 and homology-based clustered genes.

optional arguments:
  -h, --help            show this help message and exit
  -i GENE_PATH, --gene-path GENE_PATH
                        a directory path to ppanini_gene_caller outputs which includes txt, gff, and faa files for each sample.
  -o OUTPUT, --output OUTPUT
                        Path for outputs
  -r, --resume          bypass commands if the output files exist
  --threads THREADS     number of threads/processes
                        [DEFAULT: 1]
  --scale {rpk,count}   scale the abundance table
  --memory CD_HIT_MEMORY
                        memory for -M option in CD-Hit 
```

**Input files description:**

The input for `ppanini_press` is the output directory from ppanini_gene_caller:

* **GENE_COUNTS:** featureCounts output files for all samples from the ppanini_gene_caller step. For each sample, we have $SAMPLE_FILENAME.txt under ppanini_gene_caller output.

* **mapping file for sufficiently mapped genes to UniRef90:** These files are under `hits` directory in the ppanini_gene_caller output.

* **GENE_SEQUNCES:** prodigal output files for genes that didn't map sufficiently to UniRef90 in ppanini_gene_caller. These files are under `no_hits` directory. For each sample, we have `$SAMPLE_FILENAME_no_hits.faa` under `no_hits` directory under ppanini_gene_caller output.


** Basic usage:**

```$ ppanini_press -i ppanini_gene_caller_output -o ppanini_press_output --threads 4```

** Gene families table as output:**

The final result of ppanini_press is a table for gene families, where rows are gene families (UniRefs or homology-based clusters), and columns are sample IDs.
To see the table you can use:

```$ less -S ppanini_press_output/gene_families_table.txt ``` 

Which yields (abbreviated):

```
# header                SRS015051      SRS015269      SRS015290      SRS015430      SRS015450      SRS015584      SRS015640      SRS015752      SRS015793      SRS015937      SRS015996
Cluster 0               0              100.103199174  0              0              0              0              0              0              0              0              0
Cluster 10003           0              0              0              0              0              0              0              78.5440613027  0              0              0
Cluster 1001            0              0              0              0              0              0              0              0              0              429.018136335  0
Cluster 10013           0              0              0              0              0              0              0              112.847222222  0              0              0
Cluster 1002            0              0              0              0              0              0              0              0              0              160.725453408  0

.
.
.

UniRef90_Y9BN98         0              0              0              0              0              0              0              0              0              76.1904761905  0
UniRef90_Z3N8E7         0              0              0              0              0              0              0              0              0              50.7246376812  0
UniRef90_Z6R8L8         0              0              0              0              0              0              0              0              0              62.7177700348  0
UniRef90_Z7U7N3         0              0              0              0              0              0              0              34.1207349081  0              0              0
UniRef90_Z8BCZ8         0              0              0              0              0              0              0              0              0              0              80.8080808081
```
This file can be used as input for PPANINI. There are other important files that could be used for further down stream analyses: 

** 1) map_cluster_gene.txt:** this file include a mapping file from artificial cluster names to its belonging genes that they have no sufficient hit to UniRef database. 

** 2) map_uniref_gene.txt:**: this file includes a mapping file from UniRef ID to its belonging genes.

** 3) no_hits.faa:** This file includes gene sequences for genes with no sufficient hits in UniRef database.

** 4) hits.faa:** This file includes gene sequences for genes with sufficient hits in UniRef database.

** 5) genes.faa:** This file includes gene sequences for all genes.

** 6) no_hits_reads.clust90  no_hits_reads.clust90.clstr:** CD-Hit output for genes with no sufficient hits in UniRef database.

** 7) genes.uniref90hits:** Diamond output of mapping genes to UniRef90.

** 8) hits.txt:** mapping UniRef IDs to genes with sufficient hits. 

** 9) no_hits:** genes with no sufficient hits to UniRef90 database.


------------------------------------------------------------------------------------------------------------------------
## Creating a gene table using HUMAnN2 ##
A gene family abundance table can be obtained from [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2) output. HUMAnN2 is a reference based approach which means we get the abundance only for genes that have a hit or reference in the gene database, and there is no information for uncharacterized genes except a row that shows the overall abundance of all reads that they haven't had a hit in the gene reference database. To use HUMAnN2 outputs for all your samples, the following steps need to be performed:

1. Run HUMAnN2 each sample to get gene families output for all your samples. [here](https://github.com/biobakery/biobakery/wiki/humann2#22-running-humann-20-the-basics) is help for how to run HUMAnN2.

2. Merge gene families abundance table for all sample using [humann2_join_tables](https://github.com/biobakery/biobakery/wiki/humann2#42-humann-20-for-multiple-samples).

3. Create unstratified gene family abundance table from gene family abundance table using [humann2_split_stratified_table](https://github.com/biobakery/biobakery/wiki/humann2#52-visualizing-stratified-humann-20-output). Use the `genefamilies_unstratified.tsv` File form the output folder. Gene families from HUMAnN2 outputs are stratified by species, but for PPANINI we need use unstratified gene families.

## Contributions ## 
Thanks go to these wonderful people:

