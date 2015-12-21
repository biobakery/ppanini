#**PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration**

PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g. prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools.

Google Group
 ppanini-users: https://groups.google.com/forum/#!forum/ppanini-users

URL
 http://huttenhower.sph.harvard.edu/ppanini

# **Citation**
  **Systematic approach to prioritization of 'important' microbial genes**, Gholamali Rahnavard, Afrah Shafquat, Bahar Sayoldin, Eric A. Franzosa, Curtis Huttenhower (under preparation)
  	

# **1. SETUP**

## **1.1 REQUIREMENTS**

### **PPANINI**
* [matplotlib](http://matplotlib.org/)
* [Python 2.7.*](https://www.python.org/download/releases/2.7/)
* [Biopython](http://biopython.org/wiki/Download)
* [Numpy 1.9.*](http://www.numpy.org/)

### **PREPPANINI**
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://samtools.sourceforge.net/)
* [usearch](http://www.drive5.com/usearch/) **or** [vsearch](https://github.com/torognes/vsearch)
* [diamond](http://omictools.com/diamond-s8955.html) **or** [usearch](http://www.drive5.com/usearch/)  **or** [rapsearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2/)

## **1.2 INSTALLATION**

To install, execute the following command in your Terminal/Commmand prompt:

```
#!cmd
hg clone http://bitbucket.org/biobakery/ppanini
cd ppanini
python setup.py install
```

The prerequisites for executing this command are: 

* [Mercurial](https://mercurial.selenic.com/wiki/Download)

------------------------------------------------------------------------------------------------------------------------------

# **2. Running PPANINI**

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


## **2.1 INPUTS**

* ``-i or --input_table`` : Gene Abundance Table containing annotated gene abundance values in CPM or counts per million
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

* ``--output_folder``: folder containing all the output files
* ``--gene_catalog``: File containing the entire genes catalog for the metagenomic niche (**REQUIRED** if uc file not provided)
* ``--uc``:  [Optional] File containing the clustering information for all the genes in input file (**REQUIRED** if gene_catalog not provided)
* ``--basename``: name prefix for all intermediate output files produced
* ``--log_level``: level of debugging information to be provided; Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
* ``--threads``: Number of threads to be used while clustering
* ``--usearch``: Runs USEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/usearch]
* ``--vsearch``: Runs VSEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/vsearch]
* ``--tshld_abund``: Percentile threshold used to prioritize genes. Default value 75th percentile of the gene abundance observed.
* ``--tshld_prev``: Prevalence cut-off used to prioritize genes. Default value 1/10 samples i.e. val - 2*Standard Error(distribution) > 0.1
* ``--quad``: The quadrant of genes to prioritize {1: High Abundance, Low Prevalence; 2: High Abundance, High Prevalence; 3: Low Abundance, High Prevalence, 4: Low Abundance, Low Prevalence}
* ``--bypass_abund_prev``: To bypass the calculation of important genes


## **2.2 OUTPUT**

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

--------------------------------------------------------------------------------------------------------------

# **3. PREPPANINI: Creating a PPANINI table**

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

# **4. VISUALIZATION**

## **4.1 Metagenomic vs. Genomic Priority plots**

```
#!cmd

usage: ppanini_visualizer [-h] -i INPUT_TABLE
                             [--original-table ORIGINAL_TABLE]
                             [--bypass-cloud] [--prev PREV] [--abund ABUND]
                             [-m MAPPER] [--write-mapper] [--zorder ZORDER]
                             [--hexplot] [--bypass-priority]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TABLE, --input_table INPUT_TABLE
                        Gene abundance table with metadata
  --original-table ORIGINAL_TABLE
                        Gene abundance table with metadata
  --bypass-cloud        To draw Abundance Prevalence Cloud
  --prev PREV           Graph will be prevalence across centroids
  --abund ABUND         Graph will be mean abundance across centroids
  -m MAPPER, --mapper MAPPER
                        GO to UniRef mapper
  --write_mapper        Gene to GO table written
  --zorder ZORDER       Zorder [1,2,3,4] [Old, UniRef, UniRef/GO, NA]
  --hexplot             Plot HEXBIN
  --bypass-priority     Generates Metagenome vs. Genome Priority plots
```

## **4.2 Histograms for metagenome hits**

```
#!cmd
usage: ppanini_plot_metagenome_genome [-h] -i INPUT_FILE [--bypass-parse]
                                 [--parse-only]
                                 [--metagenome-fasta METAGENOME_FASTA]
                                 [--bypass-hist] [--bypass-scatter]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        Gene Genomes blast results
  --bypass-parse        Input file is parsed
  --parse-only          To only parse
  --metagenome-fasta METAGENOME_FASTA
                        Metagenome FASTA file
  --bypass-hist         Generates Histogram
  --bypass-scatter      Generates Scatterplot
```

## **4.3 GraPhlAn plots**

```
#!cmd
usage: ppanini_plot_genome_hits [-h] -i INPUT_FILE --map MAP [--bypass_scatter]
                           [--bypass-stats] [--bypass-graphlan-rings]
                           [--pangenome-size PANGENOME_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        Gene Genomes blast results parsed**
  --map MAP             Gene to GO mapper from ppanini_visualizer
  --bypass-scatter      Scatter plot for genomes
  --bypass-stats        Write stats for genome gene hits
  --bypass-graphlan-rings
                        Generates graphlan rings file
  --pangenome-size PANGENOME_SIZE
                        Pangenome size mapping file
```

------------------------------------------------------------------------------------------------------------------------

# **5. TOOLS**

## **5.1 Normalization**


```
#!cmd
Usage: ppanini_normalize_table <input_table> > <normalized_table>
```

## **5.2 Join tables**


```
#!cmd
usage: ppanini_join_tables <table1> <table2> ... > merged_table.txt
```

## **5.3 Centroids extraction from gene abundance table**


```
#!cmd
Usage: ppanini_imp_centroids_prabXtract <imp_centroids_list> <centroids_abundance_matrix_file> > <imp_centroids_abundance_matrix_file>
```

##**5.4 Centroids extraction from gene catalog fasta**


```
#!cmd
Usage: ppanini_imp_centroids_extracter <imp_centroids_list> <fasta_file> <imp_centroids_fasta_file>
```

## **5.5 Create mapper file for PREPPANINI**


```
#!cmd
usage: ppanini_create_mapper [-h] [--assemblies ASSEMBLIES] [--reads READS]
                        [--gff3s GFF3S] [--sams SAMS] [--bams BAMS]
                        [--abund ABUND] [--annot ANNOT] [--faas FAAS]
                        [--fnas FNAS] [--niche NICHE] --samples SAMPLES
                        [-o OUTPUT_TABLE]

optional arguments:
  -h, --help            show this help message and exit
  --assemblies ASSEMBLIES
                        FOLDER containing ASSEMBLIES
  --reads READS         Folder containing READS files
  --gff3s GFF3S         GFF3 Folder
  --sams SAMS           SAMS Folder
  --bams BAMS           SAMS Folder
  --abund ABUND         SAMS Folder
  --annot ANNOT         SAMS Folder
  --faas FAAS           SAMS Folder
  --fnas FNAS           SAMS Folder
  --niche NICHE         GFF3 Folder
  --samples SAMPLES     GFF3 Folder
  -o OUTPUT_TABLE, --output-table OUTPUT_TABLE
                        Gene Table to write
```

## **5.6 Write mapper**

```
#!cmd
usage: ppanini_write_mapper <uniref_ids> <map uniref_go_ids> > <uniref_ids_go_select>
```
# **License**

This software is licensed under the MIT license.

Copyright (c) 2015 Gholamali Rahnavard, Afrah Shafquat, Eric A. Franzosa, Curtis Huttenhower.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------------------------------------------------

# **Maintained by**
	
[Gholamali Rahnavard](mailto:rahnavar@hsph.harvard.edu) and 
[Afrah Shafquat](mailto:shafquat@hsph.harvard.edu)

