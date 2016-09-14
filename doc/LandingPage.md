# PPANINI User Manual #

PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g. prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools..

If you use the PPANINI software, please cite our manuscript:** Systematic approach to prioritization of 'important' microbial genes**, Gholamali Rahnavard, Afrah Shafquat, Bahar Sayoldin, Eric A. Franzosa, Curtis Huttenhower (under preparation)
  	


## Contents ##

* [PREPPANINI: Creating a PPANINI table](#markdown-header-preppanini-creating-a-ppanini-table)
    * [Requirements](#markdown-header-requirements)
    * [How to run](#markdown-header-how-to-run)
* [Visualization](#markdown-header-visualization)
    * [Metagenomic vs. Genomic Priority plots](#markdown-header-metagenomic-vs-genomic-priority-plots)
    * [Histograms for metagenome hits](#markdown-header-histograms-for-metagenome-hits)
    * [GraPhlAn plots](#markdown-header-graphlan-plots)
* [Tools](#markdown-header-tools)
    * [Normalization](#markdown-header-normalization)
    * [Join tables](#markdown-header-join-tables)
    * [Centroids extraction from gene abundance table](#markdown-header-centroids-extraction-from-gene-abundance-table)
    * [Create mapper file for PREPPANINI](#markdown-header-create-mapper-file-for-preppanini)
    * [Write mapper](#markdown-header-write-mapper)

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

## Visualization ##

### Metagenomic vs. Genomic Priority plots ###

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

### Histograms for metagenome hits ###

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

## GraPhlAn plots ##

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

## Tools ##

### Normalization ###


```
#!cmd
Usage: ppanini_normalize_table <input_table> > <normalized_table>
```

### Join tables ###


```
#!cmd
usage: ppanini_join_tables <table1> <table2> ... > merged_table.txt
```

### Centroids extraction from gene abundance table ###


```
#!cmd
Usage: ppanini_imp_centroids_prabXtract <imp_centroids_list> <centroids_abundance_matrix_file> > <imp_centroids_abundance_matrix_file>
```

### Centroids extraction from gene catalog fasta ###


```
#!cmd
Usage: ppanini_imp_centroids_extracter <imp_centroids_list> <fasta_file> <imp_centroids_fasta_file>
```

### Create mapper file for PREPPANINI ###


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

### Write mapper ###

```
#!cmd
usage: ppanini_write_mapper <uniref_ids> <map uniref_go_ids> > <uniref_ids_go_select>
```
