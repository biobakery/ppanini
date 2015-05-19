Shafquat, Afrah

shafquat@hsph.harvard.edu

November 18, 2014

#**PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration**


```
#!python

usage: ppanini.py [-h] -i INPUT_TABLE [-o OUTPUT_FOLDER] [--uc UC]
                  [--usearch USEARCH] [--vsearch VSEARCH]
                  [--basename BASENAME] [--log_level LOG_LEVEL]
                  [--threads THREADS] [--tshld_abund TSHLD_ABUND]
                  [--tshld_prev TSHLD_PREV]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TABLE, --input_table INPUT_TABLE
                        REQUIRED: Gene abundance table with metadata
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Folder containing results
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
```


## **INPUTS**

* ``-i or --input_table`` : Gene Abundance Table containing annotated gene abundance values in CPM or counts per million
* * Such tables can be obtained using (i) HUMAnN2, (ii) preppanini.py or (iii) manually creating the table using samtools (idxstats) etc.
* * See the mock gene table for an example. naan/input/mock_gene_table.tsv

```
#!text
#SAMPLES SAMPLE_X  SAMPLE_Y
geneID_XYZ|UniRef90_XYZ  0.09 0.00
geneID_MNO|UniRef90_unknown  0.00 0.09
```

* * **Metadata**: *Required*

```
#!text
#FAAS	location_fasta_SAMPLE_X	location_fasta_SAMPLE_Y
#SAMPLES SAMPLE_X	SAMPLE_Y
```

* * **Metadata**: *Optional*

```
#!text
#NICHES	N/A	NICHE_SAMPLE_X	NICHE_SAMPLE_Y
```

* * * Each NICHE corresponds to the type of sample i.e. Human Stool, Skin, Soil, Rainforest etc. 
* * * This data is used to calculate the alpha- and beta- prevalence of the gene centroids i.e. prevalence within a specific niche and/or prevalence across different niches
* * * In absence of niche data, only alpha-prevalence is calculated.

* ``--output_folder``: folder containing all the output files
* ``--uc``:  [Optional] File containing the clustering information for all the genes in input file
* ``--basename``: name prefix for all intermediate output files produced
* ``--log_level``: level of debugging information to be provided; Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]
* ``--threads``: Number of threads to be used while clustering
* ``--usearch``: Runs USEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/usearch]
* ``--vsearch``: Runs VSEARCH for clustering genes using the path provided, including the name. E.g. [/n/usr/bin/vsearch]
* ``--tshld_abund``: Percentile threshold used to prioritize genes. Default value 75th percentile of the gene abundance observed.
* ``--tshld_prev``: Prevalence cut-off used to prioritize genes. Default value 1/10 samples i.e. val - 2*Standard Error(distribution) > 0.1

## **OUTPUT**

Returns a list of "important" uncharacterized genes.

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


##METHODOLOGY

1. Cluster centroids according to annotated UniRef90 IDs
2. Cluster centroids annnotated UniRef90_unknown, using 90% sequence homology (USEARCH)
3. Add all genes within the same sample that belong to a specific centroid (UniRef from STEP1 or o/w from STEP2) across samples
4. Take the mean of gene centroid abundance for each gene centroid (where centroid is present; Criteria for presence= >0)
5. Take the prevalence of each gene centroid as #no.samples the gene centroid is present in. (Criteria for presence= >0)
6. If (mean gene abundance value > percentile(distribution, TSHLD_ABUND)) AND (gene_prevalene value -2*StandardError(distribution)> TSHLD_PREV); Gene is considered IMPORTANT
7. Output the important genes in the results folder specified by user
8. Dump all the intermediate files in the ``tmp`` folder.

**For NICHE-SPECIFIC ANALYSIS**:
* Instead of mean gene abundance, we use max(mean abundance) of gene observed across niches. Percentile is taken within all the max-mean abundances observed.
* Instead of gene prevalence, we use alpha_prevalence of the gene in each niche, and if the gene satisfies the prevalence criteria described above in ANY of the niches, the gene is considered important. Distribution is of gene alpha prevalences observed per niche.

#APPENDIX

##RELEVANT DEFINITIONS

* **Alpha-prevalence**: Prevalence of gene centroid in a specific niche of samples (i.e. HumanSkin, HumanStool, DesertSoil etc.) 
* **Beta-prevalence**: Median of gene centroids's alpha_prevalence across niches.

* **Important genes**
* * Genes that are abundant in >=10th percentile of mean abundance across samples
* * Genes that are prevalent in >=10th percentile of observed prevalence
* * * In the case of **niche-specific analysis**, genes that are >=10th percentile of alpha_prevalence observed for each niche for ANY of the niches!
* **Uncharacterized genes**: Genes that have UniRef90 annotations, but are below *** level of characterization acc. to GO Annotations
* **Unannotated genes**: Genes that lack UniRef90 annotations

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
