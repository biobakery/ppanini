Shafquat, Afrah

shafquat@hsph.harvard.edu

November 18, 2014

#**Novel Automated gene Annotation eNumeration (NAAN)**

Prioritization of functional characterization of novel and uncharacterized genes

```
#!python

usage: naan.py [-h] [-i INPUT_TABLE] [-f FASTA_FOLDER] [-u [UCLUST_FOLDER]]
               [-t TYPE] [-o OUTPUT_FOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TABLE, --input_table INPUT_TABLE
                        Gene abundance table with metadata
  -f FASTA_FOLDER, --fasta_folder FASTA_FOLDER
                        Folder containing fasta files
  -u USEARCH_FOLDER, --usearch_folder [USEARCH_FOLDER]
                        Path for USEARCH program
  -t TYPE, --type TYPE  Type of analysis Choose: [gene_table,
                        reads_assemblies]
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Folder containing results

```


## **INPUTS**

* **INPUT_TABLE: Gene Abundance Table**
* * Table containing gene abundance values in FPKM (i.e. #READS/(length of gene-100 nts)) 
* * * Such tables can be obtained using HUMAnN2 or manually creating the table using samtools (idxstats) etc.
* * * See the mock gene table for an example. naan/input/mock_gene_table.tsv

```
#!text
#GENEID ANNOTATION  SAMPLE_X  SAMPLE_Y
geneID_XYZ  UniRef90_XYZ  0.09 0.00
geneID_MNO  UniRef90_unknown  0.00 0.09
```

* * **Metadata**: *Required*

```
#!text
#FASTA	N/A	location_fasta_SAMPLE_X	location_fasta_SAMPLE_Y
#GENEID	ANNOTATION	SAMPLE_X	SAMPLE_Y
```

* * **Metadata**: *Optional*

```
#!text
#NICHE	N/A	NICHE_SAMPLE_X	NICHE_SAMPLE_Y
```

* * * Each NICHE corresponds to the type of sample i.e. Human Stool, Skin, Soil, Rainforest etc. 
* * * This data is used to calculate the alpha- and beta- prevalence of the gene centroids i.e. prevalence within a specific niche and/or prevalence across different niches
* * * In absence of niche data, only alpha-prevalence is calculated.

* **FASTA_FOLDER**: Location of the fasta files containing gene sequences (amino acids) per sample.
* * See mock fasta files under input/faa_only
* **USEARCH_FOLDER**: Location of the USEARCH PATH = version required 7.0 (?confirm)
* **TYPE**: Defines type of input. input can be gene_table or reads/assemblies data
* **OUTPUT**: Location of output directory

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
``


##METHODOLOGY

1. Cluster centroids according to annotated UniRef90 IDs
2. Cluster centroids annnotated UniRef90_unknown, using 90% sequence homology (USEARCH)
3. Add all genes within the same sample that belong to a specific centroid (UniRef from STEP1 or o/w from STEP2) across samples
4. Take the mean of gene centroid abundance for each gene centroid (where centroid is present; Criteria for presence= >0)
5. Take the prevalence of each gene centroid as #no.samples the gene centroid is present in. (Criteria for presence= >0)
6. Use 10th percentile of the mean abundance and prevalence range as the threshold to categorize importance
7. Output the important genes in the results folder specified by user
8. Dump all the intermediate files in the ``tmp`` folder.

**For NICHE-SPECIFIC ANALYSIS**:
5. Take the prevalence of each gene centroid across samples within EACH NICHE. (Alpha_prevalence)
6. Use the 10th percentile of the mean abundance and alpha_prevalence observed within EACH NICHE as the threshold to categorize importance
* * If the centroid alpha_prevalence within ANY NICHE is >= 10th percentile within that NICHE, then count it as important


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
