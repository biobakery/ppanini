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
  -u [UCLUST_FOLDER], --uclust_folder [UCLUST_FOLDER]
                        Path for UCLUST program
  -t TYPE, --type TYPE  Type of analysis Choose: [gene_table,
                        reads_assemblies]
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Folder containing results

```


## **Inputs**

### **INPUT_TABLE**

See the mock gene table for an example. naan/input/mock_gene_table.tsv

**Metadata**

*Required*

```
#!text
#FASTA	N/A	location_fasta_SAMPLE_X	location_fasta_SAMPLE_Y
#GENEID	ANNOTATION	SAMPLE_X	SAMPLE_X
```

*Optional*


```
#!text
#NICHE	N/A	NICHE_SAMPLE_X	NICHE_SAMPLE_Y
```

* Each NICHE corresponds to the type of sample i.e. Human Stool, Skin, Soil, Rainforest etc. 
* This data is used to calculate the alpha- and beta- prevalence of the gene centroids i.e. prevalence within a specific niche and/or prevalence across different niches
* In absence of niche data, only alpha-prevalence is calculated.



