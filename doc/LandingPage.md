#PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration#
PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g. prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools.
----

 * Download the PPANINI software ([ppanini.tar.gz](https://bitbucket.org/biobakery/ppanini/downloads/biobakery-ppanini-0.6.0.tar.gz)) then follow the [steps to install and run](https://bitbucket.org/biobakery/ppanini/).
 * For additional information, please see the [PPANINI User Manual](https://bitbucket.org/biobakery/ppanini/). [PPANINI wiki page](https://bitbucket.org/biobakery/biobakery/wiki/PPANINI) provides a tutorial for a quick start. 
 * Please direct questions to the [PPANINI google group](https://groups.google.com/forum/#!forum/ppanini-users) (subscribe to receive PPANINI news).
 * If you use the PPANINI software, please cite our manuscript: **Systematic approach to prioritization of 'important' microbial genes**, Gholamali Rahnavard, Afrah Shafquat, Bahar Sayoldin, Eric A. Franzosa, Curtis Huttenhower (under preparation).
  	

----
#PPANINI Workflow#
![ppanini_workflow.png](https://bitbucket.org/repo/rnag7z/images/5336251-ppanini_workflow.png)
----
#Demo#
##Communities from Human Microbiome Project (HMP) to start with##

[A genes abundances table](https://www.dropbox.com/s/drxvgs42iyvk5k0/stool_gene_centroids_table.txt?dl=0) for 93 stool samples with a
[UCLUST file](https://www.dropbox.com/s/b8ufu3ryiyuo3ax/stool_gene_clusters.uc?dl=0) containing centroids of genes. The UCLUST file is used to collapse unannotated
genes into artifical clusters. This step could be bypassed.
```
$ppanini -i stool_gene_centroids_table.txt --uc stool_gene_clusters.uc -o OUTPUT 

for bypassing clustering unannotated genes:
$ppanini -i stool_gene_centroids_table.txt --bypass-clustering -o OUTPUT 
```
The output is a table of prioritized important genes with their prevelance, abundances, and ppanini score.

[An genes abundances table](https://www.dropbox.com/s/lnpef7hixuimm62/AN_gene_table.txt?dl=0) for 70 Anterior nares samples with a 
[gene catalog fasta file](https://www.dropbox.com/s/2mohfte3lkplqsy/AN_centroids_for_clustering.fasta?dl=0) which is used to cluster unannoted genes.
```
$ppanini -i AN_gene_table.txt --gene-catalog AN_centroids_for_clustering.fasta -o OUTPUT 

for bypassing clustering unannotated genes:
$ppanini -i AN_gene_table.txt --bypass-clustering -o OUTPUT 
```

##Simple Demo##
###input###
Download [Gene abaundances table](https://www.dropbox.com/s/utrjt28sxn16glu/genetable.txt?dl=0) and a [FASTA file](https://www.dropbox.com/s/2bgyid79rf97lg0/samples.fasta?dl=0) for clustering unannotated genes

###Running Command###
```
$ppanini -i genetable.txt --gene-catalog samples.fasta -o OUTPUT --vsearch /path/to/vsearch

or clustring step could be bypassed by:
$ppanini -i genetable.txt --bypass-clustering -o OUTPUT 
```

###Output###
[List](https://www.dropbox.com/s/c1q7zw90uuekx2k/genetable_imp_centroid_prev_abund.txt?dl=0) of importnat genes(centroids) with prevelence, abundance, and ppanini score is the output. 

----
# PPANINI Evaluation#
An evaluation of PPANINI in terms of sensitivity and specificity has been used. Stool samples from HMP datasets are used to find prioritized important genes and a union of two essential genes sets is used as golden standard for the ROC plot. 
The result shows PPANINI prioritized important genes very well (AUC= .87) even we used a specific gene set as golden standard. 

![roc_plot_ppanini_union_of_essential_gene_datasets.png](https://bitbucket.org/repo/49y6o9/images/3568610095-roc_plot_ppanini_union_of_essential_gene_datasets.png)