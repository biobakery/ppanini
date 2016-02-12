#PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration#
PPANINI provides a computational pipeline to prioritize microbial genes based on their metagenomic properties (e.g. prevalence and abundance). The resulting prioritized list of gene candidates can then be analyzed further using our visualization tools.
----

 * Download the PPANINI software ([ppanini.tar.gz](https://bitbucket.org/biobakery/ppanini/downloads/biobakery-ppanini-0.6.0.tar.gz)) then follow the [steps to install and run](#markdown-header-getting-started-with-ppanini).

 * For additional information, please see the [PPANINI User Manual](https://bitbucket.org/biobakery/ppanini/).

 * Please direct questions to the [PPANINI google group](https://groups.google.com/forum/#!forum/ppanini-users) (subscribe to receive PPANINI news).

 * If you use the PPANINI software, please cite our manuscript: **Systematic approach to prioritization of 'important' microbial genes**, Gholamali Rahnavard, Afrah Shafquat, Bahar Sayoldin, Eric A. Franzosa, Curtis Huttenhower (under preparation).
  	

----
#PPANINI Workflow#
<center>![ppanini_workflow.png](https://bitbucket.org/repo/rnag7z/images/994033213-ppanini_workflow.png)</center>

----
#Communities from Human Microbiome Project (HMP) to start with#
[An abundances table](https://www.dropbox.com/s/drxvgs42iyvk5k0/stool_gene_centroids_table.txt?dl=0) for genes from stool bodysite samples with a [UCLUST file](https://www.dropbox.com/s/b8ufu3ryiyuo3ax/stool_gene_clusters.uc?dl=0) containing centroids of genes. 

[An abundances table](https://www.dropbox.com/s/lnpef7hixuimm62/AN_gene_table.txt?dl=0) for genes from Anterior nares bodysite samples  with a [gene catalog fasta file](https://www.dropbox.com/s/2mohfte3lkplqsy/AN_centroids_for_clustering.fasta?dl=0)

*Simple Demo*
input 
Download [Gene abaundances table](https://www.dropbox.com/s/utrjt28sxn16glu/genetable.txt?dl=0)
Download [FASTA file](https://www.dropbox.com/s/2bgyid79rf97lg0/samples.fasta?dl=0) for clustering unannotated genes
```
ppanini -i Mock_FASTA/mockgenetable.txt --gene-catalog Mock_FASTA/all.fasta -o OUTPUT --vsearch ./vsearch/vsearch-1.1.3-osx-x86_64
```
----
# PPANINI Evaluation#
An evaluation of PPANINI using stool sample against union of two essentail gense sets using ROC plot.

<center>![roc_plot_ppanini_union_of_essential_gene_datasets.png](https://bitbucket.org/repo/49y6o9/images/3568610095-roc_plot_ppanini_union_of_essential_gene_datasets.png)</center>