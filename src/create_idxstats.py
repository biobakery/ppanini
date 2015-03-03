import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing

def generate_abundance_viabwt2(assembly_x_withpath, reads_x, sample):
  '''Calculate Genes/Contigs abundance from Contig_assemblies and reads
  Input: assembly_x_withpath= path_to_assemblies_file
             reads_x = path_to_reads_file To go from ASSEMBLIES, READS to IDXSTATS
             sample = sample_name

  Output: generate_abundance_viabam(assembly_x_bam, sample)

  where assembly_x_bam = path_to_bam_file, 
            sample = sample_name'''
  #assembly_x_withpath = mapper[sample]['#ASSEMBLIES'] #Full path inlcuded to Assembly X
  assembly_x = assembly_x_withpath.rpartition('/')[-1] #Name of Assembly X
  #reads_x = mapper[sample]['#READS'] #Path included to Reads of X
  assembly_x_index = 'tmp/' + assembly_x + '_index' # 
  assembly_x_bam = 'tmp/' + assembly_x + '.bam'

  os.system('bowtie2-build --quiet ' + assembly_x_withpath + ' ' + assembly_x_index)
  os.system('tar -xOvf ' + reads_x + ' | \
                     bowtie2 -x ' + assembly_x_index + ' -U - --no-unal --very-sensitive | \
                     samtools view -bS - > ' + assembly_x_bam)
  return generate_abundance_viabam(assembly_x_bam, sample)


def generate_abundance_viasam(assembly_x_sam_withpath, sample):
  '''Calculate Genes/Contigs abundance from SAM file
  Input: assembly_x_sam_withpath= path_to_sam_file
             sample = sample_name

  Output: generate_abundance_viabam(assembly_x_bam, sample)

  where assembly_x_bam = path_to_bam_file, 
            sample = sample_name'''

  assembly_x_sam = assembly_x_sam_withpath.rpartition('/')[-1] #Name of Assembly X
  assembly_x_bam = 'tmp/' + assembly_x_sam + '.bam'

  os.system('samtools view -bS ' + assembly_x_sam_withpath + ' > ' + assembly_x_bam)

  return generate_abundance_bam(assembly_x_bam, sample)

def generate_abundance_viabam(assembly_x_bam_withpath, sample):
  '''Calculate Genes/Contigs abundance from BAM file
  Input: assembly_x_bam_withpath= path_to_bam_file
           sample = sample_name

  Output: [assembly_x_stats, sample]

  where assembly_x_stats = path_to_samtools_abundance_file, 
          sample = sample_name'''

  try:
    os.mkdir('tmp/idxstats')
  except:
    pass

  assembly_x_bam = assembly_x_bam_withpath.rpartition('/')[-1]
  assembly_x_bam_presort = 'tmp/' + assembly_x_bam + '.sorted'
  assembly_x_bam_sorted = 'tmp/' + assembly_x_bam + '.sorted.bam'
  #assembly_x_stats = 'tmp/' + assembly_x_bam + '.txt'
  assembly_x_stats = 'tmp/idxstats/' + sample + '.txt'

  os.system('samtools sort ' + assembly_x_bam_withpath + ' ' + assembly_x_bam_presort)
  os.system('samtools index ' + assembly_x_bam_sorted)
  os.system('samtools idxstats ' + assembly_x_bam_sorted + ' > ' + assembly_x_stats)

  return [assembly_x_stats, sample]
  
def read_abundance_tables(mapper):
  ''' Reads the abundance tables generated from get_abundance and
    returns a dictionary cataloging all samples and their component abundances (genes or contigs)

  Input: mapper = {sample:{'READS': path_to_reads_file, 
                 'CONTIG_ASSEMBLIES': path_to_assemblies, 
                 'FASTAS': path_to_fasta_file,
                 'SAMS': path_to_sam_file,
                 'BAMS': path_to_bam_file,
                 'GFF3S': path_to_gff3_file,
                 'NICHE': niche_metadata,
                 'abundance_file': path_to_abundance_file_from_get_abundance}, ...}

  Output: abundance_dict= {sample:{ID: abundance}, ...}
  * ID can be GENE_ID or CONTIG_ID depending on workflow '''

  abundance_dict = {}
  for sample in mapper:
    abundance_dict[sample] = {}
    foo = open(mapper[sample]['abundance_file'])
    for line in foo.readlines():
      if not (line.startswith('#') and line.startswith('*')):
        #FILE_FORMAT: Seq_Name<\t>Seq_Length<\t>Mapped_Reads<\t>Unmapped_Reads
        split_line = [re.sub('[\t\r\n]', '', i) for i in line.split('\t')]
        #abundance_dict[sample][contig_name] = [number_mapped_reads, sequence_length]
        abundance_dict[sample][split_line[0].strip()] = float(split_line[2])/(float(split_line[1])-100.0) #avg. read length?
  
  return abundance_dict

def update_abundance_dict(abundance_dict, gene_contig_mapper):
  '''Update abundance_dict from {sample:{contig:abundance}} to {sample:{gene:abundance}} via gene to contig mapping

  Input: abundance_dict= {sample:{contig: abundance}, ...}
       gene_contig_mapper={sample:{gene:contig}, ...}
  
  Output: abundance_dict_mod= {sample:{gene:abundance}, ...}'''
  
  abundance_dict_mod = {}

  for sample in gene_contig_mapper:
    abundance_dict_mod[sample] = {}
    for gene in gene_contig_mapper[sample]:
      abundance_dict_mod[sample][gene] = abundance_dict[sample][gene_contig_mapper[sample][gene]]

  return abundance_dict_mod
  
if __name__ == '__main__':
  reads = sys.argv[2]
  assembly = sys.argv[1]
  sample = sys.argv[3]
  generate_abundance_viabwt2(assembly, reads, sample)