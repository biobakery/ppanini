import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing
import Bio
from Bio import Seq

def read_fasta(fasta_filename):
    '''Reads a fasta_file and returns a fasta dict
    Input: fasta_filename = path_to_fasta_file
    Output: fasta_seq = {sequence_header: sequence, ...}'''

    fasta_file = open(fasta_filename)
    fasta_seq = {}

    name = ''

    for line in fasta_file.readlines():
        if not line.startswith('#'):
            if line.startswith('>'):
                name = line.split(' ')[0][1:].strip()
            else:
                if name not in fasta_seq:
                    fasta_seq[name] =  re.sub('[\r\t\n]','', line)
                else:
                    fasta_seq[name] +=  re.sub('[\r\t\n]','', line)
    return fasta_seq


def pullgenes_fromcontigs(contig_file,gff3_file, output_file):
    '''Return dictionary of annotations for genes from each sample's contig_assembly

    Input: mapper = {sample:{'READS': path_to_reads_file, 
                             'CONTIG_ASSEMBLIES': path_to_assemblies, 
                             'FASTAS': path_to_fasta_file,
                             'SAMS': path_to_sam_file,
                             'BAMS': path_to_bam_file,
                             'GFF3S': path_to_gff3_file,
                             'NICHE': niche_metadata}, ...}
               all_paths = {'uniref90': path_to_uniref90_index, 
                                            'uniref50': path_to_uniref50_index, 
                                            'umap90_50': path_to_uniref90_uniref50_mapping}
               nprocesses = int, Number of processes

    Output: [annotations_dict, gene_contig_mapper]
    annotations_dict = {sample :{gene: UniRef annotation}}
    gene_contig_mapper = {sample: {gene: contig}}'''

    gene_contig_mapper = {}
    contig_gene_mapper = {}
    gene_start_stop = {}

    [gene_contig_mapper, gene_start_stop, contig_gene_mapper] = read_gff3(gff3_file)

    genes_fasta = {}
    contigs_fasta_dict = read_fasta(contig_file)
    for contig in contigs_fasta_dict:
        if contig in contig_gene_mapper:
            for gene in contig_gene_mapper[contig]:
                [start_x, stop_x, strand] = gene_start_stop[gene]
                try:
                    if strand == '+':
                        genes_fasta[gene] = contigs_fasta_dict[contig][start_x-1:stop_x+1]
                    else:
                        contig_len = len(contigs_fasta_dict[contig])
                        start_minus = contig_len - stop_x
                        stop_minus = contig_len - start_x+1
                        genes_fasta[gene] = Bio.Seq.reverse_complement(contigs_fasta_dict[contig])[start_minus:stop_minus]
                except:
                    raise Exception('Circular DNA Detected')
                    new_stop_x = stop_x - len(contigs_fasta_dict[contig])+1
                    genes_fasta[gene] = contigs_fasta_dict[contig] + contigs_fasta_dict[contig][:new_stop_x]

    write_fasta(genes_fasta, output_file)
    return gene_contig_mapper


def read_gff3(filename):
    '''Reads GFF3 files and returns the relevant information
    Input: filename = path_to_gff3_file
    Output: [gene_contig_mapper, gene_start_stop, contig_gene_mapper]
    gene_contig_mapper = {gene: contig, ...}
    gene_start_stop = {gene: [start, stop, strand], ...}
    contig_gene_mapper = {contig: [List of genes], ...}'''

    gene_contig_mapper = {}
    contig_gene_mapper = {}
    gene_start_stop = {}

    with open(filename,'r') as foo:
        foo_lines = foo.readlines()
        for line in foo_lines:
            if re.match('(\w+)\t(\w+)\tgene\t(\w+)', line):
                split_line = [re.sub('[\r\t\n]','', i).strip() for i in line.split('\t')]
                gid = split_line[-1].split('=')[-1]
                gene_contig_mapper[gid] = split_line[0]
                gene_start_stop[gid] = [int(split_line[3]), int(split_line[4]), split_line[6]]

                if split_line[0] in contig_gene_mapper:
                    contig_gene_mapper[split_line[0]] += [gid]
                else:
                    contig_gene_mapper[split_line[0]] = [gid]

    return [gene_contig_mapper, gene_start_stop, contig_gene_mapper]


def write_fasta(seqs_dict, filename):
    '''Writes dictionary of fasta sequences into text file
    Input: seqs_dict = {geneID: sequence}
           filename = path_to_output_genes_fastafile
    '''
    with open(filename,'w') as foo:
        test = if_protein(seqs_dict.values()[0])
        

        if not format:
            for seq in seqs_dict:
                foo.writelines(['>'+seq+'\n'])
                foo.writelines([Bio.Seq.translate(seqs_dict[seq], to_stop=True)+'\n'])
        else:
            for seq in seqs_dict:
                foo.writelines(['>'+seq+'\n'])
                foo.writelines([seqs_dict[seq]+'\n'])

def if_protein(sequence):
    format = False  #'FNA'
    try:
        Bio.Seq.translate(sequence)
    except:
        format = True #'FAA'
    return format

if __name__ == '__main__':
    contig_file = sys.argv[1]
    gff3_file = sys.argv[2]
    output_file = sys.argv[3]
    pullgenes_fromcontigs(contig_file, gff3_file, output_file)