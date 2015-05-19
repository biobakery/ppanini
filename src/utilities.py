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

    for line in fasta_file:
        if not line.startswith('#'):
            if line.startswith('>'):
                name = re.sub('>','', line.split(' ')[0].strip())
            else:
                if name not in fasta_seq:
                    fasta_seq[name] =  re.sub('[\r\t\n]','', line)
                else:
                    fasta_seq[name] +=  re.sub('[\r\t\n]','', line)
    return fasta_seq

def pullgenes_fromcontigs(contig_file, gff3_file, fna_file, faa_file):
    '''Pulls genes from contigs using the coordinates from GFF3 file provided

    Input: contig_file: FASTA file containing contigs in FNA format.
           gff3_file: File containing coordinates of genes in GFF3 format
           fna_file: filepath for genes written in nucleotides sequences
           faa_file: filepath for genes written in amino-acid sequences'''

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

    write_fasta(genes_fasta, fna_file, False) #FNA
    write_fasta(genes_fasta, faa_file, True) #FAA


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
        for line in foo:
            if re.match('(\w+)\t(\w+)\tgene\t(\w+)', line):
                split_line = [re.sub('[\r\t\n]', '', i).strip() for i in line.split('\t')]
                gid = split_line[-1].split('=')[-1]
                gene_contig_mapper[gid] = split_line[0]
                gene_start_stop[gid] = [int(split_line[3]), int(split_line[4]), split_line[6]]

                if split_line[0] in contig_gene_mapper:
                    contig_gene_mapper[split_line[0]] += [gid]
                else:
                    contig_gene_mapper[split_line[0]] = [gid]

    return [gene_contig_mapper, gene_start_stop, contig_gene_mapper]


def write_fasta(seqs_dict, filename, to_translate):
    '''Writes dictionary of fasta sequences into text file
    Input: seqs_dict = {geneID: sequence}
           filename = path_to_output_genes_fastafile'''

    with open(filename,'w') as foo:

        test = ''
        for i in seqs_dict:
            test = seqs_dict[i]
            break
        format = is_protein(test)
        to_translate =  to_translate and not format
        
        if to_translate: # if not FAA already and to be translated
            for seq in seqs_dict:
                len_seq = str(len(seqs_dict[seq]))
                foo.writelines(['>' + seq + '|' + len_seq + '\n'])
                foo.writelines([Bio.Seq.translate(seqs_dict[seq], to_stop=True)+'\n'])
        else:
            ind = 1
            if format:
                ind = 3 #amino acids * 3 nucleotides
            for seq in seqs_dict:
                len_seq = str(len(seqs_dict[seq]*ind))
                foo.writelines(['>' + seq + '|' + len_seq + '\n'])
                foo.writelines([seqs_dict[seq] + '\n'])

def is_protein(sequence):
    '''Returns True if the sequence is protein.
    Input: (str) format sequence
    Output: (bool) True if amino acids; False if nucleotides sequence'''

    format = False  #'FNA'
    try:
        Bio.Seq.translate(sequence)
    except:
        format = True #'FAA'
    return format

def create_folders(list_folders):
    '''Creates the list of folders if they dont exist already
    Input: list_folders = [List of folders to be created]

    Output: Folders created'''

    for fname in list_folders:
        try:
            os.mkdir(fname)
        except:
            pass

def read_dict(gene_annotations_file):
    '''Reads tabulated file into a dictionary

    Input: gene_annotations_file = path_to_output_gene_annotations_table

    Output: dictX = {geneID: annotation}'''

    dictX = {}

    with open(gene_annotations_file) as foo:
        for line in foo:
            if not line.startswith('#'):
                split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split('\t')]
                dictX[split_line[0]] = split_line[1]
    return dictX

def write_dict(dictX, gene_annotations_file):
    '''Writes dictionary of genes and their annotations into text file
    Input: dictX = {geneID: annotation}
           gene_annotations_file = path_to_output_gene_annotations_table'''

    with open(gene_annotations_file, 'w') as foo:
        foo.writelines('#GENEID\tANNOTATION')
        for i in dictX:
            foo.writelines(['\t'.join([i, dictX[i]])+'\n'])

if __name__ == '__main__':
    pass