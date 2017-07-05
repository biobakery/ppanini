
from __future__ import print_function # PYTHON 2.7+ REQUIREDimport os
import re
import sys
import os
import shutil
import csv
import pdb
import Bio
import numpy
import logging
import argparse
import subprocess
import multiprocessing

from Bio import Seq
import pandas as pd

logger = logging.getLogger(__name__)

# constants
# ---------------------------------------------------------------

c_strat_delim     = "|"
c_taxon_delim     = "."
c_name_delim      = ": "
c_multiname_delim = ";"
c_str_unknown     = "NO_NAME"
c_ungrouped       = "UNGROUPED"
c_unmapped        = "UNMAPPED"
c_unintegrated    = "UNINTEGRATED"
c_many_bytes      = 1e8
c_zip_multiplier  = 10
# the last line in the file with this indicator is the header
GENE_TABLE_COMMENT_LINE="#"

# the extension used for biom files
BIOM_FILE_EXTENSION=".biom"

c_topsort = {
    c_unmapped:0,
    c_ungrouped:1,
    c_unintegrated:2,
    "UniRef50_unknown":3,
    "UniRef90_unknown":4,
}

#from zopy.utils import iter_rows, tprint, warn

#==== Helper function from zopy by Eric Franzosa ========
def warn ( *args ):
    script = "?"
    if sys.argv[0] != "":
        script = os.path.split( sys.argv[0] )[1].upper()
    args = ["WARNING ({}):".format( script )] + list( args )
    print >>sys.stderr, " ".join( map( str, args ) )

def iter_rows( path ):
    """ easy table loading """
    lens = []
    with try_open( path ) as fh:
        for row in reader( fh ):
            lens.append( len( row ) )
            yield row
    if len( set( lens ) ) != 1:
        warn( "rows didn't all have equal lengths:", set( lens ) )

def tprint( *args, **kwargs ):
    """ coerce list of items to strings then print with tabs between """
    print >>kwargs.get( "file", sys.stdout ), "\t".join( map( str, args ) )

def try_open( path, *args ):
    """ open an uncompressed or gzipped file; fail gracefully """
    fh = None
    try:
        if re.search( r".gz$", path ):
            print >>sys.stderr, "Treating", path, "as gzipped file"
            fh = gzip.GzipFile( path, *args )
        else:
            fh = open( path, *args )
    except:
        die( "Problem opening", path )
    return fh   
# ---------------------------------------------------------------
# text manipulation
# ---------------------------------------------------------------

def reader ( file_handle ):
    """ my favorite options for csv reader """
    for aItems in csv.reader( file_handle, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE ):
        yield aItems

def make_directory(output_dir):
    if not os.path.isdir(output_dir):
        try:
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")
    else:
        try:
            print("Removing the old output directory: " + output_dir)
            shutil.rmtree(output_dir)
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")
        
    
    if not os.access(output_dir, os.W_OK):
        sys.exit("CRITICAL ERROR: The output directory is not " + 
            "writeable. This software needs to write files to this directory.\n" +
            "Please select another directory.")
        
    print("Output files will be written to: " + output_dir) 



def read_ppanini_imp_genes_table_dead(filename):
	gene_table = pd.read_csv(filename, sep='\t', index_col=0)
	#print gene_table.columns.values
	#print gene_table['abundance']
	ppanini_table = {'genes': list(gene_table.index), 'abundance_rank': gene_table['abundance'], 'prevalence_rank':gene_table['prevalence'], 'ppanini_score':gene_table['ppanini_score'] }
	return ppanini_table
def read_fasta(fasta_filename):
	'''Reads a fasta_file and returns a fasta dict
	Input: fasta_filename = path_to_fasta_file
	Output: fasta_seq = {sequence_header: sequence, ...}'''
	logger.debug('read_fasta '+fasta_filename)

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
	#print fasta_seq
	return fasta_seq
def parse_table(m8_filename, fasta_filename):
	'''Parse the BLAST results to give gene hits to genomes
	Input: 
	m8_filename = filename of blast results
	fasta_filename = filename of corresponding fasta file
	
	Output: 
	table = {gene: [List of genomes]}'''

	fasta_dict = read_fasta(fasta_filename)

	for seq in fasta_dict:
		fasta_dict[seq] = float(len(fasta_dict[seq]))
	table = {}
	foo = open(m8_filename)
	for line in foo:
		split_i = line.split('\t')
		try:
			threshold = float(split_i[2])*float(split_i[3])/fasta_dict[split_i[0]]
		except:
			raise Exception('Gene '+split_i[0]+' not found in fasta: '+fasta_filename)
		if threshold > 90.0:
			split_sp = split_i[1].split('|')
			sp = [i for i in split_sp if 'g__' in i and '.s__' in i]
			if split_i[0] not in table:
				table[split_i[0]] = sp
			elif sp not in table[split_i[0]]:
				table[split_i[0]] += sp
	return table
def count_genomes(m8_filename):
	table = pd.DataFrame.from_csv(m8_filename, sep='\t', index_col=None, header =None)
	table.columns =['gene', 'genome']
	gene_count_genome = table['gene'].value_counts()
	return gene_count_genome
def gene2genomes(m8_filename):
	'''Read parsed table for {gene: genomes}
	Input: 
	m8_filename = filename of blast results

	Output: 
	table = {gene: [List of genomes]}'''
	
	table ={}
	
	foo = open(m8_filename)

	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		try:
			table[split_i[0]] += [split_i[1]]
		except:
			table[split_i[0]] = [split_i[1]]
	return table
def number_of_unique_genomes(mg_file):
	metagenomic_table  =  pd.DataFrame.from_csv(mg_file, sep='\t', index_col=None, header =None)
	metagenomic_table.columns =['gene', 'genome']
	uniq_genomes = []
	for gene in metagenomic_table:
		for genome in metagenomic_table[gene]:
			if genome not in uniq_genomes:
				uniq_genomes +=[genome]
	no_uniq_genomes = len(uniq_genomes)
	#print 'No. of unique genomes: '+str(no_uniq_genomes)
	#ppanini_output = pd.DataFrame.from_csv(ppanini_output_file, sep='\t', index_col=0, header =0)#read_ppanini_imp_genes_table(ppanini_output_file)
	return no_uniq_genomes 
def read_data(mg_file, ppanini_output_file):
	metagenomic_table  =  pd.DataFrame.from_csv(mg_file, sep='\t', index_col=None, header =None)
	metagenomic_table.columns =['gene', 'genome']
	uniq_genomes = []
	for gene in metagenomic_table:
		for genome in metagenomic_table[gene]:
			if genome not in uniq_genomes:
				uniq_genomes +=[genome]
	no_uniq_genomes = len(uniq_genomes)
	print ('No. of unique genomes: '+str(no_uniq_genomes))
	ppanini_output = pd.DataFrame.from_csv(ppanini_output_file, sep='\t', index_col=0, header =0)#read_ppanini_imp_genes_table(ppanini_output_file)
	return metagenomic_table, ppanini_output, no_uniq_genomes 
def read_abund_prev(filename):
	keys = {'abund':0, 'alpha':0,'beta':0}
	abund = []
	alpha = []
	genes = []
	with open(filename) as foo:
		for line in foo:
			split_line = [re.sub('[\r\t\n]','', i) for i in line.split('\t')]
			if line.startswith('Centroids'):
				for i, val in enumerate(split_line):
					if 'abund' in val:
						keys['abund'] = i
					elif 'alpha' in val:
						keys['alpha'] = i
					elif 'beta' in val:
						keys['beta'] = i
			else:
				# pdb.set_trace()
				genes +=[split_line[0]]
				abund +=[float(split_line[keys['abund']])]
				alpha +=[float(split_line[keys['alpha']])]
	abund_prev = {'genes': genes, 'abundance': abund, 'prevalence': alpha}
	return abund_prev
def pullgenes_fromcontigs(contig_file, gff3_file, fna_file, faa_file):
	'''Pulls genes from contigs using the coordinates from GFF3 file provided

	Input: contig_file: FASTA file containing contigs in FNA format.
		   gff3_file: File containing coordinates of genes in GFF3 format
		   fna_file: filepath for genes written in nucleotides sequences
		   faa_file: filepath for genes written in amino-acid sequences'''
	logger.debug('pullgenes_fromcontigs '+contig_file+' '+gff3_file)

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
						start_minus = -1*stop_x
						stop_minus = -1*(start_x-1)
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
	
	logger.debug('read_gff3 '+filename)

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
	logger.debug('write_fasta to '+filename)

	with open(filename,'w') as foo:

		test = ''
		for i in seqs_dict:
			test = seqs_dict[i]
			break
		format = is_protein(test)
		to_translate =  to_translate and not format
		
		if to_translate: # if not FAA already and to be translated
			for seq in seqs_dict:
				t_seq = Bio.Seq.translate(seqs_dict[seq], to_stop=True)
				len_seq = str(len(t_seq)*3)
				foo.writelines(['>' + seq + '|' + len_seq + '\n'])
				foo.writelines([t_seq+'\n'])
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

	logger.debug('is_protein')

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

	logger.debug('create_folders '+'\t'.join(list_folders))
	
	for fname in list_folders:
		try:
			os.stat(fname)
		except:
		    os.mkdir(fname)
		

def read_dict(gene_annotations_file):
	'''Reads tabulated file into a dictionary

	Input: gene_annotations_file = path_to_output_gene_annotations_table

	Output: dictX = {geneID: annotation}'''

	logger.debug('read_dict '+gene_annotations_file)

	dictX = {}

	with open(gene_annotations_file) as foo:
		for line in foo:
			if not line.startswith('#'):
				split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split('\t')]
				dictX[split_line[0]] = split_line[1]
	return dictX

def read_dict_num(gene_annotations_file):
	'''Reads tabulated file into a dictionary

	Input: gene_annotations_file = path_to_output_gene_annotations_table

	Output: dictX = {geneID: annotation}'''

	logger.debug('read_dict)num '+gene_annotations_file)

	dictX = {}

	with open(gene_annotations_file) as foo:
		for line in foo:
			if not line.startswith('#'):
				split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split('\t')]
				if len(split_line) <2:
					split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split(' ')]
				try:
					dictX[split_line[0]] = float(split_line[1])
				except:
					pdb.set_trace()
	return dictX

def write_dict(dictX, gene_annotations_file):
	'''Writes dictionary of genes and their annotations into text file
	Input: dictX = {geneID: annotation}
		   gene_annotations_file = path_to_output_gene_annotations_table'''

	logger.debug('write_dict to '+gene_annotations_file)

	with open(gene_annotations_file, 'w') as foo:
		foo.writelines('#GENEID\tANNOTATION\n')
		for i in dictX:
			foo.writelines(['\t'.join([i, dictX[i] ])+'\n'])

def is_present(metadata, meta_type):
	'''Returns True if meta_type is present in metadata extracted from mappert_file

	Input: metadata = [metadata strings]; Rows with # as first character in table
		   meta_type = Type of metadata that you are querying e.g. FASTAS, NICHE etc.

	Output: [line, ind]
			line = The corresponding line from metadata, [] if not present
			ind = The index of the line in the metadata sequence, [] if not present'''

	logger.debug('is_present '+meta_type)

	line = []
	ind = []
	for i, val in enumerate(metadata):
		if val.upper().startswith(meta_type):
			line = val
			ind = i
			break
	return [line, ind]
def generate_gene_table(abundance_dict, annotations_dict, niche_flag, mapper, output_table, samples):
	'''Input: abundance_dict: {sample:{gene:abundance,...},...}
			  annotation_dict: {gene: annotation},
			  niche_flag: True if NICHE exists
			  mapper'''
	logger.debug('generate_gene_table '+output_table)
#	samples = abundance_dict.keys()
	# fasta_row = [mapper[i]['FAAS'] for i in samples]

	with open(output_table, 'w') as foo:
		if niche_flag:
			niche_row = [mapper[i]['NICHE'] for i in samples]
			foo.writelines([str.join('\t', ['#NICHE']+niche_row)+'\n']) #header
		# foo.writelines([str.join('\t', ['#FAAS']+fasta_row)+'\n']) #header
		foo.writelines([str.join('\t', ['#SAMPLES']+samples)+'\n']) #header
		
		#for i, sample in enumerate(samples):
		for gene in abundance_dict:
			#abund_x_i = abundance_dict[sample][gene]
			#data_row = numpy.zeros(len(samples))
			#data_row[i] = abund_x_i
			str_data_row = [str(ele) for ele in abundance_dict[gene]]
			if gene in annotations_dict: 
				annot_x_i = annotations_dict[gene]
				if annot_x_i.startswith('UniRef90'):
					umap_i = 'UniRef50_unknown'
					annot_x = annot_x_i + '|' + umap_i
				else:
					annot_x = 'UniRef90_unknown|' + annot_x_i 
			else:
				annot_x = 'UniRef90_unknown|UniRef50_unknown'
			foo.writelines([str.join('\t', [gene+'|'+annot_x]+str_data_row)+'\n'])

def extract_mapping():
	'''Uses centroids from first file and 
	produces mapped centroids to genes through mapper'''
	foo = open(sys.argv[1])
	centroids = [re.sub('[\t\r\n]','',i).strip() for i in foo]
	
	mapper = open(sys.argv[2])
	mapping = {}
	for line in mapper:
		split_line = [re.sub('[\r\n\t]','', i) for i in line.split('\t')]
		try:
			mapping[split_line[0]] += [split_line[1]]
		except:
			mapping[split_line[0]] = [split_line[1]]
	
	for i in centroids:
		if i in mapping:
			for j in mapping[i]:
				print ('\t'.join([i,j]))
		else:
			print (i+'\tNA')
			

# ---------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------
def read_map(map_obj):
	csv_map = csv.reader(open(map_obj), csv.excel_tab)
	uniref_go = {}
	for line in csv_map:
		for uid in line[1:]:
			uniref_go[uid.strip()] = line[0]
	return uniref_go

def size_warn( path ):
    m = 1 if ".gz" not in path else c_zip_multiplier
    if m * os.path.getsize( path ) > c_many_bytes:
		print( "  This is a large file, one moment please...", file=sys.stderr )

def try_zip_open( path, write=None ):
    """ 
    open an uncompressed or gzipped file; fail gracefully 
    """
    fh = None

    # set the open mode
    if write:
        open_mode = "w"
    elif path.endswith(".bz2"):
        open_mode = "r"
    else:
        open_mode = "rt"

    try:
        if path.endswith(".gz"):
            fh = gzip.open( path, open_mode )
        elif path.endswith(".bz2"):
            fh = bz2.BZ2File( path, open_mode )
        else:
            fh = open( path, open_mode )
    except EnvironmentError:
        sys.exit( "Problem opening file: " + path)
    return fh

def read_biom_table( path ):
    """
    return the lines in the biom file
    """

    try:
        import biom
    except ImportError:
        sys.exit("Could not find the biom software."+
            " This software is required since the input file is a biom file.")
        
    try:
        tsv_table = biom.load_table( path ).to_tsv().split("\n")
    except (EnvironmentError, TypeError):
        sys.exit("ERROR: Unable to read biom input file.")
        
    return tsv_table

def gzip_bzip2_biom_open_readlines( path ):
    """
    return the lines in the opened file for tab delimited text, gzip, bzip2 and biom files
    """

    # if the file is biom, convert to text and return lines
    if path.endswith(BIOM_FILE_EXTENSION):
        for line in read_biom_table(path):
            yield line
    else:
        with try_zip_open( path ) as file_handle:
            for line in file_handle:
                if path.endswith(".bz2"):
                    # convert the line to text from binary
                    yield line.decode('utf-8').rstrip()
                else:
                    yield line.rstrip()
def process_gene_table_with_header(gene_table, allow_for_missing_header=None):
    """
    Process through the header portion of the gene table file
    """
    
    # try to open the file
    try:
        lines=gzip_bzip2_biom_open_readlines( gene_table )
    except EnvironmentError:
        sys.exit("Unable to read file: " + gene_table)
            
    # find the headers
    header=""
    first_data_line=""
    for line in lines:
        if line[0] == GENE_TABLE_COMMENT_LINE:
            header = line
        else:
            first_data_line = line
            break
            
    if not header and not allow_for_missing_header:
        sys.exit("File does not have a required header: " + gene_table +
            " . Please add a header which includes the indicator: " +
            GENE_TABLE_COMMENT_LINE)
        
    # provide the header, if one was found
    if header:
        yield header
    
    # provide the first data line
    yield first_data_line
    
    # now provide the remaining lines
    for line in lines:
        yield line

def rev_uniref_mapper ( path_in= '' , path_out ='' , start=0, skip=None, allowed_keys=None, allowed_values=None ):
    """
    Load a file like:
    A 1 2
    B 1
    B 3
    C 1 2 4
    To a nested dict structure:
    {A:{1:1, 2:1}, B:{1:1, 3:1}, C:{1:1, 2:2, 4:1}
    Inner values are not important (set to 1)
    """
    if len(sys.argv)<3:
        sys.exit('Please provide an input and output path')
    path_in = sys.argv[1]
    path_out = sys.argv[2]
    polymap_all = {}
    print( "Loading mapping file from:", path_in, file=sys.stderr )
    size_warn( path_in )
    for line in gzip_bzip2_biom_open_readlines( path_in ):
        row = line.split("\t")
        key = row[start]
        if allowed_keys is None or key in allowed_keys:
            for i, value in enumerate( row ):
                if i != start and (skip is None or i not in skip):
                    if allowed_values is None or value in allowed_values:
                        polymap_all.setdefault( value, {} )[key] = 1 #polymap.setdefault( key, {} )[value] = 1
    
    with gzip.open(path_out+'_dict.txt.gz', 'wt') as csv_file:
        writer = csv.writer(csv_file, delimiter='\t')
        for key, values in polymap_all.items():
           writer.writerow([key, ";".join(values)])                                          
    print("Mappping Uniref90 to annotation is done")

def load_polymap ( path, start=0, skip=None, allowed_keys=None, allowed_values=None ):
    """
    Loads a file to a dictionry 
    INPUT:
    for each row in file format is:
    UniRefID    GOID1;GOID2
    OUTPUT:
    Dictionary : UniRefID is key and value is GOID1;GOID2
    """
    polymap_all = {}
    print( "Loading mapping file from:", path, file=sys.stderr )
    size_warn( path )
    for line in gzip_bzip2_biom_open_readlines( path ):
        row = line.split("\t")
        key = row[start]
        polymap_all[key] = row[1]
    '''with open(path) as csv_file:
        reader = csv.reader(csv_file)
        polymap_all = dict(reader)'''
    return polymap_all
			


if __name__ == '__main__':
	pass
