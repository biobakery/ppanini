#!/usr/bin/env python

"""
Join a set of gene tables into a single table

This module will join gene tables output by ppanini. 

Dependencies: Biom (only required if running with .biom files)

To Run: 
$ ./join_tables.py -i <input_dir> -o <gene_table.{tsv,biom}>

"""

import argparse
import sys
import tempfile
import os
import shutil
import re
from ..utilities import gzip_bzip2_biom_open_readlines, process_gene_table_with_header, rev_load_polymap, load_polymap_dic
from .. import config

GENE_TABLE_DELIMITER="\t"
# ---------------------------------------------------------------
# utilities used by the split and join tables scripts
# ---------------------------------------------------------------

# indicator of a comment line or the header
# the last line in the file with this indicator is the header
GENE_TABLE_COMMENT_LINE="#"

# the extension used for biom files
BIOM_FILE_EXTENSION=".biom"

        
        
def join_gene_tables(gene_tables,output,verbose=None, mapper= None, scale = None):
    """
    Join the gene tables to a single gene table
    """
    
    gene_table_data={}
    start_column_id=""
    samples=[]
    file_basenames=[]
    index=0
    find_count = 0
    miss_count = 0
    for gene_table in gene_tables:
        
        if verbose:
            print("Reading file: " + gene_table)
        
        lines=process_gene_table_with_header(gene_table, allow_for_missing_header=True)
        header=next(lines)
        
        # get the basename of the file
        file_basename='.'.join(os.path.basename(gene_table).split('.')[:-1])
        file_basenames.append(file_basename)
        
        if False:
            header_info=header.split(GENE_TABLE_DELIMITER)
            if not start_column_id:
                start_column_id=header_info[0]
            # allow for multiple samples
            sample_names=header_info[1:]
        else:
            # if there is no header in the file then use the file name as the sample name
            sample_names=[file_basename]
        
        for line in lines:
            data=line.split(GENE_TABLE_DELIMITER)
            # Skip the header line or incomplete lines
            if len(data) < 7 or data[6].startswith('/'):
                continue
            try:
                # Normalize counts
                if scale == 'rpk':
                    if int(data[5]) >0:
                        data_points = [str(int(data[6]) * 1000.0/int(data[5]))]#hits * 1000/len 
                    else:
                        # use row counts if the count is zero
                        data_points = [data[6]]
                elif scale == 'count':
                    # no scale, use the raw counts
                    data_points = [data[6]]
                else:system.exit("scale is not valid!")
                #print (data[1]+'_'+data[0])
                
                # bind sample name and gene id to use as gene name
                gene = data[1]+'_'+data[0].split("_")[1]
                
                # add the gene abundance to its cluster and use its cluster name
                if mapper and gene in mapper:
                    gene = mapper[gene].keys()[0]
                    #print ("mapping cluster found for:", gene )
                    find_count += 1

                else:
                    print ("No mapping cluster found for:", gene )
                    miss_count += 1
                    continue
            except IndexError:
                gene=""

            if gene:
                current_data=gene_table_data.get(gene,"")
                fill = index - current_data.count(GENE_TABLE_DELIMITER)
                if fill > 0:
                    # fill in zeros for samples without data then add data point
                    gene_table_data[gene]=current_data + GENE_TABLE_DELIMITER.join(["0"]*fill) + GENE_TABLE_DELIMITER + GENE_TABLE_DELIMITER.join(data_points) + GENE_TABLE_DELIMITER
                elif fill < 0:
                    # add data point to other data point from the same sample
                    current_data_points=current_data.split(GENE_TABLE_DELIMITER)
                    for i, point in enumerate(data_points):
                        store_index=len(data_points)*-1-1+i
                        current_data_points[store_index]=str(float(current_data_points[store_index])+float(point))
                    gene_table_data[gene] = GENE_TABLE_DELIMITER.join(current_data_points)
                else:
                    # add data point to end of list
                    gene_table_data[gene] = current_data + GENE_TABLE_DELIMITER.join(data_points) + GENE_TABLE_DELIMITER
        samples+=sample_names
        index+=len(sample_names)
    # if all of the header names for the files are the same
    # then use the file names as headers

    if samples.count(samples[0]) == len(samples):
        samples=file_basenames
                
    # write the joined gene table
    if not start_column_id:
        start_column_id="# header "
    sample_header=[start_column_id]+samples
    total_gene_tables=len(samples)
    sorted_gene_list=sorted(list(gene_table_data))
    try:
        file_handle=open(output,"w")
        file_handle.write(GENE_TABLE_DELIMITER.join(sample_header)+"\n")
    except EnvironmentError:
        sys.exit("Unable to write file: " + output)  
        
    for gene in sorted_gene_list:
        # extend gene data for any gene that is not included in all samples
        current_data=gene_table_data[gene]
        fill = total_gene_tables - current_data.count(GENE_TABLE_DELIMITER)
        if fill:
            current_data=current_data + GENE_TABLE_DELIMITER.join(["0"]*fill) + GENE_TABLE_DELIMITER
        file_handle.write(gene+GENE_TABLE_DELIMITER+current_data.rstrip(GENE_TABLE_DELIMITER)+"\n")
    
    print "Matched genes to clusters : ", find_count, " Unmatched genes: ", miss_count
    file_handle.close()

def get_args():
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Join gene, pathway, or taxonomy tables\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i","--input",
        help="the directory of tables\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the table to write\n",
        required=True)
    parser.add_argument(
        "--file_name",
        help="only join tables with this string included in the file name")
    parser.add_argument(
        "-s","--search-subdirectories", 
        help="search sub-directories of input folder for files\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--mapping-uniref",
        dest= 'mapping_uniref', 
        help="Mapping file: gene to uniref90 file\n", 
        default='')
    parser.add_argument(
        "--mapping-cluster",
        dest= 'mapping_cluster', 
        help="Mapping file: cluster to genes file\n", 
        default='')
    parser.add_argument(
        "-r","--resume", 
        help="bypass commands if the output files exist\n", 
        action="store_true",
        default=config.resume)
    parser.add_argument(
        "--scale",
        dest= 'scale', 
        help="scale the abundance table\n",
        choices=["rpk","count"], 
        default='rpk')
    

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=get_args()
    config.resume = args.resume
    # check for format of the gene tables
    input_dir=os.path.abspath(args.input)
    
    # check the directory exists
    if not os.path.isdir(input_dir):
        sys.exit("The input directory provided can not be found." + 
            "  Please enter a new directory.")
    polymap = None
    if args.mapping_uniref != '':
        print ("Loading mapping uniref-gene file ...")
        polymap =  rev_load_polymap ( path_in= args.mapping_uniref , path_out ='' , 
                                     start=0, skip=None, allowed_keys=None, allowed_values=None, write_output = False, sep = '\t' )
        #print("UniRef Mapper: ",polymap)
    if args.mapping_cluster != '':
        print ("Loading mapping cluster-genes file ...")
        temp_map = rev_load_polymap ( path_in= args.mapping_cluster , path_out ='' , 
                                     start=0, skip=None, allowed_keys=None, allowed_values=None, write_output = False, sep = ';' )
        #print("Mapper :",temp_map)
        polymap.update(temp_map)

    gene_tables=[]
    file_list=os.listdir(input_dir)
    # add in files in subdirectories, if set
    if args.search_subdirectories:
        for possible_folder in os.listdir(input_dir):
            if os.path.isdir(os.path.join(input_dir,possible_folder)):
                try:
                    file_list+=[os.path.join(possible_folder, file) for file in os.listdir(os.path.join(input_dir,possible_folder))]
                except EnvironmentError:
                    pass
    
    # filter out files which do not meet the name requirement if set
    biom_flag=False
    reduced_file_list=[]
    if args.file_name:
        for file in file_list:
            if re.search(args.file_name,file):
                reduced_file_list.append(file)
                if file.endswith(BIOM_FILE_EXTENSION):
                    biom_flag=True
    else: 
        for file in file_list:
            # ignore dot files, like ".DS_Store" on Apple OS X
            if file[0] == ".":
                if args.verbose:
                    print("Not including file in input folder: " + file)
            else:
                reduced_file_list.append(file)
                if file.endswith(BIOM_FILE_EXTENSION):
                    biom_flag=True
    file_list=reduced_file_list
            
    # Check for the biom software if running with a biom input file
    if biom_flag:
        if not util.find_exe_in_path("biom"):
            sys.exit("Could not find the location of the biom software."+
            " This software is required since the input file is a biom file.")       
    
    args.output=os.path.abspath(args.output)
    output_dir=os.path.dirname(args.output)
    
    # Create a temp folder for the biom conversions
    if biom_flag:
        temp_dir=tempfile.mkdtemp( 
            prefix='ppanini_split_gene_tables_', dir=output_dir)
        if args.verbose:
            print("Temp folder created: " + temp_dir)
    
    for file in file_list:
        if file.endswith(BIOM_FILE_EXTENSION):
            # create a new temp file
            file_out, new_file=tempfile.mkstemp(dir=temp_dir)
            os.close(file_out)
        
            # convert biom file to tsv
            if args.verbose:
                print("Processing file: " + os.path.join(input_dir,file))
            util.biom_to_tsv(os.path.join(input_dir,file),new_file)
            gene_tables.append(new_file)
        elif file.endswith('.txt') or file.endswith('.tsv') or file.endswith('.csv'):
            gene_tables.append(os.path.join(input_dir,file))
            
    # sort the gene tables so they are in the same order on all platforms
    gene_tables.sort()
        
    # split the gene table
    if gene_tables:
        if args.verbose:
            print("Joining gene table")
            
        if biom_flag:
            # create a new temp file
            file_out, new_file=tempfile.mkstemp(dir=temp_dir)
            os.close(file_out)
            join_gene_tables(gene_tables,new_file)
            util.tsv_to_biom(new_file, args.output)
        else:
            join_gene_tables(gene_tables,args.output,verbose=args.verbose, mapper = polymap, scale = args.scale)
                
        # deleting temp folder with all files
        if biom_flag:
            if args.verbose:
                print("Deleting temp files in temp folder: " + temp_dir)
            shutil.rmtree(temp_dir)
        if polymap:
            print("Gene families table created: " + args.output)
        else:
            print("Gene table created: " + args.output)
    else:
        print("Zero gene tables were found to join.")

if __name__ == "__main__":
    main()
