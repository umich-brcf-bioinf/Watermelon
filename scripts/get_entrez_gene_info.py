#!/usr/bin/python

import sys
import os
import csv
import datetime
import time
from collections import defaultdict

def time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def log(message):
    print('{}|{}'.format(time_stamp(), message), file=sys.stderr)

#?? add code to validate the command line argument; error handling
gene_info = sys.argv[1]         #'NCBI_annotation_2016_06_09/NCBI_gene_info_2016_07_12'
gene_expr = sys.argv[2]         #'gene_diffexp.txt' #C.Plus_v_Six.Plus_gene.foldchange.1.5.txt'
tax_id = sys.argv[3]            #'9606'
outfile_tag = os.path.basename(gene_expr.replace('.txt', '_annot.txt')) # instead of writing to a file write to stdout
outfile_name = os.path.abspath(outfile_tag)

log('reading gene info')
gene_details = defaultdict(list)
with open(gene_info,'r') as geneinfo_file:
    next(geneinfo_file) # skip header
    reader=csv.reader(geneinfo_file,delimiter='\t')
    for row in reader:
        (geneinfo_tax_id, geneinfo_gene_id, geneinfo_gene_symbol, geneinfo_gene_description) = (row[0], row[1], row[2], row[8])
        if geneinfo_tax_id == tax_id:
            gene_details[geneinfo_gene_symbol].append(geneinfo_gene_id)
            gene_details[geneinfo_gene_symbol].append(geneinfo_gene_description)

with open(gene_expr, 'r') as file_to_annotate, open(outfile_name, 'w') as annotated_file: 
    reader=csv.reader(file_to_annotate,delimiter='\t')

    for row in reader:
        row = [col.strip() for col in row]
        if row[0] == '#test_id':
            row[1] = 'gene_symbol'
            row[2] = 'gene_id'
            row.insert(3, 'gene_desc')
            print('\t'.join(row), file=annotated_file)
        elif row[0].startswith('#'):
            print('\t'.join(row), file=annotated_file)
        else:
            left = row[0:2]
            right = row[3:]
            gene_name = row[1]
            gene_id_symbol = gene_details.get(gene_name, ['.','.']) # If key is not available then return default value ('.') 
            outline = "\t".join(left + gene_id_symbol + right)  
            print(outline, file=annotated_file)

log('done')
