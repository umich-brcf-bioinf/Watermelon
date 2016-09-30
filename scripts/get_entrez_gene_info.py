#!/usr/bin/python

import sys
import os
import csv
import datetime
import time
import argparse
from collections import defaultdict

def time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def log(message):
    print('{}|{}'.format(time_stamp(), message), file=sys.stderr)

def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description='Script adds NCBI gene_info annotations to cuffdiff output')

    parser.add_argument(
        '-g', '--geneinfo', type=str, help='path to input NCBI gene_info file', required=True)
    parser.add_argument(
        '-e', '--diffexp', type=str, help='path to input differential_expression_file', required=True)
    parser.add_argument(
        '-t', '--taxid', type=str, help='organism taxonomy ID', required=True)
    parser.add_argument(
        '-o', '--outdir', type=str, help='output directory name', required=True)
    # Array for all arguments passed to script
    args = parser.parse_args()
    return args #gene_info, diffexp, tax_id, output_dir
    
args = parse_command_line_args()
gene_info = args.geneinfo
gene_expr = args.diffexp
tax_id = args.taxid
outdir = args.outdir
print(gene_info + ', ' + gene_expr + ', ' + tax_id+ ', ' + outdir)

outfile_tag = os.path.basename(gene_expr.replace('.txt', '_annot.txt')) # instead of writing to a file write to stdout
outfile_name = os.path.join(outdir, outfile_tag)
print(outfile_name)

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
