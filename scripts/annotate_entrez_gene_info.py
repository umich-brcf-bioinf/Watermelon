#!/usr/bin/python

# adds annotations to output of: scripts/flag_diffex.py

from __future__ import print_function, absolute_import, division 
import sys
import os
import csv
import datetime
import time
import argparse
from collections import defaultdict

WARNING_PERCENTAGE_ANNOTATED_CUTOFF = 50


def time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def log(message):
    print('{}|{}'.format(time_stamp(), message), file=sys.stderr)

def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description='Adds NCBI gene_info annotations to diffex output of scripts/Cuffdiff_out_format_v6.pl')

    parser.add_argument(
        '-i', '--geneinfo', type=str, help='path to input NCBI gene_info file', required=True)
    parser.add_argument(
        '-e', '--diffexp', type=str, help='path to input differential_expression_file', required=True)
    parser.add_argument(
        '-g', '--genome', type=str, help='organism taxonomy ID', required=True)  
    parser.add_argument(
        '-o', '--outdir', type=str, help='output directory name', required=True)
    args = parser.parse_args()
    return args


def check_geneinfo_matches(all_lines,
                           matching_lines,
                           warning_percent_annotated_cutoff=WARNING_PERCENTAGE_ANNOTATED_CUTOFF):
    percentage_matches = int((matching_lines/all_lines)*100)
    msg_format = 'Annotations added for {} ({}%) of rows in diffexp file'
    msg = msg_format.format(matching_lines, percentage_matches)
    log(msg)
    if percentage_matches < warning_percent_annotated_cutoff:
        msg_format = 'WARNING: Low number of annotations! (< {}% of genes were annotated)'
        msg = msg_format.format(warning_percent_annotated_cutoff)
        log(msg)


args = parse_command_line_args()
gene_info = args.geneinfo
gene_expr = args.diffexp
outdir = args.outdir

taxonomy = { 'mm10': '10090', 'hg19': '9606' }

def get_taxid(genome):
    if genome in taxonomy:
        return taxonomy.get(genome)
    else:
        msg_fmt = 'ERROR: Taxonomy id not found for -g/--genome:[{}]. Available taxonomy-genome mappings are: {}'
        msg = msg_fmt.format(genome, taxonomy)
        print(msg, file=sys.stderr)
        sys.exit()

tax_id = get_taxid(args.genome)

outfile_tag = os.path.basename(gene_expr.replace('.txt', '.annot.txt'))
outfile_name = os.path.join(outdir, outfile_tag)

log('reading entrez gene info')



#### def get_geneinfo_lines():
    
gene_details = defaultdict(list)
with open(gene_info,'r') as geneinfo_file:
    next(geneinfo_file) # skip header
    reader=csv.reader(geneinfo_file,delimiter='\t')
    for row in reader:
        (geneinfo_tax_id, geneinfo_gene_id, geneinfo_gene_symbol, geneinfo_gene_description) = (row[0], row[1], row[2], row[8])
        if geneinfo_tax_id == tax_id:
            gene_details[geneinfo_gene_symbol].append(geneinfo_gene_id)
            gene_details[geneinfo_gene_symbol].append(geneinfo_gene_description)


#### def add_annotation(args.diffexp, outfile_name):
all_lines = 0
matching_gene_symbol_count = 0
with open(gene_expr, 'r') as file_to_annotate, open(outfile_name, 'w') as annotated_file: 
    reader=csv.reader(file_to_annotate,delimiter='\t')
    for row in reader:
        all_lines += 1
        row = [col.strip() for col in row]
        if row[0] == 'test_id':
            row[1] = 'gene_symbol'
            row[2] = 'gene_id'
            row.insert(3, 'gene_desc')
            print('\t'.join(row), file=annotated_file)
        else:
            left = row[0:2]
            right = row[3:]
            gene_name = row[1]
            if gene_name in gene_details:
                gene_id_symbol = gene_details[gene_name]
                matching_gene_symbol_count += 1
            else:
                gene_id_symbol = ['.','.']
            outline = "\t".join(left + gene_id_symbol + right)  
            print(outline, file=annotated_file)

check_geneinfo_matches(all_lines, matching_gene_symbol_count)
log('done')
