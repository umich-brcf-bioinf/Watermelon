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

TAXONOMY = {'hg19': '9606', 'GRCh37': '9606', 'hg38': '9606', 'GRCh38': '9606',   # human
            'mm10': '10090', 'mm9': '10090',    # mouse
            'rn5': '10116', 'rn6': '10116',     # rat
            'ce10': '6239', 'ce11': '6239', 'WS220': '6239', 'WBS235': '6239',   # c. elegans
            'GRCz10' : '7955', 'Zv9': '7955',   # zebra fish
            'TAIR9': '3702', 'TAIR10': '3702',  # arabidopis
            'MSU6': '4530',                     # rice
            'ecoMG1655': '511145', 'ecoUTI89': '364106', # ecoli
            'dm6': '7227', # fly
            }

def time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def log(message):
    print('{}|annotate_entrez_gene_info |{}'.format(time_stamp(), message), file=sys.stderr)

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

def get_taxid(genome):
    if genome in TAXONOMY:
        return TAXONOMY.get(genome)
    else:
        msg_fmt = 'ERROR: Taxonomy id not found for -g/--genome:[{}]. Available taxonomy-genome mappings are: {}'
        msg = msg_fmt.format(genome, TAXONOMY)
        print(msg, file=sys.stderr)
        sys.exit()

def main():
    args = parse_command_line_args()
    gene_info = args.geneinfo
    gene_expr = args.diffexp
    outdir = args.outdir

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
                gene_details[geneinfo_gene_symbol] = [geneinfo_gene_id, geneinfo_gene_description]

    #### def add_annotation(args.diffexp, outfile_name):
    all_lines = 0
    matching_gene_symbol_count = 0
    with open(gene_expr, 'r') as file_to_annotate, open(outfile_name, 'w') as annotated_file:
        reader=csv.reader(file_to_annotate,delimiter='\t')
        for row in reader:
            all_lines += 1
            row = [col.strip() for col in row]
            if row[0] == 'id':
                row[0] = 'gene_symbol'
                row.insert(1, 'gene_id')
                row.insert(2, 'gene_desc')
                print('\t'.join(row), file=annotated_file)
            elif row[0] == 'gene_id':
                row[0] = 'gene_symbol'
                row[1] = 'tx_id'
                row.insert(2, 'gene_id')
                row.insert(3, 'gene_desc')
                print('\t'.join(row), file=annotated_file)
            else:
                if len(row) == 8:
                    left = row[0:1]
                    right = row[1:]
                else:
                    left = row[0:2]
                    right = row[2:]

                gene_name = row[0]
                if gene_name in gene_details:
                    gene_id_symbol = gene_details[gene_name]
                    matching_gene_symbol_count += 1
                else:
                    gene_id_symbol = ['.','.']

                outline = "\t".join(left + gene_id_symbol + right)
                print(outline, file=annotated_file)

    check_geneinfo_matches(all_lines, matching_gene_symbol_count)
    log('done')

if __name__ == '__main__':
    main()
