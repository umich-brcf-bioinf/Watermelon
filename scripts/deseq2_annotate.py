#!/usr/bin/python

from __future__ import print_function, absolute_import, division
import sys
import os
import datetime
import time
import argparse
import pandas as pd

def time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def log(message):
    print('{}| annotate_gene_info |{}'.format(time_stamp(), message), file=sys.stderr)

def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description='Adds EntrezGeneIDs and Descrtiptions to diffex output')

    parser.add_argument(
        '-m', '--mapping_file', type=str, help='Path to mapping file with required columns "gene_id", "entrezgene_id", "external_gene_name", "description"', required=True)
    parser.add_argument(
        '-i', '--input_file', type=str, help='Path to input file', required=True)
    parser.add_argument(
        '-o', '--output_filename', type=str, help='Path to annotated output', required=True)
    parser.add_argument(
        '--mapping_idx', type=str, help='String matching the column to use as the index of the mapping file', required=True)
    parser.add_argument(
        '--input_idx', type=str, help='String matching the column to use as the index of the input file', required=True)
    args = parser.parse_args()
    return args

def main():
    args = parse_command_line_args()

    log('Reading ID mapping info')

    # Read in mapping info
    attr_table = pd.read_csv(args.mapping_file, sep="\t", dtype=object)
    try:
        #TWS FIXME? - There are instances of biomaRt queries returning multiple results
        #These are cases where a query ID matches equally well to multiple items from another database, so biomaRt returns them all
        #Right now, just keeping the first occurence - default of drop_duplicates
        attr_table = attr_table.drop_duplicates(args.mapping_idx).set_index(args.mapping_idx)
    except KeyError:
        print("{} is not a column in {}".format(args.mapping_idx, args.mapping_file))
        sys.exit(os.EX_CONFIG)

    # Read in input file
    try:
        to_annotate = pd.read_csv(args.input_file, sep="\t", dtype=object, index_col=args.input_idx)
    except ValueError:
        print("{} is not a column in {}".format(args.input_idx, args.input_file))
        sys.exit(os.EX_CONFIG)

    #Add annotation columns
    output = attr_table.join(to_annotate).reset_index()

    #Remove redundant external_gene_name column if it's identical to gene_id
    if output[args.mapping_idx].equals(output['external_gene_name']):
        output.drop(columns='external_gene_name', inplace=True)

    output.to_csv(args.output_filename, sep="\t", na_rep="NA", index=False)

    log('done')

if __name__ == '__main__':
    main()
