#!/usr/bin/env/python

import argparse
import datetime
import sys
import time
import os

import pandas as pd

def time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def log(message):
    print('{}| annotate_gene_info |{}'.format(time_stamp(), message), file=sys.stderr)

def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description='Adds EntrezGeneIDs and Descrtiptions to diffex annotated_source_df')

    parser.add_argument(
        '-m', '--mapping_file', type=str, help='Path to mapping file with required columns "gene_id", "entrezgene_id", "external_gene_name", "description"', required=True)
    parser.add_argument(
        '-i', '--input_file', type=str, help='Path to input file', required=True)
    parser.add_argument(
        '-o', '--annotated_source_df_filename', type=str, help='Path to annotated annotated_source_df', required=True)
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
    #The mapping table should not have duplicate values in the index
    #The onus is on whoever creates the mapping table to do it correctly
    #There are instances of biomaRt queries returning multiple results
    #These are cases where a query ID matches equally well to multiple
    #items from another database, so biomaRt returns them all. Current solution
    #is to have the conflicting attributes separated by commas
    try:
        annotation_mapping_df = pd.read_csv(args.mapping_file, sep='\t', dtype=object, index_col=args.mapping_idx)
    except ValueError:
        msg = '{} is not a column in {}'.format(args.mapping_idx, args.mapping_file)
        raise ValueError(msg)

    # Read in input file
    try:
        source_df = pd.read_csv(args.input_file, sep='\t', dtype=object, index_col=args.input_idx)
    except ValueError:
        msg = '{} is not a column in {}'.format(args.input_idx, args.input_file)
        raise ValueError(msg)

    #Add annotation columns
    annotated_source_df = annotation_mapping_df.join(source_df, how='right').reset_index()

    #Remove redundant external_gene_name column if it's identical to gene_id
    if annotated_source_df[args.mapping_idx].equals(annotated_source_df['external_gene_name']):
        annotated_source_df.drop(columns='external_gene_name', inplace=True)

    annotated_source_df.to_csv(args.annotated_source_df_filename, sep='\t', na_rep='.', index=False)

    log('done')

if __name__ == '__main__':
    main()
