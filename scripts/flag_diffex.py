#!/bin/env python

import argparse
import sys

import numpy as np
import pandas as pd

def parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(
        description='adds key columns to gene|isoform.diff file')

    parser.add_argument(
        '-f', '--foldchange', type=float, help='linear foldchange cut-off', required=True)
    parser.add_argument(
        'input_filepath', type=str, help='path to input differential_expression_file')
    parser.add_argument(
        'output_filepath', type=str, help='path to output file')

    args = parser.parse_args(sys_argv)
    return args 

def main(sys_argv):
    args = parse_command_line_args(sys_argv)
    input_df = pd.read_csv(args.input_filepath, sep='\t')
    input_df['significant'] = 1
    input_df['linear_fold_change'] = np.exp2(input_df['log2(fold_change)'])
    input_df['diff_exp'] = 1
    input_df.to_csv(args.output_filepath, index=False, sep='\t')

if __name__ == '__main__':
    main(sys.argv[1:])
