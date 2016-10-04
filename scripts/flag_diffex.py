#!/bin/env python

import argparse
import math
import sys

import numpy as np
import pandas as pd

FDR_CUTOFF = 0.05

def _parse_command_line_args(sys_argv):
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

def _add_columns(input_df, linear_fold_change):
    log_fold_change = math.log2(linear_fold_change)
    df = input_df.copy()
    df['linear_fold_change'] = np.exp2(df['log2(fold_change)'])
    status_ok = lambda r: r['status'] == 'OK'
    significant = lambda r: r['qvalue'] <= FDR_CUTOFF
    effect_size_ok = lambda r: abs(r['log2(fold_change)']) >= log_fold_change
    yes_no = lambda x: 'Yes' if x else 'No'
    diff_exp = lambda r: yes_no(status_ok(r) and significant(r) and effect_size_ok(r))
    df['diff_exp'] = df.apply(diff_exp, axis=1)
    return df

def _validate_inputs(input_df, args):
    required_fields = set(['log2(fold_change)', 'qvalue', 'status'])
    actual_fields = set(input_df.columns.values.tolist())
    missing_fields = sorted(required_fields - actual_fields)
    if missing_fields:
        msg = 'Input file [{}] is missing required field(s) [{}].'
        raise ValueError(msg.format(args.input_file, ','.join(missing_fields)))

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    df = pd.read_csv(args.input_filepath, sep='\t')
    df = _add_columns(df, args.foldchange)
    df.to_csv(args.output_filepath, index=False, sep='\t')

if __name__ == '__main__':
    main(sys.argv[1:])
