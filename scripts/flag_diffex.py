#!/bin/env python

import argparse
import math
import sys

import numpy as np
import pandas as pd

FDR_CUTOFF = 0.05
LOG2_FC = 'log2(fold_change)'
LINEAR_FC = 'linear_fold_change'
STATUS = 'status'
STATUS_OK = 'OK'
Q_VALUE = 'q_value'
DIFF_EXP = 'diff_exp'

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
    log2_fold_change = math.log2(linear_fold_change)
    df = input_df.copy()
    df[LINEAR_FC] = np.exp2(df[LOG2_FC])
    status_ok = lambda r: r[STATUS] == STATUS_OK
    significant = lambda r: r[Q_VALUE] <= FDR_CUTOFF
    effect_size_ok = lambda r: abs(r[LOG2_FC]) >= log2_fold_change
    yes_no = lambda x: 'Yes' if x else 'No'
    diff_exp = lambda r: yes_no(status_ok(r) and significant(r) and effect_size_ok(r))
    df[DIFF_EXP] = df.apply(diff_exp, axis=1)
    return df

def _validate_required_fields(input_df, args):
    required_fields = set([LOG2_FC, Q_VALUE, STATUS])
    actual_fields = set(input_df.columns.values.tolist())
    missing_fields = sorted(required_fields - actual_fields)
    if missing_fields:
        msg_fmt = 'Input file [{}] is missing required field(s) [{}].'
        msg = msg_fmt.format(args.input_file, ','.join(missing_fields))
        raise ValueError(msg)

def _validate_log2fc_numeric(input_df, args):
    def is_not_float(x):
        is_float = True
        try:
            float(x)
        except:
            is_float = False
        return not is_float
    non_numeric_index = input_df[LOG2_FC].apply(is_not_float)
    non_numeric_log2fc = input_df[non_numeric_index][LOG2_FC]
    non_numeric_excerpt = [(str(value) + ':' + str(line + 2)) for line, value in non_numeric_log2fc.head(n=5).to_dict().items()]
    if len(non_numeric_log2fc):
        msg_fmt = ('Input file [{}]: {} log2(fold_change) value(s) '
                   'are not numeric: [{}]')
        msg = msg_fmt.format(args.input_file,
                             len(non_numeric_log2fc),
                             ','.join(non_numeric_excerpt))
        raise ValueError(msg)

def _validate_fold_change_threshold(input_df, args):
    if args.foldchange < 1:
        msg_fmt = 'Specified foldchange cutoff [{}] must be >= 1'
        msg = msg_fmt.format(args.foldchange)
        raise ValueError(msg)

def _validate_inputs(input_df, args):
    validations = [_validate_required_fields,
                   _validate_log2fc_numeric,
                   _validate_fold_change_threshold,]
    for validation in validations:
        validation(input_df, args)

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    df = pd.read_csv(args.input_filepath, sep='\t')
    df = _add_columns(df, args.foldchange)
    df.to_csv(args.output_filepath, index=False, sep='\t')

if __name__ == '__main__':
    main(sys.argv[1:])
