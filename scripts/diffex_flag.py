#!/bin/env python
from __future__ import print_function, absolute_import, division

import argparse
import datetime
import math
import sys
import time

import numpy as np
import pandas as pd

FDR_CUTOFF = 0.05
STATUS_OK = 'OK'

REQUIRED_FIELDS = argparse.Namespace(
    log2_fold_change='log2(fold_change)',
    status='status',
    q_value='q_value',
)

NEW_FIELDS = argparse.Namespace(
    linear_fold_change = 'linear_fold_change',
    diff_exp = 'diff_exp',
)

def _time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def _log(message):
    print('{}|flag_diffex|{}'.format(_time_stamp(), message), file=sys.stderr)

def _add_columns(input_df, linear_fold_change):
    log2_fold_change = math.log(linear_fold_change, 2)
    df = input_df.copy()
    df[NEW_FIELDS.linear_fold_change] = np.exp2(df[REQUIRED_FIELDS.log2_fold_change])
    status_ok = lambda r: r[REQUIRED_FIELDS.status] == STATUS_OK
    significant = lambda r: r[REQUIRED_FIELDS.q_value] <= FDR_CUTOFF
    effect_size_ok = lambda r: abs(r[REQUIRED_FIELDS.log2_fold_change]) >= log2_fold_change
    yes_no = lambda x: 'Yes' if x else 'No'
    diff_exp = lambda r: yes_no(status_ok(r) and significant(r) and effect_size_ok(r))
    df[NEW_FIELDS.diff_exp] = df.apply(diff_exp, axis=1)
    return df

def _validate_required_fields(input_df, args):
    required_fields = set(vars(REQUIRED_FIELDS).values())
    actual_fields = set(input_df.columns.values.tolist())
    missing_fields = sorted(required_fields - actual_fields)
    if missing_fields:
        msg_fmt = 'Input file [{}] is missing required field(s) [{}].'
        msg = msg_fmt.format(args.input_filepath, ','.join(missing_fields))
        raise ValueError(msg)

def _validate_log2fc_numeric(input_df, args):
    def is_not_float(x):
        is_float = True
        try:
            float(x)
        except:
            is_float = False
        return not is_float
    non_numeric_index = input_df[REQUIRED_FIELDS.log2_fold_change].apply(is_not_float)
    non_numeric_log2fc = input_df[non_numeric_index][REQUIRED_FIELDS.log2_fold_change]
    non_numeric_excerpt = [(str(value) + ':' + str(line + 2)) for line, value in non_numeric_log2fc.head(n=5).to_dict().items()]
    if len(non_numeric_log2fc):
        msg_fmt = ('Input file [{}]: {} log2(fold_change) value(s) '
                   'are not numeric: [{}]')
        msg = msg_fmt.format(args.input_filepath,
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

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    _log('reading {}'.format(args.input_filepath))
    df = pd.read_csv(args.input_filepath, sep='\t')
    _validate_inputs(df, args)
    _log('adding flags')
    df = _add_columns(df, args.foldchange)
    _log('saving {}'.format(args.output_filepath))
    df.to_csv(args.output_filepath, index=False, sep='\t')

if __name__ == '__main__':
    main(sys.argv[1:])
