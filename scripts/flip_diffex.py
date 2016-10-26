#!/bin/env python
from __future__ import print_function, absolute_import, division

import argparse
import os
import sys

import pandas as pd

REQUIRED_FIELDS = argparse.Namespace(
    sample_1='sample_1',
    sample_2='sample_2',
    log2_fold_change='log2(fold_change)',
    status='value_1',
    significant='value_2',
)

def _log(msg):
    print(msg, file=sys.stderr)

def _flip_comparisons(dataframe, comparisons):
    _COMPARISON_NAME_FMT = '{}_{}'
    _COLUMN_NAMES = ['sample_1', 'sample_2', 'value_1', 'value_2', 'log2(fold_change)']

    def _flip_comparison(row):
        this_comparison = _COMPARISON_NAME_FMT.format(row['sample_2'], row['sample_1'])
        if this_comparison in comparisons:
            return row['sample_2'], row['sample_1'], row['value_2'], row['value_1'], -row['log2(fold_change)']
        else:
            return row['sample_1'], row['sample_2'], row['value_1'], row['value_2'], row['log2(fold_change)']

    dataframe[_COLUMN_NAMES] = dataframe.apply(_flip_comparison, axis=1).apply(pd.Series)


def _validate_required_fields(input_df, args):
    required_field_set = set(vars(REQUIRED_FIELDS).values())
    actual_fields = set(input_df.columns.values.tolist())
    missing_fields = sorted(required_field_set - actual_fields)
    if missing_fields:
        msg_fmt = 'Input file [{}] is missing required field(s) [{}].'
        msg = msg_fmt.format(args.input_filepath, ','.join(missing_fields))
        raise ValueError(msg)

def _validate_inputs(input_df, args):
    validations = [_validate_required_fields]
    for validation in validations:
        validation(input_df, args)

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(
        description=('sorts diffex file splitting into pairwise comparison files based '
                     'on distinct values of sample_1 and sample_2'))
    parser.add_argument(
        'input_filepath', type=str, help='path to input gene or isoform file')
    parser.add_argument(
        'output_filepath', type=str, help='path to output file')
    parser.add_argument(
        'comparisons',
        type=str,
        help=('commma separated list of comparisons; comparisons found in the input file but '
              'not in this list will be passed through unmodified'))

    args = parser.parse_args(sys_argv)
    return args 

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    _log('reading {}'.format(args.input_filepath))
    df = pd.read_csv(args.input_filepath, sep='\t')
    _validate_inputs(df, args)
    comparisons = args.comparisons.split(',')
    _log('flipping comparisons')
    _flip_comparisons(df, comparisons)
    _log('saving {}'.format(args.output_filepath))
    df.to_csv(args.output_filepath, index=False, sep='\t')
    _log('done')

if __name__ == '__main__':
    main(sys.argv[1:])
