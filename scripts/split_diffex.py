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
    status='status',
    significant='significant',
    diff_exp='diff_exp',
)


class _GroupHandler(object):
    def __init__(self, target_dir, suffix):
        self._target_dir = target_dir
        self._suffix = suffix

    def handle(self, group_name, group_df):
        output_filename = '{}{}'.format(group_name, self._suffix)
        output_filepath = os.path.join(self._target_dir, output_filename)
        group_df.to_csv(output_filepath, index=False, sep='\t')


def _log(msg):
    print(msg, file=sys.stderr)

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(
        description=('sorts diffex file splitting into pairwise comparison files based '
                     'on distinct values of sample_1 and sample_2'))
    parser.add_argument(
        'input_filepath', type=str, help='path to input gene or isoform file')
    parser.add_argument(
        'output_dir', type=str, help='path to existing output dir')

    args = parser.parse_args(sys_argv)
    return args 

def _get_sort_value(row):
    yes_no = lambda x: 0 if x=='yes' else 1
    ok_no_test = lambda x: 0 if x=='OK' else 1
    status_sort = ok_no_test(row[REQUIRED_FIELDS.status])
    significant_sort = yes_no(row[REQUIRED_FIELDS.significant])
    diff_exp_sort = yes_no(row[REQUIRED_FIELDS.diff_exp])
    log2_fold_change = row[REQUIRED_FIELDS.log2_fold_change]
    return (status_sort, significant_sort, diff_exp_sort, -log2_fold_change)

def _split_groups_by_sample(df, group_handler, log):
    group_by_cols = [REQUIRED_FIELDS.sample_1, REQUIRED_FIELDS.sample_2]

    def concat_values(group):
        values = []
        for cols in group_by_cols:
            values.append(group[cols])
        return '_'.join(values)

    unique_groups_df = df[group_by_cols].drop_duplicates()
    total_groups = len(unique_groups_df)
    group_count = 0
    for index, group in unique_groups_df.iterrows():
        group_name = concat_values(group)
        group_count += 1
        log('processing [{}] ({}/{})'.format(group_name, group_count, total_groups))
        in_group = (df[REQUIRED_FIELDS.sample_1] == group[REQUIRED_FIELDS.sample_1]) & \
                   (df[REQUIRED_FIELDS.sample_2] == group[REQUIRED_FIELDS.sample_2])
        group_df = df[in_group]
        group_handler(group_name, group_df)

def _sort(input_df):
    df = pd.DataFrame(input_df)
    df['_sort_value'] = input_df.apply(_get_sort_value, axis=1)
    df.sort_values('_sort_value', axis=0, inplace=True)
    df = df = df.drop('_sort_value', 1)
    return df

def _validate_required_fields(input_df, args):
    required_field_set = set(vars(REQUIRED_FIELDS).values())
    actual_fields = set(input_df.columns.values.tolist())
    missing_fields = sorted(required_field_set - actual_fields)
    if missing_fields:
        msg_fmt = 'Input file [{}] is missing required field(s) [{}].'
        msg = msg_fmt.format(args.input_file, ','.join(missing_fields))
        raise ValueError(msg)

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    _log('reading {}'.format(args.input_filepath))
    df = pd.read_csv(args.input_filepath, sep='\t')
    _validate_required_fields(df, args)
    _log('sorting')
    df = _sort(df)
    group_handler = _GroupHandler(args.output_dir, '.csv')
    _split_groups_by_sample(df, group_handler.handle, _log)
    _log('done')

if __name__ == '__main__':
    main(sys.argv[1:])
