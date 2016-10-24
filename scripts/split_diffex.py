#!/bin/env python
from __future__ import print_function, absolute_import, division

import argparse
from collections import defaultdict
import os
import sys

import pandas as pd

OUTPUT_FILE_EXTENSION = '.txt'

REQUIRED_FIELDS = argparse.Namespace(
    sample_1='sample_1',
    sample_2='sample_2',
    log2_fold_change='log2(fold_change)',
    status='status',
    significant='significant',
    diff_exp='diff_exp',
)

GROUP_BY_COLUMNS= [REQUIRED_FIELDS.sample_1, REQUIRED_FIELDS.sample_2]

COMPARISON_NAME_COLUMN = '_comparison_name'

def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)

class _FilteringHandler(object):
    def __init__(self, included_comparisons, handler, log):
        self._included_comparisons = included_comparisons
        self._base_handler = handler
        self._log = log
        self._stats = defaultdict(int)
        self._processed_count = 0

    def handle(self, comparison_name, comparison_df):
        if comparison_name in self._included_comparisons:
            self._stats['split'] += 1
            self._log('split comparison [{}]'.format(comparison_name))
            self._base_handler.handle(comparison_name, comparison_df)
        else:
            self._stats['skipped'] += 1
            self._log('skipped comparison [{}]'.format(comparison_name))

    def end(self):
        self._base_handler.end()
        total_count = sum([x for x in self._stats.values()])
        self._log('total comparisons: {}'.format(total_count))
        for stat_name, stat_count in self._stats.items():
            self._log('{} comparisons: {}'.format(stat_name, stat_count))

class _ComparisonHandler(object):
    def __init__(self, target_dir, suffix):
        self._target_dir = target_dir
        self._suffix = suffix

    def handle(self, comparison_name, comparison_df):
        output_filename = '{}{}'.format(comparison_name, self._suffix)
        output_filepath = os.path.join(self._target_dir, output_filename)
        comparison_df.to_csv(output_filepath, index=False, sep='\t')

    def end(self):
        pass

def _log(msg):
    print(msg, file=sys.stderr)

def _get_sort_value(row):
    yes_no = lambda x: 0 if x=='yes' else 1
    ok_no_test = lambda x: 0 if x=='OK' else 1
    status_sort = ok_no_test(row[REQUIRED_FIELDS.status])
    significant_sort = yes_no(row[REQUIRED_FIELDS.significant])
    diff_exp_sort = yes_no(row[REQUIRED_FIELDS.diff_exp])
    log2_fold_change = row[REQUIRED_FIELDS.log2_fold_change]
    return (status_sort, significant_sort, diff_exp_sort, -log2_fold_change)

def _get_unique_comparison_df(group_by_columns, df):
    def _build_comparison_name(row):
        values = []
        for cols in group_by_columns:
            values.append(row[cols])
        return '_'.join(values)
    unique_comparison_df = df[group_by_columns].drop_duplicates()
    unique_comparison_df[COMPARISON_NAME_COLUMN] = unique_comparison_df.apply(_build_comparison_name,
                                                                              axis=1)
    return unique_comparison_df

def _split_comparisons(df, group_by_columns, comparison_handler, log):
    unique_comparison_df = _get_unique_comparison_df(group_by_columns, df)
    total_comparisons = len(unique_comparison_df)
    comparison_count = 0
    for index, comparison in unique_comparison_df.iterrows():
        comparison_name = comparison[COMPARISON_NAME_COLUMN]
        comparison_count += 1
#        log('processing [{}] ({}/{})'.format(comparison_name, comparison_count, total_comparisons))
        comparison_df = df.copy()
        for cols in group_by_columns:
            in_comparison = comparison_df[cols] == comparison[cols]
            comparison_df = comparison_df[in_comparison]
        comparison_handler.handle(comparison_name, comparison_df)
    comparison_handler.end()
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
        msg = msg_fmt.format(args.input_filepath, ','.join(missing_fields))
        raise ValueError(msg)

def _validate_included_comparisons_present(input_df, args):
    included_comparisons = set(args.included_comparisons.split(','))
    unique_comparison_df = _get_unique_comparison_df(GROUP_BY_COLUMNS, input_df)
    found_comparisons = set(unique_comparison_df[COMPARISON_NAME_COLUMN])
    missing_comparisons = sorted(included_comparisons - found_comparisons)
    if missing_comparisons:
        msg_fmt = 'Input file [{}] is missing requested comparison(s) [{}].'
        msg = msg_fmt.format(args.input_filepath, ','.join(missing_comparisons))
        raise ValueError(msg)

def _validate_inputs(input_df, args):
    validations = [_validate_required_fields, _validate_included_comparisons_present]
    for validation in validations:
        validation(input_df, args)

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(
        description=('sorts diffex file splitting into pairwise comparison files based '
                     'on distinct values of sample_1 and sample_2'))
    parser.add_argument(
        'input_filepath', type=str, help='path to input gene or isoform file')
    parser.add_argument(
        'output_dir', type=str, help='path to existing output dir')
    parser.add_argument(
        'included_comparisons',
        type=str,
        help=('commma separated list of comparisons; comparisons found in the input but '
              'not in this list will be ignored'))

    args = parser.parse_args(sys_argv)
    return args 

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    _log('reading {}'.format(args.input_filepath))
    df = pd.read_csv(args.input_filepath, sep='\t')
    _validate_inputs(df, args)
    _log('sorting')
    df = _sort(df)
    comparison_handler = _ComparisonHandler(args.output_dir, OUTPUT_FILE_EXTENSION)
    included_comparison_list = args.included_comparisons.split(',')
    filtering_handler = _FilteringHandler(included_comparison_list,
                                          comparison_handler,
                                          _log)
    _mkdir(args.output_dir)
    _split_comparisons(df, GROUP_BY_COLUMNS, filtering_handler, _log)
    _log('done')

if __name__ == '__main__':
    main(sys.argv[1:])
