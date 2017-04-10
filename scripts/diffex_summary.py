#!/bin/env python
from __future__ import print_function, absolute_import, division

import argparse
from collections import OrderedDict
import datetime
import os
import re
import sys
import time

import pandas as pd

DEFAULT_ANNOTATION_COLUMN = 'gene_id'
DEFAULT_ANNOTATION_NULL = '.'
DEFAULT_TRIM_SUFFIX = '.annot.txt'

def _time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def _log(message):
    print('{}|diffex_summary|{}'.format(_time_stamp(), message), file=sys.stderr)

def _get_stats(args, df, name):
    count_annotated = _count_annotated(args.annotation_column,
                                                 args.annotation_null,
                                                 df)
    count_diffex_pass = _count_diffex_pass(args.diffex_call_column,
                                           args.diffex_call_pass,
                                           df)
    percent_annotated = round(100.0 * count_annotated / len(df), 2)
    stats = OrderedDict([('comparison', name),
                         ('total_count', len(df)),
                         ('count_diff_expressed', count_diffex_pass),
                         ('count_annotated', count_annotated),
                         ('percent_annotated', percent_annotated),])
    return stats

def _count_diffex_pass(diffex_call_column, diffex_call_pass, df):
    return sum(df[diffex_call_column] == diffex_call_pass)

def _count_annotated(annotation_column, annotation_null, df):
    return sum(df[annotation_column] != annotation_null)

def _simplify_file_names(file_paths, suffix):
    simplified_names = {}
    common_prefix = os.path.commonprefix(file_paths)
    for path in file_paths:
        name = path
        if suffix:
            name = re.sub(suffix+'$', '', name)
        if common_prefix:
            name = re.sub('^'+common_prefix, '', name)
        simplified_names[path] = name
    return simplified_names 

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(
        description=('accepts a list of tab-separated files and emits a tab-separated summary of DE and annotated genes counts for each input file'))
    parser.add_argument(
        '--annotation_column',
        type=str,
        help='={} column to check for annotation values'.format(DEFAULT_ANNOTATION_COLUMN),
        default=DEFAULT_ANNOTATION_COLUMN)
    parser.add_argument(
        '--annotation_null',
        type=str,
        help='={} column value indicates gene was not annotated'.format(DEFAULT_ANNOTATION_NULL),
        default=DEFAULT_ANNOTATION_NULL)
    parser.add_argument(
        '--diffex_call_column',
        type=str,
        required=True,
        help='column to check for diffex call passed')
    parser.add_argument(
        '--diffex_call_pass',
        type=str,
        required=True,
        help='column value indicates diffex call passed')
    parser.add_argument(
        '--output_file',
        type=str,
        required=True,
        help='output path and filename')
    parser.add_argument(
        '--output_xlsx',
        type=str,
        help='optional Excel output filename')
    parser.add_argument(
        '--trim_suffix',
        type=str,
        help='={} Output will remove this string when listing the file names'.format(DEFAULT_TRIM_SUFFIX),
        default=DEFAULT_TRIM_SUFFIX)
    parser.add_argument(
        'input_files',
        type=str,
        help='paths to input files',
        nargs='+')

    args = parser.parse_args(sys_argv)
    return args 

def main(sys_argv, log=_log):
    args = _parse_command_line_args(sys_argv)
    log('summarizing {} diffex files'.format(len(args.input_files)))
    names = _simplify_file_names(args.input_files, args.trim_suffix)
    all_stats = []
    for input_file in sorted(args.input_files):
        input_df = pd.read_table(input_file, header=0)
        stats = _get_stats(args, input_df, names[input_file])
        all_stats.append(stats)
    output_df = pd.DataFrame(all_stats, columns=stats.keys())
    output_df.to_csv(args.output_file, sep='\t', index=False)
    if args.output_xlsx:
        output_df.to_excel(args.output_xlsx, index=False)
    log('done')

if __name__ == '__main__':
    main(sys.argv[1:])
