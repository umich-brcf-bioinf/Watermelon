#!/usr/bin/env python
'''Accepts list of tab-delim text summary files and merges into excel and text output files.
'''
from __future__ import print_function, absolute_import, division
import argparse
from os.path import basename, splitext
import sys

import pandas as pd

_DEFAULT_SUMMARY_BASENAME = 'combined_summary'
_SOURCE_COLUMN = 'source'

def _log(msg):
    print(msg, file=sys.stderr)

def _base_name(file_path):
    return splitext(basename(file_path))[0]

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(
        description='Accepts list of summary files and merges to Excel and text files.')

    parser.add_argument(
        '-o', '--output_base',
        type=str,
        default=_DEFAULT_SUMMARY_BASENAME,
        help='={} path and basename of output files'.format(_DEFAULT_SUMMARY_BASENAME),
        required=False)
    parser.add_argument(
        'input_summary_paths',
        nargs='+',
        type=str,
        help='path(s) to input summaries (tab-delim text files)')
    args = parser.parse_args(sys_argv)
    return args

def main(sys_argv, log=_log):
    args = _parse_command_line_args(sys_argv)
    summaries = {}
    headers = set()
    combined_df = None
    for file_name in sorted(args.input_summary_paths):
        log('reading: {}'.format(file_name))
        base_name = _base_name(file_name)
        df = pd.read_table(file_name, delimiter="\t")
        if len(df) < 1:
            log('WARNING: file [{}] is empty'.format(file_name))

        df.insert(0, _SOURCE_COLUMN, value=base_name)
        headers.add('|'.join(df.columns))
        summaries[base_name] = df
        combined_df = pd.concat([combined_df, df], ignore_index=True, sort=True)

    if len(headers) > 1:
        log('WARNING: incoming headers were not consistent')
    else:
        combined_df = combined_df[list(df.columns)]

    combined_df = combined_df.fillna('')
    cols = list(combined_df)
    cols.insert(0, cols.pop(cols.index(_SOURCE_COLUMN)))
    combined_df = combined_df.ix[:, cols]

    file_name = args.output_base + '.txt'
    combined_df.to_csv(file_name, index=False, sep='\t')
    log('created combined file: {}'.format(file_name))

    file_name = args.output_base + '.xlsx'
    writer = pd.ExcelWriter(file_name)
    combined_df.to_excel(writer,
                         sheet_name=_base_name(args.output_base),
                         index=False)
    for name, df in sorted(summaries.items()):
        df.to_excel(writer, sheet_name=name, index=False)
    writer.save()
    log('created combined file: {}'.format(file_name))
    log('done')

if __name__ == '__main__':
    main(sys.argv[1:])
