#!/bin/env python
from __future__ import print_function, absolute_import, division

import argparse
import datetime
import glob
import os
import sys
import time

import numpy as np
import pandas as pd

COUNTS_INDEX_HEADER = 'GeneId'
STAT_COLUMNS = ['__ambiguous', '__no_feature', '__too_low_aQual', '__not_aligned', '__alignment_not_unique']
STATS_INDEX_HEADER = 'readcounts'
STATS_FEATURE_ASSIGNED = 'feature_assigned'
STATS_TOTAL = 'total'

DEFAULT_WORKING_DIR = os.getcwd()
DEFAULT_SUFFIX = '_counts.txt'
DEFAULT_COUNTS_FILENAME = 'htseq_merged.txt'
DEFAULT_STATS_FILENAME = 'htseq_merged_stats.txt'

DESCRIPTION = \
'''Merges individual sample htseq count files into a single, merged file with
a row for every gene and a column for every sample. HTSeq stats are saved into
a separate file.'''

def _time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def _log(message, *args):
    if args:
        message = message.format(*[i for i in args])
    print('{}|deseq2_htseq_merge|{}'.format(_time_stamp(), message), file=sys.stderr)

def _build_sample_dataframes(args, log=_log):
    sample_count_files = glob.glob(os.path.join(args.htseq_dir,
                                                '*' + args.suffix))
    if len(sample_count_files) == 0:
        msg = ('Found no htseq sample files (working_dir=[{}], suffix=[{}]); '
               'review inputs and try again').format(args.htseq_dir,
                                                     args.suffix)
        raise ValueError(msg)
    log('read {} sample count files', len(sample_count_files))

    strip_suffix = lambda x: os.path.basename(x[:x.rfind(args.suffix)])
    htseq_file_to_df = lambda x: pd.read_table(x,
                                               sep='\t',
                                               header=None,
                                               index_col=0,
                                               names=[COUNTS_INDEX_HEADER,
                                                      strip_suffix(x)])
    return list(map(htseq_file_to_df, sample_count_files))

def _merge_dataframes(sample_data_frames):
    df = pd.concat(sample_data_frames, axis=1)

    counts_df = df.loc[~df.index.isin(STAT_COLUMNS)].copy()
    counts_df.index.name=COUNTS_INDEX_HEADER
    counts_df.fillna(value=0, inplace=True)
    counts_df = counts_df.applymap(np.int64)

    stats_df = df.loc[df.index.isin(STAT_COLUMNS)].copy()
    stats_df.loc[STATS_FEATURE_ASSIGNED,:] = counts_df.sum()
    stats_df.loc[STATS_TOTAL,:] = stats_df.sum()
    stats_df.fillna(value=0, inplace=True)
    stats_df = stats_df.applymap(np.int64)
    stats_df.index.name=STATS_INDEX_HEADER

    return counts_df, stats_df

def _gene_count_from_sample_df(df):
    return len(set(df.index) - set(STAT_COLUMNS))

def _validate_dataframes(sample_dfs, counts_df, stats_df):
    expected_sample_count = len(sample_dfs)
    actual_sample_count = counts_df.shape[1]
    if expected_sample_count != actual_sample_count:
        msg = ('Expected [{}] sample columns in merged file, but '
               'found [{}]').format(expected_sample_count, actual_sample_count)
        raise ValueError(msg)

    min_row_count = _gene_count_from_sample_df(sample_dfs[0])
    print(min_row_count, len(counts_df))
    if  min_row_count > len(counts_df):
        msg = ('Expected at least [{}] rows in merged file, but '
               'found [{}]').format(min_row_count, len(counts_df))
        raise ValueError(msg)

    if counts_df.shape[1] != stats_df.shape[1]:
        msg = ('Expected columns in counts [{}] to match columns in stats [{}]'
               ).format(counts_df.shape[1], stats_df.shape[1])
        raise ValueError(msg)

    if len(stats_df) < 2 :
        msg = ('Expected at least 2 rows in stats file, but '
               'found [{}]').format(len(stats_df))
        raise ValueError(msg)

def _save_merged_dataframes(args, counts_df, stats_df, log=_log):
    counts_df.to_csv(args.counts_filename, sep='\t')
    stats_df.to_csv(args.stats_filename, sep='\t')
    log('merged [{}] samples to count file [{}] and stats file [{}]',
         counts_df.shape[1],
         args.counts_filename,
         args.stats_filename)

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--htseq_dir',
        type=str,
        required=False,
        help='={}; path to individual htseq count files'.format(DEFAULT_WORKING_DIR),
        default=DEFAULT_WORKING_DIR)
    parser.add_argument(
        '--suffix',
        type=str,
        required=False,
        help=('={}; simple text suffix that identifies HTSeq sample '
              'files').format(DEFAULT_SUFFIX),
        default=DEFAULT_SUFFIX)
    parser.add_argument(
        '--counts_filename',
        type=str,
        required=False,
        help='={}; output of merged counts '.format(DEFAULT_COUNTS_FILENAME),
        default=DEFAULT_COUNTS_FILENAME)
    parser.add_argument(
        '--stats_filename',
        type=str,
        required=False,
        help='={}; output of merged stats '.format(DEFAULT_STATS_FILENAME),
        default=DEFAULT_STATS_FILENAME)
    args = parser.parse_args(sys_argv)
    args.htseq_dir = os.path.join(args.htseq_dir, '')
    return args

def main(sys_argv, log=_log):
    args = _parse_command_line_args(sys_argv)
    sample_data_frames = _build_sample_dataframes(args, log)
    counts_df, stats_df = _merge_dataframes(sample_data_frames)
    _validate_dataframes(sample_data_frames, counts_df, stats_df)
    _save_merged_dataframes(args, counts_df, stats_df, log)
    _log('done')

if __name__ == '__main__':
    main(sys.argv[1:])
