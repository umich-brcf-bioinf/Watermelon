#!/bin/env python
from __future__ import print_function, absolute_import, division

import argparse
import datetime
import glob
import os
import time

import numpy as np
import pandas as pd

COUNT_INDEX_HEADER = 'GeneId'
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
    message = message.format(args)
    print('{}|deseq2_htseq_merge|{}'.format(_time_stamp(), message), file=sys.stderr)


def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--htseq_dir',
        type=str,
        required=True,
        help='= {}; path to individual htseq count files'.format(DEFAULT_WORKING_DIR),
        default=DEFAULT_WORKING_DIR)
    parser.add_argument(
        '--suffix',
        type=str,
        required=False,
        help=('= {}; simple text suffix that identifies HTSeq sample '
              'files').format(DEFAULT_SUFFIX),
        default=DEFAULT_SUFFIX)
    parser.add_argument(
        '--counts_filename',
        type=str,
        required=False,
        help=('= {}; output of merged counts '.format(DEFAULT_COUNTS_FILENAME),
        default=DEFAULT_COUNTS_FILENAME)
    parser.add_argument(
        '--stats_filename',
        type=str,
        required=False,
        help=('= {}; output of merged stats '.format(DEFAULT_STATS_FILENAME),
        default=DEFAULT_STATS_FILENAME)

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)

    strip_suffix = lambda x: os.path.basename(x[:x.rfind(args.suffix)])
    htseq_file_to_df = lambda x: pd.read_table(x,
                                               sep='\t',
                                               header=None,
                                               index_col=0,
                                               names=[COUNT_INDEX_HEADER,
                                                      strip_suffix(x)])
    sample_count_files = glob.glob(os.path.join(args.htseq_dir,
                                                '*' + args.suffix))
    if len(sample_count_files) == 0:
        msg = ('Found no htseq sample files (working_dir=[{}], suffix=[{}]); '
               'review inputs and try again').format(args.htseq_dir,
                                                     args.suffix)
        raise ValueError(msg)
    _log('read {} sample count files', len(sample_count_files))
    df = pd.concat(map(htseq_file_to_df, sample_count_files), axis=1)

    counts = df.loc[~df.index.isin(STAT_COLUMNS)].copy()

    stats = df.loc[df.index.isin(STAT_COLUMNS)].copy()
    stats.loc[STATS_FEATURE_ASSIGNED,:] = counts.sum()
    stats.loc[STATS_TOTAL,:] = stats.sum()
    stats = stats.applymap(np.int64)
    stats.index.name=STATS_INDEX_HEADER

    expected_sample_count = len(sample_count_files)
    actual_sample_count = counts.shape[1]
    if expected_sample_count != actual_sample_count:
        msg = ('Expected [{}] sample columns in merged file, but '
               'found [{}]').format(expected_sample_count, actual_sample_count)
        raise ValueError(msg)

    counts.to_csv(args.counts_filename, sep='\t')
    stats.to_csv(args.stats_filename, sep='\t')
    _log(('merged [{}] sample count files to count file '
          '[{}] and stats file [{}]'),
          len(sample_count_files),
          args.counts_filename,
          args.stats_filename)
    _log('done')

if __name__ == '__main__':
    main(sys.argv[1:])
