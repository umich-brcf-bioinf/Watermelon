#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division
from argparse import Namespace
import os
from os.path import join
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd
from pandas.util.testing import assert_frame_equal
from testfixtures.tempdirectory import TempDirectory

import scripts.deseq2_htseq_merge as deseq2_htseq_merge

class Deseq2HtseqMergeTest(unittest.TestCase):
    def test_parse_command_line_args_defaults(self):
        command_line_args = []
        args = deseq2_htseq_merge._parse_command_line_args(command_line_args)
        self.assertEqual(os.path.join(os.getcwd(), ''), args.htseq_dir)
        self.assertEqual('_counts.txt', args.suffix)
        self.assertEqual('htseq_merged.txt', args.counts_filename)
        self.assertEqual('htseq_merged_stats.txt', args.stats_filename)
        self.assertEqual(len(vars(args)), 4)

    def test_parse_command_line_args_explicitValues(self):
        command_line_args = ['--htseq_dir', '/dir',
                             '--suffix', '_bar.baz',
                             '--counts_filename', 'hoopy.txt',
                             '--stats_filename', 'frood.txt']
        args = deseq2_htseq_merge._parse_command_line_args(command_line_args)
        self.assertEqual('/dir/', args.htseq_dir)
        self.assertEqual('_bar.baz', args.suffix)
        self.assertEqual('hoopy.txt', args.counts_filename)
        self.assertEqual('frood.txt', args.stats_filename)
        self.assertEqual(len(vars(args)), 4)

    def test_merge_dataframes_counts(self):
        s1_df = pd.DataFrame(data={'s1':1}, index=['gene1'])
        s2_df = pd.DataFrame(data={'s2':2}, index=['gene1','gene2'])
        s3_df = pd.DataFrame(data={'s3':4}, index=['gene2','gene3'])
        dfs = [s3_df, s1_df, s2_df]

        counts_df, _ = deseq2_htseq_merge._merge_dataframes(dfs)

        df_contents = StringIO(\
'''GeneId|s1|s2|s3
gene1|1|2|0
gene2|0|2|4
gene3|0|0|4''')
        expected_df =  pd.read_csv(df_contents, sep='|', index_col=0)
        assert_frame_equal(expected_df, counts_df, check_dtype=False)

    def test_merge_dataframes_countsIgnoresStats(self):
        s1_df = pd.DataFrame(data={'s1':1}, index=['gene1', '__ambiguous'])
        s2_df = pd.DataFrame(data={'s2':2}, index=['gene1','gene2', '__no_feature'])
        s3_df = pd.DataFrame(data={'s3':4}, index=['gene2','gene3', '__not_aligned'])
        dfs = [s1_df, s2_df, s3_df]

        counts_df, _ = deseq2_htseq_merge._merge_dataframes(dfs)

        df_contents = StringIO(\
'''GeneId|s1|s2|s3
gene1|1|2|0
gene2|0|2|4
gene3|0|0|4''')
        expected_df =  pd.read_csv(df_contents, sep='|', index_col=0)
        assert_frame_equal(expected_df, counts_df, check_dtype=False)

    def test_merge_dataframes_stats(self):
        s1_df = pd.DataFrame(data={'s1':1}, index=['gene1', '__ambiguous'])
        s2_df = pd.DataFrame(data={'s2':2}, index=['gene1','gene2', '__no_feature'])
        s3_df = pd.DataFrame(data={'s3':4}, index=['gene2','gene3', '__not_aligned'])
        dfs = [s1_df, s2_df, s3_df]

        _, stats_df = deseq2_htseq_merge._merge_dataframes(dfs)

        df_contents = StringIO(\
'''readcounts|s1|s2|s3
__ambiguous|1|0|0
__no_feature|0|2|0
__not_aligned|0|0|4
feature_assigned|1|4|8
total|2|6|12''')
        expected_df =  pd.read_csv(df_contents, sep='|', index_col=0)
        assert_frame_equal(expected_df, stats_df, check_dtype=False)

    def test_validate_dataframes_ok(self):
        sample1_df = pd.DataFrame(data={'s1':1}, index=['gene1','gene2'])
        sample2_df = pd.DataFrame(data={'s2':1}, index=['gene1','gene2'])
        counts_df = pd.DataFrame(data={'s1':1, 's2':2}, index=['gene1','gene2'])
        stats_df = pd.DataFrame(data={'s1':1, 's2':2}, index=range(7))
        deseq2_htseq_merge._validate_dataframes([sample1_df, sample2_df],
                                                counts_df,
                                                stats_df)
        self.assertEqual(1,1)

    def test_validate_dataframes_wrongCountsRows(self):
        sample1_df = pd.DataFrame(data={'s1':1}, index=['gene1','gene2'])
        sample2_df = pd.DataFrame(data={'s2':1}, index=['gene1','gene2'])
        counts_df = pd.DataFrame(data={'s1':1, 's2':2}, index=['gene1'])
        stats_df = pd.DataFrame(data={'s1':1, 's2':2}, index=range(7))
        self.assertRaisesRegexp(ValueError,
                                'Expected.*but found',
                                deseq2_htseq_merge._validate_dataframes,
                                [sample1_df, sample2_df],
                                counts_df,
                                stats_df)

    def test_validate_dataframes_wrongCountColumns(self):
        sample1_df = pd.DataFrame(data={'s1':1}, index=['gene1','gene2'])
        counts_df = pd.DataFrame(data={'s1':1, 's2':2, 's3':3}, index=['gene1','gene2'])
        stats_df = pd.DataFrame(data={'s1':1, 's2':2}, index=range(7))
        self.assertRaisesRegexp(ValueError,
                                'Expected.*but found',
                                deseq2_htseq_merge._validate_dataframes,
                                [sample1_df],
                                counts_df,
                                stats_df)

    def test_validate_dataframes_wrongStatsRows(self):
        sample1_df = pd.DataFrame(data={'s1':1}, index=['gene1','gene2'])
        sample2_df = pd.DataFrame(data={'s2':1}, index=['gene1','gene2'])
        counts_df = pd.DataFrame(data={'s1':1, 's2':2}, index=['gene1','gene2'])
        stats_df = pd.DataFrame(data={'s1':1, 's2':2}, index=[0])
        self.assertRaisesRegexp(ValueError,
                                'Expected.*but found',
                                deseq2_htseq_merge._validate_dataframes,
                                [sample1_df, sample2_df],
                                counts_df,
                                stats_df)

    def test_validate_dataframes_wrongStatsColumns(self):
        sample1_df = pd.DataFrame(data={'s1':1}, index=['gene1','gene2'])
        sample2_df = pd.DataFrame(data={'s2':1}, index=['gene1','gene2'])
        counts_df = pd.DataFrame(data={'s1':1, 's2':2}, index=['gene1','gene2'])
        stats_df = pd.DataFrame(data={'s1':1, 's2':2, 's3': 3}, index=range(7))
        self.assertRaisesRegexp(ValueError,
                                'Expected columns.*to match',
                                deseq2_htseq_merge._validate_dataframes,
                                [sample1_df, sample2_df],
                                counts_df,
                                stats_df)

    def test_save_merged_dataframes(self):
        def mock_log(message, *args):
            pass
        df_contents = StringIO(\
'''GeneId|s1|s2|s3
gene1|1|3|4''')
        counts_df =  pd.read_csv(df_contents, sep='|', index_col=0)

        df_contents = StringIO(\
'''readcounts|s1|s2|s3
feature_assigned|1|3|4
total|1|3|4''')
        stats_df =  pd.read_csv(df_contents, sep='|', index_col=0)

        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            args = Namespace(counts_filename = join(temp_dir_path, 'counts.txt'),
                             stats_filename = join(temp_dir_path, 'stats.txt'))

            deseq2_htseq_merge._save_merged_dataframes(args,
                                                       counts_df,
                                                       stats_df,
                                                       mock_log)

            with open(args.counts_filename, 'r') as actual_file:
                counts_lines = actual_file.readlines()
            with open(args.stats_filename, 'r') as actual_file:
                stats_lines = actual_file.readlines()

        line = iter(counts_lines)
        self.assertEqual('GeneId\ts1\ts2\ts3\n', next(line))
        self.assertEqual('gene1\t1\t3\t4\n', next(line))
        self.assertRaises(StopIteration, next, line)

        line = iter(stats_lines)
        self.assertEqual('readcounts\ts1\ts2\ts3\n', next(line))
        self.assertEqual('feature_assigned\t1\t3\t4\n', next(line))
        self.assertEqual('total\t1\t3\t4\n', next(line))
        self.assertRaises(StopIteration, next, line)

    def test_build_sample_dataframes(self):
        def mock_log(message, *args):
            pass
        count_s1 = \
'''g1|1
g2|2
g3|3'''.replace('|', '\t')
        count_s2 = \
'''g1|10
g2|20'''.replace('|', '\t')
        with TempDirectory() as temp_dir:
            temp_dir.write('sample1_suffix.txt', count_s1.encode())
            temp_dir.write('sample2_suffix.txt', count_s2.encode())
            temp_dir.write('ignored.txt', 'foo'.encode())
            args = Namespace(htseq_dir = temp_dir.path,
                             suffix = '_suffix.txt')

            actual_dfs = deseq2_htseq_merge._build_sample_dataframes(args, mock_log)

        self.assertEqual((3,1), actual_dfs[0].shape)
        self.assertEqual(['sample1'], list(actual_dfs[0].columns))
        self.assertEqual('GeneId', actual_dfs[0].index.name)
        self.assertEqual((2,1), actual_dfs[1].shape)
        self.assertEqual(['sample2'], list(actual_dfs[1].columns))
        self.assertEqual('GeneId', actual_dfs[1].index.name)
        self.assertEqual(2, len(actual_dfs))

    def test_main(self):
        def mock_log(message, *args):
            pass
        count_s1 = \
'''g1|1
g2|2
g3|3
__ambiguous|4'''.replace('|', '\t')
        count_s2 = \
'''g1|10
g2|20
__ambiguous|40'''.replace('|', '\t')
        with TempDirectory() as temp_dir:
            counts_filename = join(temp_dir.path, 'counts.txt')
            stats_filename = join(temp_dir.path, 'stats.txt')
            temp_dir.write('sample1_htseq.txt', count_s1.encode())
            temp_dir.write('sample2_htseq.txt', count_s2.encode())
            args = ['--htseq_dir=' + temp_dir.path,
                    '--suffix=_htseq.txt',
                    '--counts_filename=' + counts_filename,
                    '--stats_filename=' + stats_filename]
            deseq2_htseq_merge.main(args, mock_log)

            with open(counts_filename, 'r') as actual_file:
                counts_lines = actual_file.readlines()
            with open(stats_filename, 'r') as actual_file:
                stats_lines = actual_file.readlines()

        line = iter(counts_lines)
        self.assertEqual('GeneId\tsample1\tsample2\n', next(line))
        self.assertEqual('g1\t1\t10\n', next(line))
        self.assertEqual('g2\t2\t20\n', next(line))
        self.assertEqual('g3\t3\t0\n', next(line))
        self.assertRaises(StopIteration, next, line)

        line = iter(stats_lines)
        self.assertEqual('readcounts\tsample1\tsample2\n', next(line))
        self.assertEqual('__ambiguous\t4\t40\n', next(line))
        self.assertEqual('feature_assigned\t6\t30\n', next(line))
        self.assertEqual('total\t10\t70\n', next(line))
        self.assertRaises(StopIteration, next, line)
