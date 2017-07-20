#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import sys
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd

from testfixtures.tempdirectory import TempDirectory

from scripts import read_stats_bbmap_single

class ReadStatsBbmapSingleTest(unittest.TestCase):

    def setUp(self):
        self.stderr_saved = sys.stderr
        self.stderr = StringIO()
        sys.stderr = self.stderr

    def tearDown(self):
        sys.stderr = self.stderr_saved


    def test_parse_command_line_args(self):
        actual_args = read_stats_bbmap_single._parse_command_line_args('--sample_id sampleA --output_dir out_dir --input_dir in_dir'.split())
        self.assertEqual('sampleA', actual_args.sample_id)
        self.assertEqual('out_dir', actual_args.output_dir)
        self.assertEqual('in_dir', actual_args.input_dir)

    def test_parse_command_line_args_missingRequiredArg(self):
        self.assertRaises(SystemExit,
                          read_stats_bbmap_single._parse_command_line_args,
                          '--output_dir out_dir --input_dir in_dir'.split())
        self.assertRaises(SystemExit,
                          read_stats_bbmap_single._parse_command_line_args,
                          '--sample_id sampleA --input_dir in_dir'.split())
        self.assertRaises(SystemExit,
                          read_stats_bbmap_single._parse_command_line_args,
                          '--sample_id sampleA --output_dir out_dir'.split())

    def test_calc_mean_basecase(self):
        read_data = \
'''Length|Count
1|1
2|2
4|3
8|4'''
        read_df = pd.read_csv(StringIO(read_data), sep='|')
        actual_mean = read_stats_bbmap_single._calc_mean(read_df)

        self.assertEqual(4.9, actual_mean)

    def test_calc_mean_invertedHeaders(self):
        read_data = \
'''Count|Length
1|1
2|2
3|4
4|8'''
        read_df = pd.read_csv(StringIO(read_data), sep='|')
        actual_mean = read_stats_bbmap_single._calc_mean(read_df)

        self.assertEqual(4.9, actual_mean)

    def test_get_mean_stdev_from_ihist_basecase(self):
        insert_data = \
'''
#Mean|214.904
#Median|176
#Mode|151
#STDev|147.733
#PercentOfPairs|99.958
#InsertSize|Count
'''.replace('|', '\t')
        (actual_mean,
         actual_stdev) = read_stats_bbmap_single._get_mean_stdev_from_ihist('file.txt', StringIO(insert_data))

        self.assertEqual(214.904, actual_mean)
        self.assertEqual(147.733, actual_stdev)

    def test_get_mean_stdev_from_ihist_missingOrInvalidMean(self):
        insert_data = \
'''
#Mode|151
#STDev|147.733
#PercentOfPairs|99.958
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#Mean\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))
        insert_data = \
'''
#Mean
#Mode|151
#PercentOfPairs|99.958
#STDev|147.733
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#Mean\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))
        insert_data = \
'''
#Mean|foo
#Mode|151
#PercentOfPairs|99.958
#STDev|147.733
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#Mean\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))


    def test_get_mean_stdev_from_ihist_missingOrInvalidSTDev(self):
        insert_data = \
'''
#Mean|42
#Mode|151
#PercentOfPairs|99.958
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#STDev\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))
        insert_data = \
'''
#Mean|42
#Mode|151
#PercentOfPairs|99.958
#STDev
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#STDev\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))
        insert_data = \
'''
#Mean|42
#Mode|151
#PercentOfPairs|99.958
#STDev|foo
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#STDev\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))


    def test_get_mean_stdev_from_ihist_missingSTDevAndMeanLabels(self):
        insert_data = \
'''
#Mode|151
#PercentOfPairs|99.958
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#Mean, #STDev\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))
        insert_data = \
'''
#Mean
#STDev
#Mode|151
#PercentOfPairs|99.958
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#Mean, #STDev\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))
        insert_data = \
'''
#Mean|foo
#STDev|bar
#Mode|151
#PercentOfPairs|99.958
'''.replace('|', '\t')
        self.assertRaisesRegex(ValueError,
                               r'File \[file\.txt\] is missing label or invalid value for \[#Mean, #STDev\]',
                               read_stats_bbmap_single._get_mean_stdev_from_ihist,
                               'file.txt',
                               StringIO(insert_data))

    def test_main_basecase(self):
        ihist_data = \
'''
#Mean|214.904
#Median|176
#Mode|151
#STDev|147.733
#PercentOfPairs|99.958
#InsertSize|Count
'''.replace('|', '\t')
        lhist_data = \
'''Length|Count
1|1
2|2
4|3
8|4'''.replace('|', '\t')
        with TempDirectory() as in_dir, TempDirectory() as out_dir:
            in_dir.write('sampleA_ihist.txt', ihist_data.encode())
            in_dir.write('sampleA_lhist.txt', lhist_data.encode())
            args = ('--sample_id sampleA '
                    '--output_dir {} '
                    '--input_dir {}').format(out_dir.path, in_dir.path).split()
            read_stats_bbmap_single.main(args)
            actual_lines = out_dir.read('sampleA_read_stats.txt').split(b'\n')
        line_iter = iter(actual_lines)
        self.assertEqual(b'#sample\tinsert_mean\tinsert_std_dev\tread_mean\tinner_mate_dist', next(line_iter))
        self.assertEqual(b'sampleA\t214.904\t147.733\t4.9\t205.104', next(line_iter))
        self.assertEqual(b'', next(line_iter))
        self.assertRaises(StopIteration, next, line_iter)


    def test_main_missingLhist(self):
        with TempDirectory() as temp_dir:
            temp_dir.write('sampleA_ihist.txt', b'stuff')
            args = '--sample_id sampleA --output_dir out_dir --input_dir {}'.format(temp_dir.path).split()
            self.assertRaisesRegex(OSError,
                                   r"No such file or directory.*/sampleA_lhist.txt",
                                   read_stats_bbmap_single.main,
                                   args)

    def test_main_missingIhist(self):
        with TempDirectory() as temp_dir:
            temp_dir.write('sampleA_lhist.txt', b'stuff')
            args = '--sample_id sampleA --output_dir out_dir --input_dir {}'.format(temp_dir.path).split()
            self.assertRaisesRegex(OSError,
                                   r"No such file or directory.*/sampleA_ihist.txt",
                                   read_stats_bbmap_single.main,
                                   args)





    # def test_main_lhistEmpty(self):
    #     self.assertEqual(1, 2)
    #
    # def test_main_lhistInvalidHeader(self):
    #     self.assertEqual(1, 2)
