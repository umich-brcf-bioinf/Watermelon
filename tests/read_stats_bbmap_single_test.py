#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import sys
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd

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
#Mean	214.904
#Median	176
#Mode	151
#STDev	147.733
#PercentOfPairs	99.958
#InsertSize	Count
'''
        (actual_mean, actual_stdev) = read_stats_bbmap_single._get_mean_stdev_from_ihist(StringIO(insert_data))

        self.assertEqual(214.904, actual_mean)
        self.assertEqual(147.733, actual_stdev)


    # def test_main_lhistEmpty(self):
    #     self.assertEqual(1, 2)
    #
    # def test_main_lhistInvalidHeader(self):
    #     self.assertEqual(1, 2)
