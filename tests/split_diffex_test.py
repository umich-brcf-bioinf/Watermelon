#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import os
import subprocess
import sys
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd
from testfixtures.tempdirectory import TempDirectory

import scripts.split_diffex as split_diffex

class MockHandler(object):
    def __init__(self):
        self._comparisons = {}

    def handle(self, comparison_name, comparison_df):
        self._comparisons[comparison_name] = comparison_df

    def end(self):
        pass

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')
class SplitDiffexTest(unittest.TestCase):
    def test_raisesValueErrorIfMissingRequiredColumns(self):
        args = Namespace(input_filepath='input.txt')
        df_contents = StringIO('field1\n')
        df = pd.read_csv(df_contents)
        self.assertRaisesRegexp(ValueError,
                                (r'Input file \[input.txt\] is missing required '
                                 r'field\(s\) \[diff_exp,log2\(fold_change\),sample_1,'
                                 r'sample_2,significant,status\].'),
                                 split_diffex._validate_required_fields,
                                 df,
                                 args)

    def test_sort(self):
        df_contents = StringIO(\
'''expected_row|status|significant|diff_exp|log2(fold_change)
19|FAIL|no|no|1
18|FAIL|no|no|2
17|FAIL|no|yes|1
16|FAIL|no|yes|2
15|FAIL|yes|no|1
14|FAIL|yes|no|2
13|FAIL|yes|yes|1
12|FAIL|yes|yes|2
11|OK|no|no|1
10|OK|no|no|2
9|OK|no|yes|1
8|OK|no|yes|2
7|OK|yes|no|-2
6|OK|yes|no|-1
5|OK|yes|no|1
4|OK|yes|no|2
3|OK|yes|yes|-Inf
2|OK|yes|yes|-2
1|OK|yes|yes|-1
0|OK|yes|yes|Inf''')
        df = pd.read_csv(df_contents, sep='|')
        actual_df = split_diffex._sort(df)

        expected_row = 0
        for index, row in actual_df.iterrows():
            self.assertEqual(expected_row, row['expected_row'])
            expected_row += 1

    def test_get_sort_value(self):
        def row(status, significant, diff_exp, log2_fold_change):
            df = pd.DataFrame.from_dict([{'status': status,
                                          'significant' : significant,
                                          'diff_exp' : diff_exp,
                                          'log2(fold_change)' : log2_fold_change}])
            return df.iloc[0]

        def assertLessThan(rowA, rowB):
            sortA = split_diffex._get_sort_value(rowA)
            sortB = split_diffex._get_sort_value(rowB)
            self.assertTrue(sortA < sortB, str(sortA) + ' not < ' + str(sortB))

        base = row('OK', 'yes', 'yes', 2)
        equals = row('OK', 'yes', 'yes', 2)
        lesser_fold_change = row('OK', 'yes', 'yes', 1)
        lesser_diff_exp = row('OK', 'yes', 'no', 2)
        lesser_significant = row('OK', 'no', 'yes', 2)
        lesser_status = row('FAIL', 'yes', 'yes', 2)

        self.assertEqual(split_diffex._get_sort_value(base),
                         split_diffex._get_sort_value(equals))
        assertLessThan(base, lesser_fold_change)
        assertLessThan(base, lesser_diff_exp)
        assertLessThan(base, lesser_significant)
        assertLessThan(base, lesser_status)

    def test_split_comparisons_callsHandleForEachDistinctSetOfValues(self):
        df_contents = StringIO(\
'''sample_1|sample_2|foo|bar
E|F|1|2
A|B|1|2
C|D|1|2
E|F|1|2
E|F|1|2
C|D|1|2
E|F|1|2''')
        df = pd.read_csv(df_contents, sep='|')
    
        mock_handler = MockHandler()
        mock_logger = lambda x: None
        group_by_cols = ['sample_1', 'sample_2']
        split_diffex._split_comparisons('_', df, group_by_cols, mock_handler, mock_logger)

        self.assertEqual(3, len(mock_handler._comparisons))
        self.assertEqual((1,4), mock_handler._comparisons['A_B'].shape)
        self.assertEqual((2,4), mock_handler._comparisons['C_D'].shape)
        self.assertEqual((4,4), mock_handler._comparisons['E_F'].shape)

    def test_validate_included_comparisons_present_raisesExceptionIfRequestedComparisonMissing(self ):
        args = Namespace(input_filepath='input.txt',
                         included_comparisons='A_B,E_F,C_D',
                         comparison_infix='_')
        df_contents = StringIO(\
'''sample_1|sample_2
A|B''')
        df = pd.read_csv(df_contents, sep='|')
        self.assertRaisesRegexp(ValueError,
                                (r'Input file \[input.txt\] is missing requested '
                                 r'comparison\(s\) \[C_D,E_F\].'),
                                 split_diffex._validate_included_comparisons_present,
                                 df,
                                 args)

class FilteringHandlerTest(unittest.TestCase):
    def test_handle_passthroughIncludedComparisons(self):
        df = 'data_frame'
        base_handler = MockHandler()
        mock_log = []
        handler = split_diffex._FilteringHandler(['A_B'], base_handler, mock_log.append)

        handler.handle('A_B', df)

        self.assertEqual(1, len(base_handler._comparisons))
        self.assertEqual(df, base_handler._comparisons['A_B'])
        self.assertEqual(1, len(mock_log))
        self.assertEqual('split comparison [A_B]', mock_log[0])

    def test_handle_ignoreNonIncludedComparisons(self):
        df = pd.DataFrame()
        base_handler = MockHandler()
        mock_log = []
        handler = split_diffex._FilteringHandler(['A_B'], base_handler, mock_log.append)

        handler.handle('C_D', df)

        self.assertEqual(0, len(base_handler._comparisons))
        self.assertEqual(1, len(mock_log))
        self.assertEqual('skipped comparison [C_D]', mock_log[0])


class ComparisonHandlerTest(unittest.TestCase):
    def test_handle_createsNewFile(self):
        df_contents = StringIO(\
'''sample_1|sample_2|foo|bar
A|B|1|2
A|B|3|4''')
        df = pd.read_csv(df_contents, sep='|')

        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            group_name = 'foo'
            suffix = '.suffix.txt'
            handler = split_diffex._ComparisonHandler(temp_dir_path, suffix)
            handler.handle(group_name, df)

            expected_filename = os.path.join(temp_dir_path, 'foo.suffix.txt')
            with open(expected_filename, 'r') as input_file:
                actual_file_data=input_file.readlines()

        self.assertEqual(3, len(actual_file_data))
        self.assertEqual('sample_1\tsample_2\tfoo\tbar\n', actual_file_data[0])
        self.assertEqual('A\tB\t1\t2\n', actual_file_data[1])
        self.assertEqual('A\tB\t3\t4\n', actual_file_data[2])

class SplitDiffexFunctoinalTest(unittest.TestCase):
    def execute(self, command):
        exit_code = 0
        try:
            actual_output = subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as e:
            exit_code = e.returncode
            actual_output = e.output
        return exit_code, str(actual_output)

    def test_commandReturnsCorrectRowAndColumnCount(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            script_name = os.path.join(SCRIPTS_DIR, 'split_diffex.py')
            input_filename = os.path.join(temp_dir_path, 'input.txt')
            output_dir = temp_dir_path

            input_file_contents = \
'''sample_1|sample_2|status|significant|diff_exp|log2(fold_change)
E|F|OK|yes|yes|2
A|B|OK|yes|yes|2
C|D|OK|yes|yes|2
E|F|OK|yes|yes|2
E|F|OK|yes|yes|2
C|D|OK|yes|yes|2
E|F|OK|yes|yes|2'''.replace('|', '\t')

            with open(input_filename, 'w') as input_file:
                input_file.write(input_file_contents)

            suffix_option = '--output_file_suffix=.gene.txt'
            included_comparisons = 'C^D,E^F'
            comparison_infix_option = '--comparison_infix ^'
            redirect_log = '2>/dev/null'
            command = 'python {} {} {} {} {} {} {}'.format(script_name,
                                                        suffix_option,
                                                        comparison_infix_option,
                                                        input_filename,
                                                        output_dir,
                                                        included_comparisons,
                                                        redirect_log)
            exit_code, command_output = self.execute(command)

            self.assertEqual(0, exit_code, command_output)
            actual_files = sorted(os.listdir(output_dir))
            actual_CD_df = pd.read_csv(os.path.join(output_dir, 'C^D.gene.txt'), sep='\t')
            actual_EF_df = pd.read_csv(os.path.join(output_dir, 'E^F.gene.txt'), sep='\t')
        self.assertEquals(['C^D.gene.txt', 'E^F.gene.txt', 'input.txt'], actual_files)
        self.assertEqual((2,6), actual_CD_df.shape)
        self.assertEqual((4,6), actual_EF_df.shape)
