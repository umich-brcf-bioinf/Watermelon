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

class MockGroupHandler(object):
    def __init__(self):
        self._groups = {}

    def handle(self, group_name, group_df):
        self._groups[group_name] = group_df

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')
class SplitDiffexTest(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.stderr = StringIO()
        self._saved_stderr = sys.stderr
        sys.stderr = self.stderr

    def tearDown(self):
        self.stderr.close()
        sys.stderr = self._saved_stderr
        unittest.TestCase.tearDown(self)

    def test_raisesValueErrorIfMissingRequiredColumns(self):
        args = Namespace(input_file='input.txt', foldchange=1)
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

    def test_split_groups_by_sample_callsHandleForEachDistinctSetOfValues(self):
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
    
        mock_handler = MockGroupHandler()
        mock_logger = lambda x: None
        split_diffex._split_groups_by_sample(df, mock_handler.handle, mock_logger)

        self.assertEqual(3, len(mock_handler._groups))
        self.assertEqual((1,4), mock_handler._groups['A_B'].shape)
        self.assertEqual((2,4), mock_handler._groups['C_D'].shape)
        self.assertEqual((4,4), mock_handler._groups['E_F'].shape)

    def test_handler_createsNewFileForEachGroup(self):
        df_contents = StringIO(\
'''sample_1|sample_2|foo|bar
A|B|1|2
A|B|3|4''')
        df = pd.read_csv(df_contents, sep='|')

        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            group_name = 'foo'
            suffix = '.csv'
            handler = split_diffex._GroupHandler(temp_dir_path, suffix)
            handler.handle(group_name, df)

            expected_filename = os.path.join(temp_dir_path, 'foo.csv')
            with open(expected_filename, 'r') as input_file:
                actual_file_data=input_file.readlines()

        self.assertEqual(3, len(actual_file_data))
        self.assertEqual('sample_1\tsample_2\tfoo\tbar\n', actual_file_data[0])
        self.assertEqual('A\tB\t1\t2\n', actual_file_data[1])
        self.assertEqual('A\tB\t3\t4\n', actual_file_data[2])

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

            command = 'python {} {} {}'.format(script_name,
                                               input_filename,
                                               output_dir)
            exit_code, command_output = self.execute(command)

            self.assertEqual(0, exit_code, command_output)
            actual_AB_df = pd.read_csv(os.path.join(output_dir, 'A_B.csv'), sep='\t')
            actual_CD_df = pd.read_csv(os.path.join(output_dir, 'C_D.csv'), sep='\t')
            actual_EF_df = pd.read_csv(os.path.join(output_dir, 'E_F.csv'), sep='\t')
        self.assertEqual((1,6), actual_AB_df.shape)
        self.assertEqual((2,6), actual_CD_df.shape)
        self.assertEqual((4,6), actual_EF_df.shape)
