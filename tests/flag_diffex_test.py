#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import os
import subprocess
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd
from testfixtures.tempdirectory import TempDirectory

import scripts.flag_diffex as flag_diffex

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')
class FlagDiffexTest(unittest.TestCase):

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

            script_name = os.path.join(SCRIPTS_DIR, 'flag_diffex.py')
            linear_fold_change_threshold = 2
            input_filename = os.path.join(temp_dir_path, 'input.txt')
            output_filename = os.path.join(temp_dir_path, 'output.txt')

            input_file_contents = \
'''field1|log2(fold_change)|q_value|status
value-2|1|1|OK
value-1|1|1|OK'''.replace('|', '\t')

            with open(input_filename, 'w') as input_file:
                input_file.write(input_file_contents)


            redirect_output = ' 2>/dev/null '
            command = 'python {} --foldchange {} {} {} {}'.format(script_name,
                                        linear_fold_change_threshold,
                                        input_filename,
                                        output_filename,
                                        redirect_output)
            exit_code, command_output = self.execute(command)

            self.assertEqual(0, exit_code, command_output)
            actual_df = pd.read_csv(output_filename, sep='\t')
        actual_rows, actual_columns = actual_df.shape
        input_row_count = len(input_file_contents.split('\n'))
        self.assertEqual(input_row_count,
                         actual_rows + 1)
        input_columns = len(input_file_contents.split('\n')[0].split('\t'))
        self.assertEqual(input_columns + 2,
                         actual_columns)

    def test_add_columns_addsLinearFoldChange(self):
        df_contents = StringIO(\
'''field1|log2(fold_change)|q_value|status
value-2|-2|1|OK
value-1|-1|1|OK
value0|0|1|OK
value1|1|1|OK
value2|2|1|OK
valueInf|Inf|1|OK
valueNegInf|-Inf|1|OK''')
        df = pd.read_csv(df_contents, sep='|')
        actual_df = flag_diffex._add_columns(df, linear_fold_change=4)

        self.assertIn('linear_fold_change',
                      actual_df.columns.values.tolist())

        log_to_linear = actual_df.set_index('log2(fold_change)')['linear_fold_change'].to_dict()
        self.assertEqual(0.25, log_to_linear[-2])
        self.assertEqual(0.5, log_to_linear[-1])
        self.assertEqual(1.0, log_to_linear[0])
        self.assertEqual(2.0, log_to_linear[1])
        self.assertEqual(4.0, log_to_linear[2])
        pos_inf = float("inf")
        neg_inf = float("-inf")
        self.assertEqual(pos_inf, log_to_linear[pos_inf])
        self.assertEqual(0, log_to_linear[neg_inf])

    def test_add_columns_addsDiffExpYes(self):
        df_contents = StringIO(\
'''test|status|q_value|log2(fold_change)
underExpressed-boundary|OK|0.05|-2
overExpressed-boundary|OK|0.05|2
underExpressed-betterFC|OK|0.05|-3
overExpressed-betterFC|OK|0.05|3
overExpressed-infFC|OK|0.05|Inf
underExpressed-infFC|OK|0.05|-Inf
underExpressed-betterFDR|OK|0.04|-2''')
        df = pd.read_csv(df_contents, sep='|')
        actual_df = flag_diffex._add_columns(df, linear_fold_change=4)
        self.assertIn('diff_exp',
                      actual_df.columns.values.tolist())
        actual_diff_exp = actual_df.set_index('test')['diff_exp'].to_dict()
        self.assertEqual('Yes', actual_diff_exp['underExpressed-boundary'])
        self.assertEqual('Yes', actual_diff_exp['overExpressed-boundary'])
        self.assertEqual('Yes', actual_diff_exp['underExpressed-betterFC'])
        self.assertEqual('Yes', actual_diff_exp['overExpressed-betterFC'])
        self.assertEqual('Yes', actual_diff_exp['underExpressed-infFC'])
        self.assertEqual('Yes', actual_diff_exp['overExpressed-infFC'])
        self.assertEqual('Yes', actual_diff_exp['underExpressed-betterFDR'])

    def test_add_columns_addsDiffExpNo(self):
        df_contents = StringIO(\
'''test|status|q_value|log2(fold_change)
underExpressed-badStatus|NotOK|0.05|-2
overExpressed-badStatus|NotOK|0.05|2
underExpressed-badFDR|OK|0.06|-2
underExpressed-infFDR|OK|Inf|-2
underExpressed-badFC|OK|0.05|-1
overExpressed-badFC|OK|0.05|1''')
        df = pd.read_csv(df_contents, sep='|')
        actual_df = flag_diffex._add_columns(df, linear_fold_change = 4)
        self.assertIn('diff_exp',
                      actual_df.columns.values.tolist())
        actual_diff_exp = actual_df.set_index('test')['diff_exp'].to_dict()
        self.assertEqual('No', actual_diff_exp['underExpressed-badStatus'])
        self.assertEqual('No', actual_diff_exp['overExpressed-badStatus'])
        self.assertEqual('No', actual_diff_exp['underExpressed-badFDR'])
        self.assertEqual('No', actual_diff_exp['underExpressed-infFDR'])
        self.assertEqual('No', actual_diff_exp['overExpressed-badFC'])
        self.assertEqual('No', actual_diff_exp['underExpressed-badFC'])

    def test_raisesValueErrorIfMissingRequiredColumns(self):
        args = Namespace(input_filepath='input.txt', foldchange=1)
        df_contents = StringIO('field1\n')
        df = pd.read_csv(df_contents)
        self.assertRaisesRegexp(ValueError,
                                (r'Input file \[input.txt\] is missing required '
                                 r'field\(s\) \[log2\(fold_change\),q_value,status\].'),
                                 flag_diffex._validate_inputs,
                                 df,
                                 args)

    def test_failsIfFoldchangeNotFloat(self):
        args = Namespace(input_filepath='input.txt', foldchange=1)
        df_contents = StringIO(\
'''field1|status|q_value|log2(fold_change)
A|OK|1|Hello
B|OK|1|5
C|OK|1|World
D|OK|1|10
E|OK|1|Inf''')
        df = pd.read_csv(df_contents, sep='|')
        self.assertRaisesRegexp(ValueError,
                                (r'Input file \[input.txt\]: 2 log2\(fold_change\) '
                                 r'value\(s\) are not numeric: \[Hello:2,World:4\]'),
                                flag_diffex._validate_inputs,
                                df,
                                args)

    def test_failsIfFoldchangeThresholdLessThan1(self):
        args = Namespace(input_file='input.txt', foldchange=0.75)
        df_contents = StringIO(\
'''field1|status|q_value|log2(fold_change)
A|OK|1|5''')
        df = pd.read_csv(df_contents, sep='|')
        self.assertRaisesRegexp(ValueError,
                                r'Specified foldchange cutoff \[0.75\] must be >= 1',
                                flag_diffex._validate_inputs,
                                df,
                                args)

