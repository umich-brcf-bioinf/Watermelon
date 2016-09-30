#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import os
import subprocess
import unittest

from testfixtures.tempdirectory import TempDirectory

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

    def test_basecase(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path

            script_name = os.path.join(SCRIPTS_DIR, 'flag_diffex.py')
            linear_fold_change_threshold = 2
            input_filename = os.path.join(temp_dir_path, 'input.txt')
            output_filename = os.path.join(temp_dir_path, 'output.txt')

            input_file_contents = \
'''field1|log2(fold_change)
value-2|-2
value-1|-1
value0|0
value1|1
value2|2'''.replace('|', '\t')

            with open(input_filename, 'w') as input_file:
                input_file.write(input_file_contents)

            command = '{} --foldchange {} {} {}'.format(script_name,
                                        linear_fold_change_threshold,
                                        input_filename,
                                        output_filename)
            exit_code, command_output = self.execute(command)

            self.assertEqual(0, exit_code)
            with open(output_filename) as output_file:
                actual_output = output_file.readlines()

        self.assertEqual(len(input_file_contents.split('\n')), len(actual_output))

        header_fields = actual_output[0].strip().split('\t')
        expected_fields = ['field1', 'log2(fold_change)', 'significant', 'linear_fold_change', 'diff_exp']
        self.assertEqual(expected_fields, header_fields)

        row_fields = actual_output[1].strip().split('\t')
        self.assertEqual(['value-2', '-2', '1', '0.25', '1'], row_fields)

        row_fields = actual_output[2].strip().split('\t')
        self.assertEqual(['value-1', '-1', '1', '0.5', '1'], row_fields)

        row_fields = actual_output[3].strip().split('\t')
        self.assertEqual(['value0', '0', '1', '1.0', '1'], row_fields)

        row_fields = actual_output[4].strip().split('\t')
        self.assertEqual(['value1', '1', '1', '2.0', '1'], row_fields)

        row_fields = actual_output[5].strip().split('\t')
        self.assertEqual(['value2', '2', '1', '4.0', '1'], row_fields)

        
    def test_failsIfMissingRequiredColumns(self):
        pass

    def test_failsIfFoldchangeNotFloat(self):
        pass

    def test_failsIfLogFoldchangeNegative(self):
        pass
