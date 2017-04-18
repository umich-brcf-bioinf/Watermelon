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

import scripts.tuxedo_flip as tuxedo_flip

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')
class TuxedoFlipTest(unittest.TestCase):
    def test_raisesValueErrorIfMissingRequiredColumns(self):
        args = Namespace(input_filepath='input.txt')
        df_contents = StringIO('field1\n')
        df = pd.read_csv(df_contents)
        self.assertRaisesRegexp(ValueError,
                                (r'Input file \[input.txt\] is missing required '
                                 r'field\(s\) \[log2\(fold_change\),sample_1,'
                                 r'sample_2,test_stat,value_1,value_2\].'),
                                 tuxedo_flip._validate_required_fields,
                                 df,
                                 args)

    def test_flip_comparisons_samplesValuesFlippedAndFoldChangeAndTestStatInverted(self):
        df_contents = StringIO(\
'''test_id|sample_1|sample_2|status|value_1|value_2|log2(fold_change)|test_stat
gene0|A|B|foo|4|1|2|3.5
gene1|A|B|foo|3|1|1.584962501|3.5
gene2|C|D|foo|1|3|-1.584962501|-3.5
gene3|C|D|foo|1|4|-2|-3.5''')
        df = pd.read_csv(df_contents, sep='|')
    
        comparison_infix = '^'
        tuxedo_flip._flip_comparisons(comparison_infix, df, ['B^A', 'D^C'])

        self.assertEquals(['gene0', 'B', 'A', 'foo', 1, 4, -2, -3.5], list(df.loc[0].values))
        self.assertEquals(['gene1', 'B', 'A', 'foo', 1, 3, -1.584962501, -3.5], list(df.loc[1].values))
        self.assertEquals(['gene2', 'D', 'C', 'foo', 3, 1, 1.584962501, 3.5], list(df.loc[2].values))
        self.assertEquals(['gene3', 'D', 'C', 'foo', 4, 1, 2, 3.5], list(df.loc[3].values))

    def test_flip_comparisons_noChangeForInputValuesNotInComparisons(self):
        df_contents = StringIO(\
'''test_id|sample_1|sample_2|status|value_1|value_2|log2(fold_change)|test_stat
gene0|A|B|foo|4|1|2|1
gene1|A|B|foo|3|1|1.584962501|1
gene2|C|D|foo|1|3|-1.584962501|-1
gene3|C|D|foo|1|4|-2|-1''')
        df = pd.read_csv(df_contents, sep='|')
    
        tuxedo_flip._flip_comparisons('_', df, ['B_A'])

        self.assertEquals(['gene2', 'C', 'D', 'foo', 1, 3, -1.584962501, -1], list(df.loc[2].values))
        self.assertEquals(['gene3', 'C', 'D', 'foo', 1, 4, -2, -1], list(df.loc[3].values))

    def test_flip_comparisons_noChangeForCorrectComparisons(self):
        df_contents = StringIO(\
'''test_id|sample_1|sample_2|status|value_1|value_2|log2(fold_change)|test_stat
gene0|A|B|foo|4|1|2|1
gene1|A|B|foo|3|1|1.584962501|1
gene2|C|D|foo|1|3|-1.584962501|-1
gene3|C|D|foo|1|4|-2|-1''')
        df = pd.read_csv(df_contents, sep='|')
    
        tuxedo_flip._flip_comparisons('_', df, ['B_A,C_D'])

        self.assertEquals(['gene2', 'C', 'D', 'foo', 1, 3, -1.584962501, -1], list(df.loc[2].values))
        self.assertEquals(['gene3', 'C', 'D', 'foo', 1, 4, -2, -1], list(df.loc[3].values))


class FlipDiffexFunctoinalTest(unittest.TestCase):
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
            script_name = os.path.join(SCRIPTS_DIR, 'tuxedo_flip.py')
            input_filename = os.path.join(temp_dir_path, 'input.txt')
            output_filename = os.path.join(temp_dir_path, 'output.txt')

            input_file_contents = \
'''test_id|sample_1|sample_2|status|value_1|value_2|log2(fold_change)|test_stat
gene0|A|B|foo|4|1|2|3.5
gene1|A|B|foo|3|1|1.584962501|3.5
gene2|A|B|foo|2|1|1|3.5
gene3|A|B|foo|1|1|0|0.5
gene4|A|B|foo|1|2|-1|-3.5
gene5|A|B|foo|1|3|-1.584962501|-3.5
gene6|A|B|foo|1|4|-2|-3.5
gene0|C|D|foo|4|1|2|3.5
gene1|C|D|foo|3|1|1.584962501|3.5
gene2|C|D|foo|2|1|1|3.5'''.replace('|', '\t')

            with open(input_filename, 'w') as input_file:
                input_file.write(input_file_contents)

            comparisons = 'B^A'
            comparison_infix_option = '--comparison_infix ^'
            redirect_output = '2>/dev/null'
            command = 'python {} {} {} {} {} {}'.format(script_name,
                                                     input_filename,
                                                     output_filename,
                                                     comparisons,
                                                     comparison_infix_option,
                                                     redirect_output)
            exit_code, command_output = self.execute(command)

            self.assertEqual(0, exit_code, command_output)
            input_df = pd.read_csv(input_filename, sep='\t')
            actual_output_df = pd.read_csv(output_filename, sep='\t')
        self.assertEqual(input_df.shape, actual_output_df.shape)
        self.assertEqual(list(input_df.columns.values),
                         list(actual_output_df.columns.values))
