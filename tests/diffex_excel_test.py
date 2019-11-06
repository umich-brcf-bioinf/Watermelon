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

import scripts.diffex_excel as diffex_excel

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')
class DiffexExcelTest(unittest.TestCase):
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
            # temp_dir_path = '/nfs/med-bfx-activeprojects/trsaari/sandbox/20191031_diffex_excel_test'

            script_name = os.path.join(SCRIPTS_DIR, 'diffex_excel.py')
            linear_fold_change_threshold = 2
            input_filename = os.path.join(temp_dir_path, 'input.txt')
            glossary_filename = os.path.join(temp_dir_path, 'glossary.txt')
            info_filename = os.path.join(temp_dir_path, 'info.txt')
            output_filename = os.path.join(temp_dir_path, 'output.xlsx')

            input_file_contents = \
'''test_id|external_gene_name|entrezgene_id|description|value_1|value_2|log2(fold_change)|test_stat|p_value
Gm11468|Gm11468|670775|desc1|0|0.358499|inf|-nan|0.0005
Kdm5d|Kdm5d|20592|desc2|0.0212861|4.58465|7.75076|8.99024|5.00E-05
Eif2s3y|Eif2s3y|26908|desc3|0.154875|30.2905|7.6|5.1|5.00E-05'''.replace('|', '\t')
            with open(input_filename, 'w') as input_file:
                input_file.write(input_file_contents)

            glossary_contents = \
'''column_name|description
foo|bar
baz|hoopy
frood|blah'''.replace('|', '\t')
            with open(glossary_filename, 'w') as glossary_file:
                glossary_file.write(glossary_contents)

            info_contents = \
'''foo
foo bar
baz     hoopy
    frood
    blah'''
            with open(info_filename, 'w') as info_file:
                info_file.write(info_contents)

            redirect_output = ' 2>/dev/null '
            command = ('python {} '
                       '-g {} -i {} '
                       '--glossary_filepath {} --info_filepath {} '
                       '{} '
                       '{} ').format(script_name,
                                     input_filename, input_filename,
                                     glossary_filename, info_filename,
                                     output_filename,
                                     redirect_output)
            exit_code, command_output = self.execute(command)

            self.assertEqual(0, exit_code, command_output)
            actual_genes_df = pd.read_excel(output_filename, 'genes')
            actual_isoforms_df = pd.read_excel(output_filename, 'isoforms')
            actual_glossary_df = pd.read_excel(output_filename, 'glossary')
            actual_info_df = pd.read_excel(output_filename, 'info', header=None)

        actual_rows, actual_columns = actual_genes_df.shape
        input_row_count = 3
        self.assertEqual(input_row_count, actual_rows)
        input_columns = 9
        self.assertEqual(input_columns + 1, actual_columns)

        self.assertEqual(actual_genes_df.shape, actual_isoforms_df.shape)

        self.assertEqual((3,2), actual_glossary_df.shape)
        self.assertEqual((5,1), actual_info_df.shape)
