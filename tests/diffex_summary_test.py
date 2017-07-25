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
from testfixtures.tempdirectory import tempdir

import scripts.diffex_summary as diffex_summary

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')
class DiffexSummarizeTest(unittest.TestCase):
    def test_parse_command_line_args(self):
        command_line = (' --output_file=OUTPUT_FILE '
                        ' --output_xlsx=OUTPUT_XLSX '
                        ' --annotation_column=ANNOT_COLUMN '
                        ' --annotation_null=ANNOT_NULL '
                        ' --diffex_call_column=DE_COLUMN '
                        ' --diffex_call_pass=CALL_PASS '
                        ' --trim_suffix=TRIM_SUFFIX '
                        ' /a/foo.baz bar.baz c/hoopy')
        args = diffex_summary._parse_command_line_args(command_line.split())

        self.assertEqual('OUTPUT_FILE', args.output_file)
        self.assertEqual('OUTPUT_XLSX', args.output_xlsx)
        self.assertEqual('ANNOT_COLUMN', args.annotation_column)
        self.assertEqual('ANNOT_NULL', args.annotation_null)
        self.assertEqual('DE_COLUMN', args.diffex_call_column)
        self.assertEqual('CALL_PASS', args.diffex_call_pass)
        self.assertEqual('TRIM_SUFFIX', args.trim_suffix)
        self.assertEqual(['/a/foo.baz','bar.baz', 'c/hoopy'], args.input_files)

    def test_parse_command_line_args_defaults(self):
        command_line = (' --output_file=OUTPUT_FILE '
                        ' --diffex_call_column=DE_COLUMN '
                        ' --diffex_call_pass=CALL_PASS '
                        ' /a/foo.baz bar.baz c/hoopy')
        args = diffex_summary._parse_command_line_args(command_line.split())

        self.assertEqual('gene_id', args.annotation_column)
        self.assertEqual('.', args.annotation_null)
        self.assertEqual('.annot.txt', args.trim_suffix)

    def test_simplify_file_names(self):
        file_paths = ['a/b/1.xlsx','a/b/2.xlsx','a/b/3']
        suffix = '.xlsx'
        actual_names = diffex_summary._simplify_file_names(file_paths, suffix)

        expected_names = {'a/b/1.xlsx':'1',
                          'a/b/2.xlsx': '2',
                          'a/b/3' : '3'}
        self.assertEqual(expected_names, actual_names)

    def test_simplify_file_names_noSuffix(self):
        file_paths = ['a/b/1.xlsx','a/b/2.xlsx','a/b/3']
        suffix = None
        actual_names = diffex_summary._simplify_file_names(file_paths, suffix)

        expected_names = {'a/b/1.xlsx':'1.xlsx',
                          'a/b/2.xlsx': '2.xlsx',
                          'a/b/3' : '3'}
        self.assertEqual(expected_names, actual_names)

    def test_simplify_file_names_noCommonPrefix(self):
        file_paths = ['a/b/1','a/b/2','c/b/1']
        suffix = None
        actual_names = diffex_summary._simplify_file_names(file_paths, suffix)

        expected_names = {'a/b/1':'a/b/1',
                          'a/b/2': 'a/b/2',
                          'c/b/1' : 'c/b/1'}
        self.assertEqual(expected_names, actual_names)

    def test_simplify_file_names_simplePair(self):
        file_paths = ['a/b/1_gene.xlsx','a/b/1_isoform.xlsx']
        suffix = None
        actual_names = diffex_summary._simplify_file_names(file_paths, suffix)

        expected_names = {'a/b/1_gene.xlsx':'1_gene.xlsx',
                          'a/b/1_isoform.xlsx':'1_isoform.xlsx'}
        self.assertEqual(expected_names, actual_names)

    def test_count_annotated(self):
        df_contents = StringIO(\
'''a|b|gene_id|c
0|0|.|0
0|0|x|0
0|0|.|0
0|0|y|0
0|0|z|0''')
        df = pd.read_csv(df_contents, sep='|')
        annotation_column='gene_id'
        annotation_null = '.'
        actual_annotated = diffex_summary._count_annotated(annotation_column,
                                                             annotation_null,
                                                             df)
        self.assertEqual(3, actual_annotated)

    def test_count_passed(self):
        df_contents = StringIO(\
'''a|b|gene_id|call
0|0|.|YES
0|0|.|0
0|0|.|YES
0|0|.|0
0|0|.|YES''')
        df = pd.read_csv(df_contents, sep='|')
        diffex_call_column='call'
        diffex_call_pass = 'YES'
        actual_pass = diffex_summary._count_diffex_pass(diffex_call_column,
                                                          diffex_call_pass,
                                                          df)
        self.assertEqual(3, actual_pass)

    def test_get_stats(self):
        df_contents = StringIO(\
'''a|b|gene_id|call
0|0|1|YES
0|0|2|0
0|0|.|YES
0|0|3|0
0|0|4|YES''')
        df = pd.read_csv(df_contents, sep='|')
        name = 'myName'
        args = Namespace(annotation_column='gene_id',
                         annotation_null='.',
                         diffex_call_column='call',
                         diffex_call_pass='YES')
        actual_stats = diffex_summary._get_stats(args, df, name)

        expected_stats = {'comparison': 'myName',
                          'total_count': 5,
                          'count_diff_expressed': 3,
                          'count_annotated': 4,
                          'percent_annotated': 80}
        self.assertEqual(expected_stats, actual_stats)

    @tempdir()
    def test_main(self, temp_dir):
        temp_dir_path = temp_dir.path
        df_contents = \
'''a|b|gene_id|call
0|0|1|YES
0|0|2|0
0|0|.|YES
0|0|3|0
0|0|4|YES'''.replace('|', '\t')
        filename_a = os.path.join(temp_dir_path,'a.annot.txt')
        filename_b = os.path.join(temp_dir_path,'b.annot.txt')
        filename_c = os.path.join(temp_dir_path,'c.annot.txt')
        with open(filename_a, 'w') as f:
            print(df_contents, file=f)
        with open(filename_b, 'w') as f:
            print(df_contents, file=f)
        with open(filename_c, 'w') as f:
            print(df_contents, file=f)
        output_filename = os.path.join(temp_dir_path,'summary.txt')
        output_xlsx = os.path.join(temp_dir_path,'summary.xlsx')
        command_line = (' --output_file={output_filename} '
                        ' --output_xlsx={output_xlsx} '
                        ' --annotation_column=gene_id '
                        ' --annotation_null=. '
                        ' --diffex_call_column=call '
                        ' --diffex_call_pass=YES '
                        ' --trim_suffix=.annot.txt '
                        ' {filename_c} {filename_a} {filename_b}'
                        ).format(output_filename=output_filename,
                                 output_xlsx=output_xlsx,
                                 filename_a=filename_a,
                                 filename_b=filename_b,
                                 filename_c=filename_c)
        mock_logger = lambda x: None
        diffex_summary.main(command_line.split(), log=mock_logger)

        actual_df = pd.read_table(output_filename, sep='\t', header=0)

        self.assertEqual(['comparison',
                          'total_count',
                          'count_diff_expressed',
                          'count_annotated',
                          'percent_annotated'],
                         list(actual_df.columns))
        self.assertEqual(3, len(actual_df))
        self.assertEqual(['a','b','c'], list(actual_df['comparison']))
        self.assertEqual(9, sum(actual_df['count_diff_expressed']))
        self.assertEqual(12, sum(actual_df['count_annotated']))

        actual_df = pd.read_excel(output_xlsx)

        self.assertEqual(['comparison',
                          'total_count',
                          'count_diff_expressed',
                          'count_annotated',
                          'percent_annotated'],
                         list(actual_df.columns))
        self.assertEqual(3, len(actual_df))
        self.assertEqual(['a','b','c'], list(actual_df['comparison']))
        self.assertEqual(9, sum(actual_df['count_diff_expressed']))
        self.assertEqual(12, sum(actual_df['count_annotated']))
