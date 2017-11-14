#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import os
from os.path import join
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import unittest

from testfixtures.tempdirectory import TempDirectory
import pandas as pd

import scripts.combine_summaries as combine_summaries

class CombineSummariesTest(unittest.TestCase):

    def test_main(self):
        mock_logger = lambda x: None
        summaryA_contents = \
'''foo|bar
Af1|Ab1
Af2|Ab2'''.replace('|', '\t')

        summaryB_contents = \
'''foo|bar
Bf1|Bb1
Bf2|Bb2'''.replace('|', '\t')

        with TempDirectory() as temp_dir:
            temp_dir.write('summaryA.txt', summaryA_contents, encoding='utf-8')
            temp_dir.write('summaryB.txt', summaryB_contents, encoding='utf-8')

            output_base = join(temp_dir.path, 'combined')
            combine_summaries.main(['--output_base',
                                    output_base,
                                    join(temp_dir.path, 'summaryA.txt'),
                                    join(temp_dir.path, 'summaryB.txt')],
                                   mock_logger)

            actual_df = pd.read_csv(output_base + '.txt', sep='\t')

        self.assertEqual(4, len(actual_df))
        self.assertEqual(['source','foo','bar'], list(actual_df.columns))
        rows = actual_df.iterrows()
        self.assertEqual(['summaryA', 'Af1', 'Ab1'], next(rows)[1].tolist())
        self.assertEqual(['summaryA', 'Af2', 'Ab2'], next(rows)[1].tolist())
        self.assertEqual(['summaryB', 'Bf1', 'Bb1'], next(rows)[1].tolist())
        self.assertEqual(['summaryB', 'Bf2', 'Bb2'], next(rows)[1].tolist())
        self.assertRaises(StopIteration, next, rows)


    def test_main_trivial(self):
        mock_logger = lambda x: None
        summaryA_contents = \
'''foo|bar
Af1|Ab1
Af2|Ab2'''.replace('|', '\t')

        with TempDirectory() as temp_dir:
            temp_dir.write('summaryA.txt', summaryA_contents, encoding='utf-8')

            output_base = join(temp_dir.path, 'combined')
            combine_summaries.main(['--output_base',
                                    output_base,
                                    join(temp_dir.path, 'summaryA.txt')],
                                   mock_logger)

            actual_df = pd.read_csv(output_base + '.txt', sep='\t')

        self.assertEqual(2, len(actual_df))
        self.assertEqual(['source','foo','bar'], list(actual_df.columns))
        rows = actual_df.iterrows()
        self.assertEqual(['summaryA', 'Af1', 'Ab1'], next(rows)[1].tolist())
        self.assertEqual(['summaryA', 'Af2', 'Ab2'], next(rows)[1].tolist())
        self.assertRaises(StopIteration, next, rows)

    def test_main_inconsistentColumns(self):
        mock_logger = lambda x: None
        summaryA_contents = \
'''foo|bar
Af1|Ab1
Af2|Ab2'''.replace('|', '\t')

        summaryB_contents = \
'''foo|baz
Bf1|Bb1
Bf2|Bb2'''.replace('|', '\t')

        with TempDirectory() as temp_dir:
            temp_dir.write('summaryA.txt', summaryA_contents, encoding='utf-8')
            temp_dir.write('summaryB.txt', summaryB_contents, encoding='utf-8')

            output_base = join(temp_dir.path, 'combined')
            combine_summaries.main(['--output_base',
                                    output_base,
                                    join(temp_dir.path, 'summaryA.txt'),
                                    join(temp_dir.path, 'summaryB.txt')],
                                   mock_logger)

            actual_df = pd.read_table(output_base + '.txt', dtype=str, keep_default_na=False)

        self.assertEqual(4, len(actual_df))
        self.assertEqual(['source','bar', 'baz', 'foo'], list(actual_df.columns))
        rows = actual_df.iterrows()
        self.assertEqual(['summaryA', 'Ab1', '', 'Af1'], next(rows)[1].tolist())
        self.assertEqual(['summaryA', 'Ab2', '', 'Af2'], next(rows)[1].tolist())
        self.assertEqual(['summaryB', '', 'Bb1', 'Bf1'], next(rows)[1].tolist())
        self.assertEqual(['summaryB', '', 'Bb2', 'Bf2'], next(rows)[1].tolist())
        self.assertRaises(StopIteration, next, rows)

    def test_main_inconsistentColumnsLogsWarning(self):
        log_messages = []
        mock_logger = lambda x: log_messages.append(x)

        summaryA_contents = \
'''foo|bar
Af1|Ab1
Af2|Ab2'''.replace('|', '\t')

        summaryB_contents = \
'''foo|baz
Bf1|Bb1
Bf2|Bb2'''.replace('|', '\t')

        with TempDirectory() as temp_dir:
            temp_dir.write('summaryA.txt', summaryA_contents, encoding='utf-8')
            temp_dir.write('summaryB.txt', summaryB_contents, encoding='utf-8')

            output_base = join(temp_dir.path, 'combined')
            combine_summaries.main(['--output_base',
                                    output_base,
                                    join(temp_dir.path, 'summaryA.txt'),
                                    join(temp_dir.path, 'summaryB.txt')],
                                   mock_logger)

        warnings = [x for x in log_messages if x.startswith('WARNING')]
        self.assertEqual(['WARNING: incoming headers were not consistent'], warnings)

    def test_main_emptyFileWarning(self):
        log_messages = []
        mock_logger = lambda x: log_messages.append(x)

        summaryA_contents = \
'''foo|bar
Af1|Ab1
Af2|Ab2'''.replace('|', '\t')

        summaryB_contents = \
'''foo|bar
'''.replace('|', '\t')

        with TempDirectory() as temp_dir:
            temp_dir.write('summaryA.txt', summaryA_contents, encoding='utf-8')
            temp_dir.write('summaryB.txt', summaryB_contents, encoding='utf-8')

            output_base = join(temp_dir.path, 'combined')
            empty_file_name = join(temp_dir.path, 'summaryB.txt')
            combine_summaries.main(['--output_base',
                                    output_base,
                                    join(temp_dir.path, 'summaryA.txt'),
                                    empty_file_name],
                                   mock_logger)

        warnings = [x for x in log_messages if x.startswith('WARNING')]
        self.assertEqual(['WARNING: file [{}] is empty'.format(empty_file_name)],
                         warnings)


    def test_main_missingHeaderError(self):
        mock_logger = lambda x: None

        summaryA_contents = \
'''foo|bar
Af1|Ab1
Af2|Ab2'''.replace('|', '\t')

        summaryB_contents = ''

        with TempDirectory() as temp_dir:
            temp_dir.write('summaryA.txt', summaryA_contents, encoding='utf-8')
            temp_dir.write('summaryB.txt', summaryB_contents, encoding='utf-8')

            output_base = join(temp_dir.path, 'combined')
            self.assertRaises(Exception,
                            combine_summaries.main,
                            ['--output_base',
                             output_base,
                             join(temp_dir.path, 'summaryA.txt'),
                             join(temp_dir.path, 'summaryB.txt')],
                            mock_logger)
