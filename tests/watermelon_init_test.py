#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import functools
import os
import os.path
import sys
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures.tempdirectory import TempDirectory

import scripts.watermelon_init as watermelon_init

# TEST_DIR = os.path.realpath(os.path.dirname(__file__))
# SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')

class WatermelonInitTest(unittest.TestCase):
    def test_build_postlude(self):
        args = Namespace(inputs_dir='INPUTS_DIR',
                         source_fastq_dir='SOURCE_FASTQ_DIR',
                         analysis_dir='ANALYSIS_DIR',
                         deliverables_dir='DELIVERABLES_DIR',
                         configfile='CONFIGFILE',
                         job_suffix='_JOB_SUFFIX')
        sample_count = 42
        file_count = 168
        actual_postlude = watermelon_init._build_postlude(args, sample_count, file_count)
        self.assertRegexpMatches(actual_postlude,
                                 r'INPUTS_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'SOURCE_FASTQ_DIR \| 42 samples \| 168 files')
        self.assertRegexpMatches(actual_postlude,
                                 r'ANALYSIS_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'DELIVERABLES_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'screen -S watermelon_JOB_SUFFIX')
        self.assertRegexpMatches(actual_postlude,
                                 r'watermelon --dry-run CONFIGFILE')

    def test_parse_args(self):
        command_line_args = ('--genome_build mm10 '
                             '--job_suffix _11_01_A '
                             '--working_dir /tmp.foo '
                             'DNASeqCore/Run_1286/rhim/Run_1286').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        expected_args = ['analysis_dir',
                         'configfile',
                         'deliverables_dir',
                         'genome_build',
                         'inputs_dir',
                         'job_suffix',
                         'source_fastq_dir',
                         'working_dir',]
        self.assertEqual(expected_args, sorted(vars(args)))
        self.assertEqual('mm10', args.genome_build)
        self.assertEqual('_11_01_A', args.job_suffix)
        self.assertEqual('DNASeqCore/Run_1286/rhim/Run_1286', args.source_fastq_dir)
        join = os.path.join
        basedir = '/tmp.foo'
        expected_analysis_filepath = join(basedir, 'analysis_11_01_A')
        self.assertEqual(expected_analysis_filepath, args.analysis_dir)
        self.assertEqual(os.path.join(expected_analysis_filepath, 'config_11_01_A.yaml'),
                         args.configfile)
        expected_deliverables_filepath = os.path.join(basedir, 'deliverables_11_01_A')
        self.assertEqual(expected_deliverables_filepath, args.deliverables_dir)
        self.assertEqual(os.path.join(basedir, 'inputs', '00-multiplexed_reads'),
                         args.inputs_dir)

    def test_parse_args_workingDirDefaultsToCwd(self):
        command_line_args = ('--genome_build mm10 '
                             '--job_suffix _11_01_A '
                             'DNASeqCore/Run_1286/rhim/Run_1286').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        self.assertEqual(os.getcwd(), args.working_dir) 

    def test_make_top_level_dirs(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            real = functools.partial(os.path.join, temp_dir_path)
            args = Namespace(analysis_dir=real('ANALYSIS_11_10_A'),
                             deliverables_dir=real('DELIVERABLES_11_10_A'),
                             inputs_dir=real('INPUTS'))
            watermelon_init._make_top_level_dirs(args)

            def is_dir(o):
                return os.path.isdir(os.path.join(temp_dir_path,o))
            actual_dirs = sorted([o for o in os.listdir(temp_dir_path) if is_dir(o)])
            self.assertEqual(['ANALYSIS_11_10_A','DELIVERABLES_11_10_A','INPUTS'],
                             actual_dirs)

# class WatermelonInitFunctoinalTest(unittest.TestCase):
#     def execute(self, command):
#         exit_code = 0
#         try:
#             actual_output = subprocess.check_output(command, shell=True)
#         except subprocess.CalledProcessError as e:
#             exit_code = e.returncode
#             actual_output = e.output
#         return exit_code, str(actual_output)
# 
#     def test_commandReturnsCorrectRowAndColumnCount(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             os.
#             script_name = os.path.join(SCRIPTS_DIR, 'watermelon_init.py')
# 
#             redirect_output = '2>/dev/null'
#             command = 'python {} {}'.format(script_name, redirect_output)
#             exit_code, command_output = self.execute(command)
# 
#         self.assertEqual(0, exit_code, command)
#         self.assertEqual('foo', command_output)
