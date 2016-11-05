#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import os
import shutil
import subprocess
import sys
import unittest

from testfixtures.tempdirectory import TempDirectory

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

TESTS_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_ROOT = os.path.dirname(TESTS_DIR)
BIN_DIR = os.path.join(WATERMELON_ROOT, 'bin')
WATERMELON_EXECUTABLE = os.path.join(BIN_DIR, 'watermelon')
class WatermelonTest(unittest.TestCase):

    def execute(self, command):
        exit_code = 0
        try:
            actual_output = subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as e:
            exit_code = e.returncode
            actual_output = e.output
        return exit_code, str(actual_output)

    def setup_tmp_dir(self, temp_dir):
        temp_dir_path = temp_dir.path  #'/tmp/foo' #
        os.chdir(temp_dir_path)
        return temp_dir_path

    def test_watermelon_basecase(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            TEST_SNAKEFILE = os.path.join(TESTS_DIR, 'test.snakefile')

            command = '{} --snakefile {} {}'.format(WATERMELON_EXECUTABLE,
                                                    TEST_SNAKEFILE,
                                                    CONFIG_FILE)
            exit_code, actual_output = self.execute(command)
            self.assertEqual(0, exit_code)

            with open('sample_A.txt') as output_file:
                actual_output_A = output_file.readlines()
            with open('sample_B.txt') as output_file:
                actual_output_B = output_file.readlines()

        self.assertEqual(['lineA\n'], actual_output_A)
        self.assertEqual(['lineB\n'], actual_output_B)


    def test_watermelon_logsResults(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            TEST_SNAKEFILE = os.path.join(TESTS_DIR, 'test.snakefile')

            command = '{} --snakefile {} {}'.format(WATERMELON_EXECUTABLE,
                                               TEST_SNAKEFILE,
                                               CONFIG_FILE)
            exit_code, actual_output = self.execute(command)

            self.assertEqual(0, exit_code)

            log_dirs = glob.glob("logs/*")
            self.assertEquals(1, len(log_dirs))
            log_dir = log_dirs[0]
            log_file_name = os.path.join(log_dir, 'watermelon.log')
            self.assertTrue(os.path.exists(log_file_name))
            with open(log_file_name) as log_file:
                actual_log = "".join(log_file.readlines())

        self.assertRegexpMatches(actual_log, r'Provided cores: 40')
        self.assertRegexpMatches(actual_log, r'3 of 3 steps \(100%\) done')
        self.assertRegexpMatches(actual_log, r'elapsed seconds: \d+')
        self.assertRegexpMatches(actual_log, r'elapsed time: 0h:0m:\d+s')
        self.assertRegexpMatches(actual_log, r'watermelon done')


    def test_watermelon_skipsLogsOnDryrun(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            TEST_SNAKEFILE = os.path.join(TESTS_DIR, 'test.snakefile')

            command = '{} --snakefile {} {} --dryrun'.format(WATERMELON_EXECUTABLE,
                                                             TEST_SNAKEFILE,
                                                             CONFIG_FILE)
            exit_code, actual_output = self.execute(command)

            self.assertEqual(0, exit_code)

            log_dirs = glob.glob("logs/*")
            self.assertEquals(0, len(log_dirs))


    def test_watermelon_copiesConfigAndSnakefileToLogDir(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE=os.path.join(TESTS_DIR, 'config.yaml')
            TEST_SNAKEFILE=os.path.join(TESTS_DIR, 'test.snakefile')

            command = '{} --snakefile {} {}'.format(WATERMELON_EXECUTABLE,
                                                    TEST_SNAKEFILE,
                                                    CONFIG_FILE)
            exit_code, actual_output = self.execute(command)

            self.assertEqual(0, exit_code)

            log_dirs = glob.glob("logs/*")
            self.assertEquals(1, len(log_dirs))
            log_dir = log_dirs[0]
            run_time = os.path.basename(log_dir)
            config_filename = os.path.join(log_dir, 'config.yaml.' + run_time)
            snakefile_filename = os.path.join(log_dir, 'test.snakefile.' + run_time)
            self.assertTrue(os.path.exists(config_filename), config_filename)
            self.assertTrue(os.path.exists(snakefile_filename), snakefile_filename)


    def test_watermelon_runsDetailedSummary(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            TEST_SNAKEFILE = os.path.join(TESTS_DIR, 'test.snakefile')

            command = '{} --snakefile {} {}'.format(WATERMELON_EXECUTABLE,
                                                    TEST_SNAKEFILE,
                                                    CONFIG_FILE)
            exit_code, actual_output = self.execute(command)

            self.assertEqual(0, exit_code)

            log_dirs = glob.glob("logs/*")
            self.assertEquals(1, len(log_dirs))
            log_dir = log_dirs[0]
            run_time = os.path.basename(log_dir)
            detailed_summary_filename = os.path.join(log_dir,
                                                     'watermelon.detailed_summary')
            self.assertTrue(os.path.exists(detailed_summary_filename),
                            detailed_summary_filename)
            with open(detailed_summary_filename) as detailed_summary_file:
                actual_detailed_summary = "".join(detailed_summary_file.readlines())

        self.assertRegexpMatches(actual_detailed_summary, r'output_file\tdate\trule')


    def test_watermelon_failsWhenInvalidConfig(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            with open('invalid_config.yaml', 'w') as config_file:
                print('an invalid config file', file=config_file)
            TEST_SNAKEFILE = os.path.join(TESTS_DIR, 'test.snakefile')

            command = '{} {} {}'.format(WATERMELON_EXECUTABLE,
                                        TEST_SNAKEFILE,
                                        'invalid_config.yaml')
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)

            log_dirs = glob.glob("logs/*")
            self.assertEquals(1, len(log_dirs))
            log_dir = log_dirs[0]
            log_file_name = os.path.join(log_dir, 'watermelon.log')
            self.assertTrue(os.path.exists(log_file_name))
            with open(log_file_name) as log_file:
                actual_log = "".join(log_file.readlines())

        self.assertRegexpMatches(actual_log, r'ERROR.*contact bfxcore support')


    def test_watermelonFailsGracefullyIfSnakefileNotFound(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE=os.path.join(TESTS_DIR, 'config.yaml')

            command = '{} --snakefile iAmNotAFile {}'.format(WATERMELON_EXECUTABLE,
                                                             CONFIG_FILE)
            exit_code, actual_output = self.execute(command)

            self.assertNotEqual(0, exit_code)
            self.assertRegexpMatches(actual_output, r'ERROR.*\[iAmNotAFile\] cannot be read')


    def test_watermelonFailsGracefullyIfSnakemakeNotInstalled(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            command = ('echo foo > snakemake; chmod u+x snakemake; '
                       'PATH=.:$PATH; '
                       '{}').format(WATERMELON_EXECUTABLE)
            exit_code, actual_output = self.execute(command)
            self.assertNotEqual(0, exit_code)
            self.assertRegexpMatches(actual_output, r'snakemake not found')


    def test_watermelon_showsUsageWhenConfigOmitted(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)

            command = WATERMELON_EXECUTABLE
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)
            self.assertRegexpMatches(actual_output, 'Usage')


    def test_watermelon_showsUsageWhenFileNotFound(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)

            command = '{} i_do_not_exist.yaml'.format(WATERMELON_EXECUTABLE)
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)
            self.assertRegexpMatches(actual_output, 'config file .* cannot be read')


    def test_watermelon_showsUsageWhenConfigIsDir(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_DIR = os.path.join(temp_dir_path, 'config_dir')
            os.mkdir(CONFIG_DIR)

            command = '{} {}'.format(WATERMELON_EXECUTABLE, CONFIG_DIR)
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)
            self.assertRegexpMatches(actual_output, 'config file .* cannot be read')


    def test_watermelon_passthroughExtraArguments(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            TEST_SNAKEFILE = os.path.join(TESTS_DIR, 'test.snakefile')

            snakemake_executable_path = os.path.join(temp_dir_path, 'snakemake')
            with open(snakemake_executable_path, 'w') as snakemake_command:
                snakemake_command_contents = '#!/bin/bash\necho $@'''
                print(snakemake_command_contents, file=snakemake_command)
            command = ('chmod u+x snakemake; '
                       'PATH=.:$PATH; '
                       '{} --snakefile {} 1 2 3 baz froody {}'.format(WATERMELON_EXECUTABLE,
                                                                      TEST_SNAKEFILE,
                                                                      CONFIG_FILE))
            exit_code, actual_output = self.execute(command)
            self.assertEqual(0, exit_code)
            self.assertRegexpMatches(actual_output,
                                     (r'--configfile {} '
                                     r'--snakefile {} '
                                      r'--cores 40 '
                                      r'-T '
                                      r'1 2 3 baz froody'). format(CONFIG_FILE,
                                                                   TEST_SNAKEFILE))

    def test_watermelon_expandsConcatenatedShortOptions(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            TEST_SNAKEFILE = os.path.join(TESTS_DIR, 'test.snakefile')

            snakemake_executable_path = os.path.join(temp_dir_path, 'snakemake')
            with open(snakemake_executable_path, 'w') as snakemake_command:
                snakemake_command_contents = '#!/bin/bash\necho $@'''
                print(snakemake_command_contents, file=snakemake_command)
            command = ('chmod u+x snakemake; '
                       'PATH=.:$PATH; '
                       '{} --snakefile {} -abcd {}'.format(WATERMELON_EXECUTABLE,
                                                                      TEST_SNAKEFILE,
                                                                      CONFIG_FILE))
            exit_code, actual_output = self.execute(command)
            self.assertEqual(0, exit_code)
            self.assertRegexpMatches(actual_output,
                                     (r'--configfile {} '
                                     r'--snakefile {} '
                                      r'--cores 40 '
                                      r'-T '
                                      r'-a -b -c -d'). format(CONFIG_FILE,
                                                              TEST_SNAKEFILE))


    def test_watermelon_defaultArgumentsSet(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            DEFAULT_SNAKEFILE = os.path.join(WATERMELON_ROOT, 'rnaseq.snakefile')
            DEFAULT_CORES = 40

            snakemake_executable_path = os.path.join(temp_dir_path, 'snakemake')
            with open(snakemake_executable_path, 'w') as snakemake_command:
                snakemake_command_contents = '#!/bin/bash\necho $@'''
                print(snakemake_command_contents, file=snakemake_command)
            command = ('chmod u+x snakemake; '
                       'PATH=.:$PATH; '
                       '{} {}'.format(WATERMELON_EXECUTABLE, CONFIG_FILE))
            exit_code, actual_output = self.execute(command)
            self.assertEqual(0, exit_code)
            self.assertRegexpMatches(actual_output,
                                     (r'--snakefile {} '
                                      r'--cores {} '
                                      r'-T'). format(DEFAULT_SNAKEFILE,
                                                     DEFAULT_CORES))

    def test_watermelon_canOverrideCores(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE = os.path.join(TESTS_DIR, 'config.yaml')
            DEFAULT_SNAKEFILE = os.path.join(WATERMELON_ROOT, 'rnaseq.snakefile')
            CUSTOM_CORES = 42

            snakemake_executable_path = os.path.join(temp_dir_path, 'snakemake')
            with open(snakemake_executable_path, 'w') as snakemake_command:
                snakemake_command_contents = '#!/bin/bash\necho $@'''
                print(snakemake_command_contents, file=snakemake_command)
            command = ('chmod u+x snakemake; '
                       'PATH=.:$PATH; '
                       '{} --cores {} {}'.format(WATERMELON_EXECUTABLE, CUSTOM_CORES, CONFIG_FILE))
            exit_code, actual_output = self.execute(command)
            self.assertEqual(0, exit_code)
            self.assertRegexpMatches(actual_output,
                                     r'--cores {}'.format(CUSTOM_CORES))
