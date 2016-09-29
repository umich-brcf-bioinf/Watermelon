#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

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

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TEST_DIR), 'scripts')
class SnakeFileTest(unittest.TestCase):

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
        WATERMELON_SCRIPT = os.path.join(SCRIPTS_DIR, 'watermelon')
        shutil.copy(WATERMELON_SCRIPT, temp_dir_path)
        return temp_dir_path

    def test_watermelon_command_basecase(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_FILE=os.path.join(TEST_DIR, 'config.yaml')
            shutil.copy(CONFIG_FILE, temp_dir_path)
            SNAKEMAKE_FILE_SOURCE=os.path.join(TEST_DIR, 'test.snakefile')
            SNAKEMAKE_FILE_DEST=os.path.join(temp_dir_path, 'rnaseq.snakefile')
            shutil.copy(SNAKEMAKE_FILE_SOURCE, SNAKEMAKE_FILE_DEST)

            command = './watermelon {}'.format('config.yaml')
            exit_code, actual_output = self.execute(command)

            self.assertEqual(0, exit_code)

            log_file_name = os.path.join(temp_dir_path, 'logs', 'watermelon.log')
            self.assertTrue(os.path.exists(log_file_name))
            with open(log_file_name) as log_file:
                actual_log = "".join(log_file.readlines())

            with open('sample_A.txt') as output_file:
                actual_output_A = output_file.readlines()
            with open('sample_B.txt') as output_file:
                actual_output_B = output_file.readlines()


        self.assertEqual(['lineA\n'], actual_output_A)
        self.assertEqual(['lineB\n'], actual_output_B)
        self.assertRegexpMatches(actual_log, r'Provided cores: 40')
        self.assertRegexpMatches(actual_log, r'3 of 3 steps \(100%\) done')

    def test_watermelon_command_fails_when_invalid_config(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            with open('invalid_config.yaml', 'w') as config_file:
                print('an invalid config file', file=config_file)
            SNAKEMAKE_FILE_SOURCE=os.path.join(TEST_DIR, 'test.snakefile')
            SNAKEMAKE_FILE_DEST=os.path.join(temp_dir_path, 'rnaseq.snakefile')
            shutil.copy(SNAKEMAKE_FILE_SOURCE, SNAKEMAKE_FILE_DEST)

            command = './watermelon {}'.format('invalid_config.yaml')
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)

            log_file_name = os.path.join(temp_dir_path, 'logs', 'watermelon.log')
            self.assertTrue(os.path.exists(log_file_name))
            with open(log_file_name) as log_file:
                actual_log = "".join(log_file.readlines())

        self.assertRegexpMatches(actual_log, r'ERROR.*contact bfxcore support')

    def test_watermelon_command_fails_gracefully_if_no_snakemake(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            command = 'module purge; ./watermelon'
            exit_code, actual_output = self.execute(command)
            self.assertNotEqual(0, exit_code)
            self.assertRegexpMatches(actual_output, r'snakemake not found')


    def test_watermelon_command_showsUsageWhenConfigOmitted(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)

            command = './watermelon'
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)
            self.assertRegexpMatches(actual_output, 'Usage')

    def test_watermelon_command_showsUsageWhenExtraArgs(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)

            command = './watermelon foo bar'  
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)
            self.assertRegexpMatches(actual_output, 'Usage')

    def test_watermelon_command_showsUsageWhenFileNotFound(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)

            command = './watermelon i_do_not_exist.yaml'
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)
            self.assertRegexpMatches(actual_output, 'config file .* cannot be read')

    def test_watermelon_command_showsUsageWhenConfigIsDir(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = self.setup_tmp_dir(temp_dir)
            CONFIG_DIR=os.path.join(temp_dir_path, 'config_dir')
            os.mkdir(CONFIG_DIR)

            command = './watermelon {}'.format(CONFIG_DIR)
            exit_code, actual_output = self.execute(command)

            self.assertEqual(1, exit_code)
            self.assertRegexpMatches(actual_output, 'config file .* cannot be read')
