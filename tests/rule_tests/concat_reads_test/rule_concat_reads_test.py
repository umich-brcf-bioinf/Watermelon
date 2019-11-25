from __future__ import print_function, absolute_import
import filecmp
from glob import glob
import io
import os
import shutil
import subprocess
import sys
import unittest
from testfixtures import TempDirectory
from tests import testing_utils #local module


TEST_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..', '..', '..'))
SNAKEFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'rnaseq.snakefile')
EXAMPLE_CONFIGFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'config', 'example_config.yaml')
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '


class ConcatReadsTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_basecase(self):

        anomalies = []
        source_working_dir = os.path.join(TEST_DIR, 'basecase', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'basecase', 'expected')
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir)
            os.chdir(tmp_actual_dir)
            #Create modified config in this temp dir, using example config and replacing values as needed
            new_input = os.path.join(tmp_actual_dir, 'inputs', '00-multiplexed_reads')
            samplesheet_path = os.path.join(source_working_dir, 'test_samplesheet.csv')
            replacement_vals = {
                'samplesheet' : samplesheet_path,
                'dirs': {
                    'input': new_input,
                    'alignment_output': 'alignment_results'
                },
                'email' : {
                    'to': 'nobody'
                }
            }
            testing_utils.create_modified_config(EXAMPLE_CONFIGFILE_PATH, 'modified_config.yaml', replacement_vals)

            command_fmt = ('snakemake -p --cores 2 '
                '--config skip_validation=True '
                '--snakefile {} '
                '--configfile {} '
                '--force alignment_results/01-raw_reads/Sample_0_R1.fastq.gz '
                'alignment_results/01-raw_reads/Sample_1_R1.fastq.gz '
                'alignment_results/01-raw_reads/Sample_1_R2.fastq.gz '
                '{}')
            command = command_fmt.format(SNAKEFILE_PATH, 'modified_config.yaml', REDIRECT_OUTPUT)
            subprocess.check_output(command, shell=True)

            for expected_filename in glob(tmp_expected_dir + '/01-raw_reads/*'):
                basename = os.path.basename(expected_filename)
                actual_filename = tmp_actual_dir + '/01-raw_reads/' + basename
                try:
                    file_matched = filecmp.cmp(expected_filename, actual_filename, shallow=True)
                    if not file_matched:
                        anomalies.append(basename + ": different than expected")
                except FileNotFoundError:
                    anomalies.append(basename + ": missing")

        self.assertEqual([], anomalies, 'some files did not match')
