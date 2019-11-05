from __future__ import print_function, absolute_import
import filecmp
from glob import glob
import os
import shutil
import subprocess
import sys
import unittest

from testfixtures import TempDirectory

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..', '..', '..'))
SNAKEFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'rnaseq.snakefile')
EXAMPLE_CONFIGFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'config', 'example_config.yaml')
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '

sys.path.append(os.path.join(WATERMELON_BASE_DIR, 'tests'))
from testing_utils import create_modified_config, gunzip #local module


class CutadaptTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    #TODO: use a string constant for "02-cutadapt"
    def _snakemake(self, source_expected_dir, source_working_dir, snakemake_flags, config_replacement_vals, config_rm_keys=None):
        anomalies = []
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir)

            os.chdir(tmp_actual_dir)

            #Create modified config in this temp dir, using example config and replacing values as needed
            configfile_path = os.path.join(temp_dir_path, 'testcase_config.yaml')
            create_modified_config(EXAMPLE_CONFIGFILE_PATH, configfile_path, config_replacement_vals, config_rm_keys)

            command_fmt = ('snakemake -p --cores 2 '
                '--use-conda '
                '--conda-prefix /nfs/med-bfx-common/pipelines/Watermelon/Watermelon-seedless/envs/built ' # FIXME: This won't work with CI - however, without it will take ages to build the env(s)
                '--config skip_validation=True '
                '--snakefile {snakefile_path} '
                '--configfile {configfile_path} '
                '{snakemake_flags} {redirect_output}')
            command = command_fmt.format(
                snakefile_path=SNAKEFILE_PATH,
                configfile_path=configfile_path,
                snakemake_flags=snakemake_flags,
                redirect_output=REDIRECT_OUTPUT)
            subprocess.check_output(command, shell=True)

            gunzip(tmp_expected_dir + '/alignment_results/02-cutadapt/*.gz')
            gunzip('alignment_results/02-cutadapt/*.gz')

            for expected_filename in glob(tmp_expected_dir + '/alignment_results/02-cutadapt/*'):
                basename = os.path.basename(expected_filename)
                actual_filename = tmp_actual_dir + '/alignment_results/02-cutadapt/' + basename
                try:
                    file_matched = filecmp.cmp(expected_filename, actual_filename, shallow=True)
                    if not file_matched:
                        anomalies.append(basename + ": different than expected")
                except FileNotFoundError:
                    anomalies.append(basename + ": missing")
        return anomalies


    def test_end_trim(self):
        source_working_dir = os.path.join(TEST_DIR, 'end_trim', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'end_trim', 'expected')
        #Replace values from the example config with these
        replacement_vals = {
            'dirs': {
                'input': '',
                'alignment_output': 'alignment_results'
            },
            'email' : {
                'to': 'nobody'
            },
            'trimming_options' : {
                'base_quality_5prime' : 0,
                'base_quality_3prime' : 0,
                'trim_length_5prime' : 3,
                'trim_length_3prime' : 10
            }
        }
        snakemake_flags = ('--force '
                          'alignment_results/02-cutadapt/Sample_0_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_1_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_2_trimmed_R1_PE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_2_trimmed_R2_PE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_3_trimmed_R2_SE.fastq.gz')
        anomalies =  self._snakemake(source_expected_dir,
                                     source_working_dir,
                                     snakemake_flags,
                                     replacement_vals)
        self.assertEqual([], anomalies, 'some files did not match')

    def test_no_trim(self):
        '''test_no_trim (rule_cutadapt_test.CutadaptTest)
        This runs cutadapt, but with trimming_options all set to 0'''
        # TWS - probably should replace this test with pseudotrim test, commented out below
        source_working_dir = os.path.join(TEST_DIR, 'no_trim', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'no_trim', 'expected')
        #Replace values from the example config with these
        replacement_vals = {
            'dirs': {
                'input': '',
                'alignment_output': 'alignment_results'
            },
            'email' : {
                'to': 'nobody'
            },
            'trimming_options' : {
                'base_quality_5prime' : 0,
                'base_quality_3prime' : 0,
                'trim_length_5prime' : 0,
                'trim_length_3prime' : 0
            }
        }
        snakemake_flags = ('--force '
                          'alignment_results/02-cutadapt/Sample_0_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_1_trimmed_R1_SE.fastq.gz ')
        anomalies =  self._snakemake(source_expected_dir,
                                     source_working_dir,
                                     snakemake_flags,
                                     replacement_vals)
        self.assertEqual([], anomalies, 'some files did not match')

    # def test_pseudotrim(self):
    #     '''test_pseudotrim (rule_cutadapt_test.CutadaptTest)
    #     This does not run cutadapt, instead links to the original fastqs (pseudotrim),
    #     if no trimming options are set, this is the outcome'''
    #     configfile_path = os.path.join(TEST_DIR, 'no_trim', 'no_trim.yaml')
    #     source_working_dir = os.path.join(TEST_DIR, 'no_trim', 'working_dir')
    #     source_expected_dir = os.path.join(TEST_DIR, 'no_trim', 'expected')
    #     #Replace values from the example config with these
    #     replacement_vals = {
    #         'dirs': {
    #             'input': '',
    #             'alignment_output': 'alignment_results'
    #         },
    #         'email' : {
    #             'to': 'nobody'
    #         }
    #     }
    #     #Remove trimming_options from the config
    #     rm_keys = ['trimming_options']
    #     snakemake_flags = ('--force '
    #                       'alignment_results/02-cutadapt/Sample_0_trimmed_R1_SE.fastq.gz '
    #                       'alignment_results/02-cutadapt/Sample_1_trimmed_R1_SE.fastq.gz ')
    #     anomalies =  self._snakemake(source_expected_dir,
    #                                  source_working_dir,
    #                                  snakemake_flags,
    #                                  replacement_vals,
    #                                  rm_keys)
    #     self.assertEqual([], anomalies, 'some files did not match')

    def test_quality_trim(self):
        source_working_dir = os.path.join(TEST_DIR, 'quality_trim', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'quality_trim', 'expected')
        #Replace values from the example config with these
        replacement_vals = {
            'dirs': {
                'input': '',
                'alignment_output': 'alignment_results'
            },
            'email' : {
                'to': 'nobody'
            },
            'trimming_options' : {
                'base_quality_5prime' : 20,
                'base_quality_3prime' : 20,
                'trim_length_5prime' : 0,
                'trim_length_3prime' : 0
            }
        }
        snakemake_flags = ('--force '
                          'alignment_results/02-cutadapt/Sample_0_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_1_trimmed_R1_SE.fastq.gz ')
        anomalies =  self._snakemake(source_expected_dir,
                                     source_working_dir,
                                     snakemake_flags,
                                     replacement_vals)
        self.assertEqual([], anomalies, 'some files did not match')
