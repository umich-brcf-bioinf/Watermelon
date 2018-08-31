from __future__ import print_function, absolute_import
import filecmp
from glob import glob
import gzip
import os
import shutil
import subprocess
import sys
import unittest

from testfixtures import TempDirectory

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
ROOT_DIR = os.path.realpath(os.path.join(TEST_DIR, '..', '..', '..'))
SCRIPTS_DIR = os.path.join(ROOT_DIR, 'scripts')
SNAKEFILE_PATH = os.path.join(ROOT_DIR, 'rnaseq.snakefile')

DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '

def _all_left_only_files(dircmp):
    left_only = sorted(dircmp.left_only)
    for dir_name, dc in sorted(dircmp.subdirs.items()):
        left_only.extend([os.path.join(dir_name, lo) for lo in _all_left_only_files(dc)])
    return left_only

def _all_diff_files(dircmp):
    diff_files = sorted(dircmp.diff_files)
    for dir_name, dc in sorted(dircmp.subdirs.items()):
        diff_files.extend([os.path.join(dir_name, lo) for lo in _all_diff_files(dc)])
    return diff_files

def missing_or_different_from_expected(expected_dir, actual_dir):
    dircmp = filecmp.dircmp(expected_dir, actual_dir)
    diffs = ['different: ' + diff_file for diff_file in _all_diff_files(dircmp)]
    diffs.extend(['only in expected: ' + expected_only for expected_only in _all_left_only_files(dircmp)])
    return diffs

class StringtieTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def _snakemake(self, configfile_path, source_expected_dir, source_working_dir):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path  # '/tmp/foo'
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir)

            os.chdir(tmp_actual_dir)
            command = '''snakemake --use-conda -p --cores 2 \
     --snakefile {} \
     --configfile {} \
     --force alignment_results/06-stringtie/sampleA.gtf alignment_results/06-stringtie/sampleB.gtf {}
'''.format(SNAKEFILE_PATH, configfile_path, REDIRECT_OUTPUT)
            subprocess.check_output(command, shell=True)

            anomalies = missing_or_different_from_expected(tmp_expected_dir,
                                                           tmp_actual_dir)
        return anomalies

    def test_basecase(self):
        configfile_path = os.path.join(TEST_DIR, 'basecase', 'basecase.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'basecase', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'basecase', 'expected')
        anomalies =  self._snakemake(configfile_path, source_expected_dir, source_working_dir)
        self.assertEqual([], anomalies, 'some files/dirs did not match expected')
