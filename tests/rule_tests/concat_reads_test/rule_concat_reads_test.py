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
SNAKEFILE_PATH = os.path.join(TEST_DIR, '..', '..', '..', 'rnaseq.snakefile')

def gunzip(source_file_pattern):
    def _gunzip_file(source_file, dest_file):
        with gzip.GzipFile(source_file, 'rb') as inF, \
             open(dest_file, 'wb') as outF:
            data = inF.read()
            outF.write(data)
            os.remove(source_file)
    for source_filename in glob(source_file_pattern):
        dest_filename = source_filename.rstrip('.gz')
        _gunzip_file(source_filename, dest_filename)

class ConcatReadsTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_basecase(self):
        anomalies = []
        configfile_path = os.path.join(TEST_DIR, 'basecase', 'basecase.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'basecase', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'basecase', 'expected')
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir) 

            os.chdir(tmp_actual_dir)
            redirect_output = '2>/dev/null'
            command = '''snakemake --cores 2 \
     --snakefile {} \
     --configfile {} \
     --force 01-raw_reads/Sample_0_R1.fastq.gz 01-raw_reads/Sample_1_R1.fastq.gz {}
'''.format(SNAKEFILE_PATH, configfile_path, redirect_output)
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
