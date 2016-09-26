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

class CutadaptTest(unittest.TestCase):
    #TODO: use a string constant for "02-cutadapt"
    def test_end_trim(self):
        anomalies = []
        configfile_path = os.path.join(TEST_DIR, 'end_trim_test', 'end_trim.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'end_trim_test', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'end_trim_test', 'expected')
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir) 

            os.chdir(tmp_actual_dir)
            command = '''snakemake --cores 2 \
     --snakefile {0} \
     --configfile {1} \
     --force 02-cutadapt/Sample_0_trimmed_R1.fastq.gz 02-cutadapt/Sample_1_trimmed_R1.fastq.gz
'''.format(SNAKEFILE_PATH, configfile_path)
            subprocess.check_output(command, shell=True)

            gunzip(tmp_expected_dir + '/02-cutadapt/*.gz')
            gunzip('02-cutadapt/*.gz')

            for expected_filename in glob(tmp_expected_dir + '/02-cutadapt/*'):
                basename = os.path.basename(expected_filename)
                actual_filename = tmp_actual_dir + '/02-cutadapt/' + basename
                try:
                    file_matched = filecmp.cmp(expected_filename, actual_filename, shallow=True)
                    if not file_matched:
                        anomalies.append(basename + ": different than expected")
                except FileNotFoundError:
                    anomalies.append(basename + ": missing")

        self.assertEqual([], anomalies, 'some files did not match')
