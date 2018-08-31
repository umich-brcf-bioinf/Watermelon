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
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '


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
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    #TODO: use a string constant for "02-cutadapt"
    def _snakemake(self, configfile_path, source_expected_dir, source_working_dir, snakemake_flags):
        anomalies = []
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir)

            os.chdir(tmp_actual_dir)
            command = '''snakemake -p --cores 2 \
     --snakefile {snakefile_path} \
     --configfile {configfile_path} \
     {snakemake_flags} {redirect_output}
'''.format(snakefile_path=SNAKEFILE_PATH,
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
        configfile_path = os.path.join(TEST_DIR, 'end_trim', 'end_trim.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'end_trim', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'end_trim', 'expected')
        snakemake_flags = ('--force '
                          'alignment_results/02-cutadapt/Sample_0_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_1_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_2_trimmed_R1_PE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_2_trimmed_R2_PE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_3_trimmed_R2_SE.fastq.gz')
        anomalies =  self._snakemake(configfile_path,
                                     source_expected_dir,
                                     source_working_dir,
                                     snakemake_flags)
        self.assertEqual([], anomalies, 'some files did not match')

    def test_no_trim(self):
        configfile_path = os.path.join(TEST_DIR, 'no_trim', 'no_trim.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'no_trim', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'no_trim', 'expected')
        snakemake_flags = ('--force '
                          'alignment_results/02-cutadapt/Sample_0_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_1_trimmed_R1_SE.fastq.gz ')
        anomalies =  self._snakemake(configfile_path,
                                     source_expected_dir,
                                     source_working_dir,
                                     snakemake_flags)
        self.assertEqual([], anomalies, 'some files did not match')

    def test_quality_trim(self):
        configfile_path = os.path.join(TEST_DIR, 'quality_trim', 'quality_trim.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'quality_trim', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'quality_trim', 'expected')
        snakemake_flags = ('--force '
                          'alignment_results/02-cutadapt/Sample_0_trimmed_R1_SE.fastq.gz '
                          'alignment_results/02-cutadapt/Sample_1_trimmed_R1_SE.fastq.gz ')
        anomalies =  self._snakemake(configfile_path,
                                     source_expected_dir,
                                     source_working_dir,
                                     snakemake_flags)
        self.assertEqual([], anomalies, 'some files did not match')