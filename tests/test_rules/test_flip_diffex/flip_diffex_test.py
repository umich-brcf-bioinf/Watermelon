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


class FlipDiffexTest(unittest.TestCase):
    def _snakemake(self, configfile_path, source_expected_dir, source_working_dir):
        anomalies = []
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path  # '/tmp/foo'
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir) 

            os.chdir(tmp_actual_dir)
            command = '''snakemake --cores 2 \
     --snakefile {0} \
     --configfile {1} \
     --force 09-flip_diffex/gene_exp_diff.flipped.txt 09-flip_diffex.isoform_diff.flipped.txt
'''.format(SNAKEFILE_PATH, configfile_path)
            subprocess.check_output(command, shell=True)

            for expected_filename in glob(tmp_expected_dir + '/09-flip_diffex/*'):
                basename = os.path.basename(expected_filename)
                actual_filename = tmp_actual_dir + '/09-flip_diffex/' + basename
                try:
                    file_matched = filecmp.cmp(expected_filename, actual_filename, shallow=True)
                    if not file_matched:
                        anomalies.append(basename + ": different than expected")
                except FileNotFoundError:
                    anomalies.append(basename + ": missing")
        return anomalies


    def test_basecase(self):
        configfile_path = os.path.join(TEST_DIR, 'basecase', 'basecase.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'end_trim', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'end_trim', 'expected')
        anomalies =  self._snakemake(configfile_path, source_expected_dir, source_working_dir)
        self.assertEqual([], anomalies, 'some files did not match')

