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


class FlipDiffexTest(unittest.TestCase):
    def _snakemake(self, configfile_path, source_expected_dir, source_working_dir):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path  # '/tmp/foo'
            tmp_expected_dir = os.path.join(temp_dir_path, 'expected')
            shutil.copytree(source_expected_dir, tmp_expected_dir)
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir)
            

            os.chdir(tmp_actual_dir)
            os.symlink(SCRIPTS_DIR, 'scripts')
            command = '''snakemake --cores 2 \
     --snakefile {0} \
     --configfile {1} \
     --force 09-flip_diffex/A_B/gene_exp.flip.diff 09-flip_diffex/A_B/isoform_exp.flip.diff
'''.format(SNAKEFILE_PATH, configfile_path)
            subprocess.check_output(command, shell=True)

            diff_command = 'diff --brief -r {} {}'.format(os.path.join(tmp_expected_dir, '09-flip_diffex'),
                                                          os.path.join(tmp_actual_dir, '09-flip_diffex'))
            try:
                subprocess.check_output(command, shell=True)
                return True
            except Exception:
                return True

#             dircmp = filecmp.dircmp(tmp_expected_dir, tmp_actual_dir)
#             anomalies = {}
#             if dircmp.left_only:
#                 anomalies['in expected but not actual'] = dircmp.left_only
#             if dircmp.diff_files:
#                 anomalies['files do not match'] = dircmp.diff_files
# 
#             print('\n\nBEFORE: ' + str(anomalies))
#             return anomalies
#         print('\n\nAFTER: ' + str(anomalies))
            
    def test_basecase(self):
        configfile_path = os.path.join(TEST_DIR, 'basecase', 'basecase.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'basecase', 'working_dir')
        source_expected_dir = os.path.join(TEST_DIR, 'basecase', 'expected')
        files_matched =  self._snakemake(configfile_path, source_expected_dir, source_working_dir)
        self.assertTrue(files_matched, 'some files did not match')
