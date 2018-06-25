from __future__ import print_function, absolute_import
import filecmp
from glob import glob
import gzip
import os
import shutil
import subprocess
import unittest

from testfixtures import TempDirectory

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SNAKEFILE_PATH = os.path.join(TEST_DIR, '..', '..', '..', 'rnaseq.snakefile')
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '

class SnakemakeDryRunTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def Xtest_basecase(self):
        configfile_path = os.path.join(TEST_DIR, 'config.yaml')
        source_working_dir = os.path.join(TEST_DIR, 'working_dir')
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            tmp_actual_dir = os.path.join(temp_dir_path, 'actual')
            shutil.copytree(source_working_dir, tmp_actual_dir)

            os.chdir(tmp_actual_dir)
            command = \
'''snakemake -p --cores 2 \
     --snakefile {} \
     --configfile {} \
     --dryrun \
    {}
'''.format(SNAKEFILE_PATH, configfile_path, REDIRECT_OUTPUT)
            actual_output = subprocess.check_output(command,
                                                    shell=True)#,
                                                    #stderr= subprocess.STDOUT)

        self.assertRegexp('foo', str(actual_output))
