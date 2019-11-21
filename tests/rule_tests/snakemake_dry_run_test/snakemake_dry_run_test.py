import os
import subprocess
import unittest

from testfixtures import TempDirectory

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..', '..', '..'))
SNAKEFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'rnaseq.snakefile')
EXAMPLE_CONFIGFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'config', 'example_config.yaml')
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '

class SnakemakeDryRunTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_dryrun_passes(self):
        with TempDirectory() as temp_dir:
            os.chdir(temp_dir.path)
            command_fmt = 'snakemake --snakefile {} --configfile {} -n --config skip_validation=True {}'
            command = command_fmt.format(SNAKEFILE_PATH, EXAMPLE_CONFIGFILE_PATH, REDIRECT_OUTPUT)
            print(command)
            try:
                #return code only available from subprocess.check_output if non-zero (raises CalledProcessError)
                return_code = 0
                actual_output = subprocess.check_output(command, shell=True)
                actual_output = actual_output.decode("utf-8")
                lines = actual_output.split('\n')
                empty_last_line = lines.pop()
                dryrun_line = lines.pop()
            except subprocess.CalledProcessError as e:
                return_code = e.returncode

        self.assertEqual(0, return_code)
        self.assertRegex(dryrun_line, 'This was a dry-run')
