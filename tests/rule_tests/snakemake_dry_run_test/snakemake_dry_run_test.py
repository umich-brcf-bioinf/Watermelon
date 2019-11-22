import os
import shutil
import subprocess
import unittest

from testfixtures import TempDirectory
from tests import testing_utils #local module

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..', '..', '..'))
SNAKEFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'rnaseq.snakefile')
EXAMPLE_CONFIGFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'config', 'example_config.yaml')
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' #if DEBUG else ' 2>/dev/null '

class SnakemakeDryRunTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_dryrun_passes(self):
        with TempDirectory() as temp_dir:
            os.chdir(temp_dir.path)

            example_data_dir = os.path.join(WATERMELON_BASE_DIR, 'data')
            tmp_data_dir = os.path.join(temp_dir.path, 'data')
            #Copy example data here
            shutil.copytree(example_data_dir, tmp_data_dir)
            #Copy samplesheet here
            shutil.copyfile(os.path.join(WATERMELON_BASE_DIR, 'config', 'example_samplesheet.csv'), 'samplesheet.csv')

            config_replacement_vals = {
                'sample_description_file': os.path.join(temp_dir.path, 'samplesheet.csv'),
                'dirs': {'input' : os.path.join(tmp_data_dir, 'sim_reads_human_chr22')},
                'references': {
                    'fasta' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.98.chr22.gtf.gz'), # These just need to exist for dry-run
                    'gtf' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.dna_sm.chr22.fa.gz') # For real-deal, unzip them first & point to non-gz
                },
                'email' : {
                    'to': 'nobody'
                },
            }

            #Create modified config in this temp dir, using example config and replacing values as needed
            configfile_path = os.path.join(temp_dir.path, 'testcase_config.yaml')
            testing_utils.create_modified_config(EXAMPLE_CONFIGFILE_PATH, configfile_path, config_replacement_vals)

            command_fmt = 'snakemake --snakefile {} --configfile {} -n --config skip_validation=True {}'
            command = command_fmt.format(SNAKEFILE_PATH, configfile_path, REDIRECT_OUTPUT)
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
                output = e.output.decode()
                print(output)

        self.assertEqual(0, return_code)
        self.assertRegex(dryrun_line, 'This was a dry-run')
