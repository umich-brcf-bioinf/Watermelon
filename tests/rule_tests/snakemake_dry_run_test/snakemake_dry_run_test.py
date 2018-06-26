from __future__ import print_function, absolute_import
# import filecmp
# from glob import glob
# import gzip
import os
import shutil
import subprocess
import unittest

from testfixtures import TempDirectory

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SNAKEFILE_PATH = os.path.join(TEST_DIR, '..',
'..',
'..',
'rnaseq.snakefile')
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '

class SnakemakeDryRunTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_basecase(self):
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
                                                    shell=True)
        actual_output = actual_output.decode("utf-8")
        lines = actual_output.split('\n')
        lines.pop()
        total_count = lines.pop().strip()
        actual_job_counts = dict([(k,v) for v,k in map(str.split, lines[-33:])])
        actual_job_counts['jobs'] = total_count
        #self.maxDiff=None
        expected_job_counts = {
            'jobs': '123',
            'align_build_tophat_sample_options': '8',
            'align_concat_reads': '8',
            'align_create_transcriptome_index': '1',
            'align_cutadapt_SE': '8',
            'align_deliverables_alignment': '1',
            'align_fastqc_tophat_align': '8',
            'align_fastqc_trimmed_reads': '8',
            'align_qc': '1',
            'align_tophat': '8',
            'all': '1',
            'deliverables_combined_summary': '1',
            'deliverables_deseq2': '1',
            'deliverables_tuxedo': '1',
            'deseq2_annotation': '7',
            'deseq2_diffex': '1',
            'deseq2_excel': '7',
            'deseq2_htseq': '8',
            'deseq2_htseq_merge': '1',
            'deseq2_metadata_contrasts': '1',
            'deseq2_run_info': '1',
            'deseq2_summary': '1',
            'tuxedo_annotate': '4',
            'tuxedo_cuffdiff': '4',
            'tuxedo_cummerbund': '4',
            'tuxedo_excel': '7',
            'tuxedo_flag': '4',
            'tuxedo_flip': '4',
            'tuxedo_group_replicates': '4',
            'tuxedo_last_split': '1',
            'tuxedo_run_info': '1',
            'tuxedo_split': '7',
            'tuxedo_summary': '1'}
        self.assertEqual(expected_job_counts, actual_job_counts)
