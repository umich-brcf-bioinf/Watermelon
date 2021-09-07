import os
import pandas as pd
import re
import subprocess

from tests import testing_utils #local module
from scripts import rnaseq_snakefile_helper as helper

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))
DEBUG = 'WATERMELON_DEBUG' in os.environ
REDIRECT_OUTPUT = ' ' if DEBUG else ' 2>/dev/null '

def test_run_wminit_align_qc(tmp_path):
    os.chdir(tmp_path)
    command_fmt = "{} --genome_build TestData --project_id test_align_qc --type align_qc --input_run_dirs {}"
    command = command_fmt.format(
        os.path.join(WATERMELON_BASE_DIR, "watermelon_init.py"),
        os.path.join(WATERMELON_BASE_DIR, "data", "sim_reads_human_chr22"))
    try:
        # 0 if good, otherwise CalledProcessError will be raised
        return_code = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        return_code = e.returncode
    assert(return_code == 0)

def test_wminit_snakemake_dryrun_align_qc(tmp_path):
    os.chdir(tmp_path)
    wminit_command_fmt = "{} --genome_build TestData --project_id test_align_qc --type align_qc --input_run_dirs {}"
    wminit_command = wminit_command_fmt.format(
        os.path.join(WATERMELON_BASE_DIR, "watermelon_init.py"),
        os.path.join(WATERMELON_BASE_DIR, "data", "sim_reads_human_chr22"))
    try:
        # 0 if good, otherwise CalledProcessError will be raised
        return_code = subprocess.check_call(wminit_command, shell=True)
    except subprocess.CalledProcessError as e:
        return_code = e.returncode
    assert(return_code == 0)
    sm_command_fmt = "snakemake --snakefile {} --configfile config_test_align_qc.yaml -n {}"
    sm_command = sm_command_fmt.format(
        os.path.join(WATERMELON_BASE_DIR, 'align_qc.smk'),
        REDIRECT_OUTPUT)
    try:
        #return code only available from subprocess.check_output if non-zero (raises CalledProcessError)
        return_code = 0
        actual_output = subprocess.check_output(sm_command, shell=True)
        actual_output = actual_output.decode("utf-8")
        lines = actual_output.split('\n')
        empty_last_line = lines.pop()
        dryrun_line = lines.pop()
    except subprocess.CalledProcessError as e:
        return_code = e.returncode
        output = e.output.decode()
    assert(return_code == 0)
    assert(re.match('This was a dry-run', dryrun_line) != None)

def test_wminit_snakemake_dryrun_diffex(tmp_path):
    os.chdir(tmp_path)
    wminit_command_fmt = "{} --genome_build TestData --project_id test_align_qc --type diffex --sample_sheet {} --count_matrix {}"
    wminit_command = wminit_command_fmt.format(
        os.path.join(WATERMELON_BASE_DIR, "watermelon_init.py"),
        os.path.join(WATERMELON_BASE_DIR, "config", "example_samplesheet.csv"),
        os.path.join(WATERMELON_BASE_DIR, "data", "expected_results_files", "gene_expected_count.txt"))
    try:
        # 0 if good, otherwise CalledProcessError will be raised
        return_code = subprocess.check_call(wminit_command, shell=True)
    except subprocess.CalledProcessError as e:
        return_code = e.returncode
    assert(return_code == 0)
    #import pdb; pdb.set_trace()
    sm_command_fmt = "snakemake --snakefile {} --configfile config_test_align_qc.yaml -n {}"
    sm_command = sm_command_fmt.format(
        os.path.join(WATERMELON_BASE_DIR, 'deseq2.smk'),
        REDIRECT_OUTPUT)
    try:
        #return code only available from subprocess.check_output if non-zero (raises CalledProcessError)
        return_code = 0
        actual_output = subprocess.check_output(sm_command, shell=True)
        actual_output = actual_output.decode("utf-8")
        lines = actual_output.split('\n')
        empty_last_line = lines.pop()
        dryrun_line = lines.pop()
    except subprocess.CalledProcessError as e:
        return_code = e.returncode
        output = e.output.decode()
        print(output)
    assert(return_code == 0)
    assert(re.match('This was a dry-run', dryrun_line) != None)

# class SnakemakeExampleDataTest(unittest.TestCase):
#     def setUp(self):
#         self.original_wd = os.getcwd()
#         self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqualEnough)
#
#     def tearDown(self):
#         os.chdir(self.original_wd)
#
#     def assertDataframeEqualEnough(self, a, b, msg):
#         #https://stackoverflow.com/a/54344148
#         # allows us to use pandas testing.assert_frame_equal within unittest.testcase
#         # Note using check_exact=False
#         try:
#             pd.testing.assert_frame_equal(a, b, check_exact=False)
#         except AssertionError as e:
#             raise self.failureException(msg) from e
#     def test_snakemake_exampledata(self):
#         with TempDirectory() as temp_dir:
#             os.chdir(temp_dir.path)
#
#             example_data_dir = os.path.join(WATERMELON_BASE_DIR, 'data')
#             tmp_data_dir = os.path.join(temp_dir.path, 'data')
#             #Copy example data here
#             shutil.copytree(example_data_dir, tmp_data_dir)
#             #Copy samplesheet here
#             shutil.copyfile(os.path.join(WATERMELON_BASE_DIR, 'config', 'example_samplesheet.csv'), 'samplesheet.csv')
#
#             #Unzip references
#             helper.gunzip_glob(os.path.join(tmp_data_dir, '*.gz'))
#
#             config_replacement_vals = {
#                 'samplesheet': os.path.join(temp_dir.path, 'samplesheet.csv'),
#                 'dirs': {
#                     'input' : os.path.join(tmp_data_dir, 'sim_reads_human_chr22'),
#                     'alignment_output' : os.path.join(temp_dir.path, 'analysis_test', 'alignment_results'),
#                     'diffex_output' : os.path.join(temp_dir.path, 'analysis_test', 'diffex_results'),
#                     'deliverables_output' : os.path.join(temp_dir.path, 'analysis_test', 'deliverables'),
#                 },
#                 'references': {
#                     'fasta' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.dna_sm.chr22.fa'),
#                     'gtf' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.98.chr22.gtf'),
#                     'annotation_tsv' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.98_annotation.tsv')
#                 },
#                 'trimming' : {
#                     'cutadapt_args' : '--nextseq-trim 13 -u 3 -m 20 --trim-n' # Not testing adapter trimming here (these reads are simulated)
#                 }
#             }
#
#             rm_keys = ['fastq_screen', 'email']
#
#             #Create modified config in this temp dir, using example config and replacing values as needed
#             configfile_path = os.path.join(temp_dir.path, 'testcase_config.yaml')
#             testing_utils.create_modified_config(EXAMPLE_CONFIGFILE_PATH, configfile_path, config_replacement_vals, rm_keys)
#
#             command_fmt = 'snakemake --cores 20 --use-singularity --singularity-args \'-B {}\' --snakefile {} --configfile {} --config skip_validation=True {}'
#             command = command_fmt.format(WATERMELON_BASE_DIR, SNAKEFILE_PATH, configfile_path, REDIRECT_OUTPUT)
#             print(command)
#
#             return_code = subprocess.call(command, shell=True)
#
#             # from nose.tools import set_trace; set_trace() # TWS DEBUG
#
#             self.assertEqual(0, return_code)
#
#             #Things to check:
#             # un-annotated rsem counts tables - FPKM, TPM (alignment_results)
#             gene_FPKM_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'alignment_results', '05-combine_counts', 'gene_FPKM.txt')
#             gene_FPKM_actual = pd.read_csv(gene_FPKM_actual_path, sep="\t")
#             gene_FPKM_expected_path = os.path.join(tmp_data_dir, 'expected_results_files', 'gene_FPKM.txt')
#             gene_FPKM_expected = pd.read_csv(gene_FPKM_expected_path, sep="\t")
#
#             gene_TPM_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'alignment_results', '05-combine_counts', 'gene_TPM.txt')
#             gene_TPM_actual = pd.read_csv(gene_TPM_actual_path, sep="\t")
#             gene_TPM_expected_path = os.path.join(tmp_data_dir, 'expected_results_files', 'gene_TPM.txt')
#             gene_TPM_expected = pd.read_csv(gene_TPM_expected_path, sep="\t")
#
#             # from nose.tools import set_trace; set_trace() # TWS DEBUG
#
#             self.assertEqual(gene_FPKM_expected, gene_FPKM_actual)
#             self.assertEqual(gene_TPM_expected, gene_TPM_actual)
#
#             # deseq2 counts tables (deliverables)
#             deseq2_raw_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'counts', 'deseq2_raw_counts.txt')
#             deseq2_raw_actual = pd.read_csv(deseq2_raw_actual_path, sep="\t")
#             deseq2_raw_expected_path = os.path.join(tmp_data_dir, 'expected_results_files', 'deseq2_raw_counts.txt')
#             deseq2_raw_expected = pd.read_csv(deseq2_raw_expected_path, sep="\t")
#
#             deseq2_rlog_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'counts', 'deseq2_rlog_normalized_counts.txt')
#             deseq2_rlog_actual = pd.read_csv(deseq2_rlog_actual_path, sep="\t")
#             deseq2_rlog_expected_path = os.path.join(tmp_data_dir, 'expected_results_files', 'deseq2_rlog_normalized_counts.txt')
#             deseq2_rlog_expected = pd.read_csv(deseq2_rlog_expected_path, sep="\t")
#
#             deseq2_depthnorm_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'counts', 'deseq2_depth_normalized_counts.txt')
#             deseq2_depthnorm_actual = pd.read_csv(deseq2_depthnorm_actual_path, sep="\t")
#             deseq2_depthnorm_expected_path = os.path.join(tmp_data_dir, 'expected_results_files', 'deseq2_depth_normalized_counts.txt')
#             deseq2_depthnorm_expected = pd.read_csv(deseq2_depthnorm_expected_path, sep="\t")
#
#             self.assertEqual(deseq2_raw_expected, deseq2_raw_actual)
#             self.assertEqual(deseq2_rlog_expected, deseq2_rlog_actual)
#             self.assertEqual(deseq2_depthnorm_expected, deseq2_depthnorm_actual)
#
#             # gene lists (deliverables, txt)
#             gene_list_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'deseq2', 'gene_lists', 'model_treatment', 'drug_v_control.annot.txt')
#             gene_list_actual = pd.read_csv(gene_list_actual_path, sep="\t")
#             gene_list_expected_path = os.path.join(tmp_data_dir, 'expected_results_files', 'drug_v_control.annot.txt')
#             gene_list_expected = pd.read_csv(gene_list_expected_path, sep="\t")
#
#             self.assertEqual(gene_list_expected, gene_list_actual)
