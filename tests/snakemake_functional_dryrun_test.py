import os
import pandas as pd
import re
import subprocess

from tests import testing_utils #local module

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
    sm_command_fmt = "snakemake --snakefile {} --configfile config_test_align_qc.yaml -n {}"
    sm_command = sm_command_fmt.format(
        os.path.join("Watermelon", 'deseq2.smk'),
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
