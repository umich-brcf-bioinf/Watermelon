import os
import pandas as pd
import pytest
import subprocess

from mock import patch
from tempfile import TemporaryDirectory
from tests import testing_utils #local module

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))


def assertDataframeEqualEnough(a, b):
    #https://stackoverflow.com/a/54608841
    # allows us to use pandas testing.assert_frame_equal with mock.patch
    # Note using check_exact=False
    try:
        pd.testing.assert_frame_equal(a, b, check_exact=False)
    except AssertionError as e:
        pytest.fail("pd.testing.assert_frame_equal with check_exact=False did not pass")
    return True


@patch("pandas.DataFrame.__eq__", assertDataframeEqualEnough)
def test_wminit_snakemake_run_align_qc(tmp_path):
    os.chdir(tmp_path)

    wminit_command_fmt = "{} --genome_build TestData --project_id test_align_qc --type align_qc --input_run_dirs {}"
    wminit_command = wminit_command_fmt.format(
        os.path.join(WATERMELON_BASE_DIR, "watermelon_init.py"),
        os.path.join(WATERMELON_BASE_DIR, "data", "sim_reads_human_chr22")
    )
    try:
        # 0 if good, otherwise CalledProcessError will be raised
        return_code = subprocess.check_call(wminit_command, shell=True)
    except subprocess.CalledProcessError as e:
        return_code = e.returncode

    assert(return_code == 0)

    rm_keys = ["email", "trimming"] # Unnecessary to trim the simulated data
    #Create modified config in this temp dir, using example config and just removing the email key
    testing_utils.create_modified_config("config_test_align_qc.yaml", "config_test_align_qc.ready", {}, rm_keys)
    command_fmt = "snakemake --snakefile {} --configfile {} --profile {}"
    # --singularity-args \'-B {}\'
    command = command_fmt.format(
        os.path.join("Watermelon", 'align_qc.smk'),
        "config_test_align_qc.ready",
        os.path.join("Watermelon", "config", "profile-greatlakes")
    )

    return_code = subprocess.call(command, shell=True)

    assert(return_code == 0)

    #Things to check:
    # un-annotated rsem counts tables - FPKM, TPM (alignment_results)
    gene_FPKM_actual_path = os.path.join(tmp_path, 'analysis_test_align_qc', 'alignment_results', '05-combine_counts', 'gene_FPKM.txt')
    gene_FPKM_actual = pd.read_csv(gene_FPKM_actual_path, sep="\t")
    gene_FPKM_expected_path = os.path.join(tmp_path, 'Watermelon', 'data', 'expected_results_files', 'gene_FPKM.txt')
    gene_FPKM_expected = pd.read_csv(gene_FPKM_expected_path, sep="\t")

    gene_TPM_actual_path = os.path.join(tmp_path, 'analysis_test_align_qc', 'alignment_results', '05-combine_counts', 'gene_TPM.txt')
    gene_TPM_actual = pd.read_csv(gene_TPM_actual_path, sep="\t")
    gene_TPM_expected_path = os.path.join(tmp_path, 'Watermelon', 'data', 'expected_results_files', 'gene_TPM.txt')
    gene_TPM_expected = pd.read_csv(gene_TPM_expected_path, sep="\t")

    assert(gene_FPKM_expected == gene_FPKM_actual)
    assert(gene_TPM_expected == gene_TPM_actual)

            # # deseq2 counts tables (deliverables)
            # deseq2_raw_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'counts', 'deseq2_raw_counts.txt')
            # deseq2_raw_actual = pd.read_csv(deseq2_raw_actual_path, sep="\t")
            # deseq2_raw_expected_path = os.path.join(temp_dir.path, 'data', 'expected_results_files', 'deseq2_raw_counts.txt')
            # deseq2_raw_expected = pd.read_csv(deseq2_raw_expected_path, sep="\t")
            #
            # deseq2_rlog_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'counts', 'deseq2_rlog_normalized_counts.txt')
            # deseq2_rlog_actual = pd.read_csv(deseq2_rlog_actual_path, sep="\t")
            # deseq2_rlog_expected_path = os.path.join(temp_dir.path, 'data', 'expected_results_files', 'deseq2_rlog_normalized_counts.txt')
            # deseq2_rlog_expected = pd.read_csv(deseq2_rlog_expected_path, sep="\t")
            #
            # deseq2_depthnorm_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'counts', 'deseq2_depth_normalized_counts.txt')
            # deseq2_depthnorm_actual = pd.read_csv(deseq2_depthnorm_actual_path, sep="\t")
            # deseq2_depthnorm_expected_path = os.path.join(temp_dir.path, 'data', 'expected_results_files', 'deseq2_depth_normalized_counts.txt')
            # deseq2_depthnorm_expected = pd.read_csv(deseq2_depthnorm_expected_path, sep="\t")
            #
            # self.assertEqual(deseq2_raw_expected, deseq2_raw_actual)
            # self.assertEqual(deseq2_rlog_expected, deseq2_rlog_actual)
            # self.assertEqual(deseq2_depthnorm_expected, deseq2_depthnorm_actual)
            #
            # # gene lists (deliverables, txt)
            # gene_list_actual_path = os.path.join(temp_dir.path, 'analysis_test', 'deliverables', 'deseq2', 'gene_lists', 'model_treatment', 'drug_v_control.annot.txt')
            # gene_list_actual = pd.read_csv(gene_list_actual_path, sep="\t")
            # gene_list_expected_path = os.path.join(temp_dir.path, 'data', 'expected_results_files', 'drug_v_control.annot.txt')
            # gene_list_expected = pd.read_csv(gene_list_expected_path, sep="\t")
            #
            # self.assertEqual(gene_list_expected, gene_list_actual)
