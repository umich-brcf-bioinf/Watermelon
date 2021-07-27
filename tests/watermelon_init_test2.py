# cd /nfs/mm-isilon/bioinfcore/ActiveProjects/cgates/Watermelon
# pytest tests/watermelon_init_test2.py

import os
import pytest
import re

import pandas as pd
from collections import namedtuple

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import watermelon_init


def test_get_analyst_name():
    # Just to be sure of what we're testing here:
    assert(watermelon_init._DEFAULT_ANALYST_INFO == "/nfs/turbo/umms-brcfpipeline/pipelines/analyst_info.csv")
    # base case
    analyst_name = watermelon_init._get_analyst_name(
        analyst_info_csv=watermelon_init._DEFAULT_ANALYST_INFO,
        user="cgates"
    )
    assert(analyst_name == "Chris Gates")
    # Missing file - warn and use default value
    with pytest.warns(UserWarning, match="Could not read"):
        analyst_name = watermelon_init._get_analyst_name("", "cgates")
        assert(analyst_name == "Analyst")
    # Unknown user key - warn and use default value
    with pytest.warns(UserWarning, match="Couldn't find"):
        analyst_name = watermelon_init._get_analyst_name(watermelon_init._DEFAULT_ANALYST_INFO, "foobar")
        assert(analyst_name == "Analyst")


def test_get_template():
    # Use this to create args object to pass in
    ArgStruct = namedtuple("ArgStruct", "type x_alt_template")
    # base case align_qc
    args = ArgStruct(type="align_qc", x_alt_template=None)
    template_cfg = watermelon_init.get_template(
        args=args,
        pipe_root=watermelon_init._WATERMELON_ROOT
    )
    template_keys = sorted(list(template_cfg.keys()))
    assert(template_keys == ['fastq_screen', 'report_info', 'trimming'])
    # Test file not found
    with pytest.raises(RuntimeError, match="Could not read"):
        template_cfg = watermelon_init.get_template(args, "foobar")


def test_validate_genomes():
    ref_dict = watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, "GRCh38", "align_qc").get("references")
    # Base case - no exceptions should be raised (see the assert False in the except)
    raised = False
    try:
        watermelon_init.validate_genomes(ref_dict, "align_qc")
    except Exception:
        raised = True
    assert(raised == False)
    # Invalid case - fake filepath
    ref_dict["fasta"] = "foo"
    raised = False
    try:
        watermelon_init.validate_genomes(ref_dict, "align_qc") # This time it will raise an exception
    except Exception:
        raised = True
    assert(raised == True)

def test_validate_fastq_dirs():
    # Assemble an in-memory samplesheet for data contained in this repo
    fastq_basedir = os.path.join(watermelon_init._WATERMELON_ROOT, "data", "sim_reads_human_chr22")
    header = ["sample,input_dir"]
    content = ["sample_0{},{}/sample_0{}".format(x, fastq_basedir, x) for x in range(1,7)]
    samplesheet_str = "\n".join(header + content)
    samplesheet_file_like = StringIO(samplesheet_str)
    ss_df = pd.read_csv(samplesheet_file_like)
    # Base case
    watermelon_init.validate_fastq_dirs(
        ss_df=ss_df,
        sample_col="sample",
        fq_col="input_dir"
    )
    # Test dir without fastq's
    ss_df.at[0, "input_dir"] = os.path.dirname(ss_df.at[0, "input_dir"])
    with pytest.raises(RuntimeError, match = "No fastq files found"):
        watermelon_init.validate_fastq_dirs(
            ss_df=ss_df,
            sample_col="sample",
            fq_col="input_dir"
        )

def test_set_up_refs_validKey():
    ref_dict = watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, "GRCh38", "align_qc")
    assert(sorted(ref_dict.keys()) == ["fastq_screen", "genome", "references"])


def test_set_up_refs_invalidKey():
    with pytest.raises(RuntimeError) as e_info:
        watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, "foobar", "align_qc")
