# cd /nfs/mm-isilon/bioinfcore/ActiveProjects/cgates/Watermelon
# pytest tests/watermelon_init_test2.py

from argparse import Namespace
import os
import os.path
import pytest
import re

import pandas as pd

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import watermelon_init

def test_generate_samplesheet_Basecase(tmp_path):
    fq_parent_dir_path = tmp_path / "fq_parent_dir"
    fq_parent_dir_path.mkdir()
    (fq_parent_dir_path / "sample1").mkdir()
    (fq_parent_dir_path / "sample2").mkdir()
    (fq_parent_dir_path / "sample3").mkdir()

    sample_sheet_df = watermelon_init.generate_samplesheet([fq_parent_dir_path])
    assert(sample_sheet_df.shape == (3,2))
    assert(list(sample_sheet_df.columns) == ['sample', 'input_dir'])
    assert(list(sample_sheet_df['sample']) == ['sample1', 'sample2', 'sample3'])
    assert(sum(sample_sheet_df['input_dir'].str.contains(str(fq_parent_dir_path))) == 3)
    assert(list(sample_sheet_df['input_dir'].apply(os.path.basename)) == list(sample_sheet_df['sample']))

def test_generate_samplesheet_UsesNaturalSorting(tmp_path):
    fq_parent_dir_path = tmp_path / "fq_parent_dir"
    fq_parent_dir_path.mkdir()
    (fq_parent_dir_path / "sample1").mkdir()
    (fq_parent_dir_path / "sample10").mkdir()
    (fq_parent_dir_path / "sample2").mkdir()
    (fq_parent_dir_path / "sample20").mkdir()
    (fq_parent_dir_path / "sample3").mkdir()
    (fq_parent_dir_path / "sample30").mkdir()

    sample_sheet_df = watermelon_init.generate_samplesheet([fq_parent_dir_path])
    assert(list(sample_sheet_df['sample']) == ['sample1', 'sample2', 'sample3', 'sample10', 'sample20', 'sample30'])


def test_generate_samplesheet_MissingParentDirRaisesRuntimeError(tmp_path):
    fq_parent_dir_path = tmp_path / "fq_parent_dir"
    fq_parent_dir = str(fq_parent_dir_path)
    with pytest.raises(RuntimeError, match=fq_parent_dir + " does not exist") as e_info:
        watermelon_init.generate_samplesheet([fq_parent_dir])

def test_generate_samplesheet_MissingOneParentDirRaisesRuntimeError(tmp_path):
    fq_parent_dir_path1 = tmp_path / "fq_parent_dir1"
    fq_parent_dir_path1.mkdir()
    (fq_parent_dir_path1 / "subdir").mkdir()
    fq_parent_dir_path2 = tmp_path / "fq_parent_dir2"
    with pytest.raises(RuntimeError, match=str(fq_parent_dir_path2) + " does not exist") as e_info:
        watermelon_init.generate_samplesheet([fq_parent_dir_path1, fq_parent_dir_path2])

def test_generate_samplesheet_ParentDirMissingSubdirsRaisesRuntimeError(tmp_path):
    fq_parent_dir_path = tmp_path / "fq_parent_dir"
    fq_parent_dir_path.mkdir()
    fq_parent_dir = str(fq_parent_dir_path)
    with pytest.raises(RuntimeError, match=fq_parent_dir + " contains no subdirectories") as e_info:
        watermelon_init.generate_samplesheet([fq_parent_dir])

def test_generate_samplesheet_OneParentDirMissingSubdirsRaisesRuntimeError(tmp_path):
    fq_parent_dir_path1 = tmp_path / "fq_parent_dir1"
    fq_parent_dir_path1.mkdir()
    (fq_parent_dir_path1 / "subdir").mkdir()
    fq_parent_dir_path2 = tmp_path / "fq_parent_dir2"
    fq_parent_dir_path2.mkdir()

    with pytest.raises(RuntimeError, match=str(fq_parent_dir_path2) + " contains no subdirectories") as e_info:
        watermelon_init.generate_samplesheet([fq_parent_dir_path1, fq_parent_dir_path2])

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


def test_get_template_ExplicitAlternateTemplate(tmp_path):
    # Make alternate template
    yaml_content=\
'''foo:
    foo1: test1
    foo2: test2
'''
    yaml_filepath = tmp_path / 'template.yaml'
    yaml_filepath.write_text(yaml_content)

    args = Namespace(x_alt_template=yaml_filepath)
    with pytest.warns(UserWarning, match='Ignoring type'):
        actual = watermelon_init.get_template(
            args=args,
            pipe_root=watermelon_init._WATERMELON_ROOT
            )

    expected = {'foo': {'foo1': 'test1', 'foo2':'test2'}}
    assert(expected == actual)

def test_get_template_ExplicitAlternateTemplateRaisesIfMissing(tmp_path):
    yaml_filepath = tmp_path / 'template.yaml'

    args = Namespace(x_alt_template=yaml_filepath)
    with pytest.warns(None), pytest.raises(RuntimeError, match='Could not read'):
        watermelon_init.get_template(
            args=args,
            pipe_root=watermelon_init._WATERMELON_ROOT
            )

def test_get_template_InferFilenameFromtype(tmp_path):
    # Make template
    yaml_content=\
'''foo:
    foo1: test1
    foo2: test2
'''
    config_dir = tmp_path / 'config'
    config_dir.mkdir()
    yaml_filepath = config_dir / 'template_foo.yaml'
    yaml_filepath.write_text(yaml_content)

    args = Namespace(type='foo', x_alt_template=None)
    actual = watermelon_init.get_template(
        args=args,
        pipe_root=tmp_path
        )

    expected = {'foo': {'foo1': 'test1', 'foo2':'test2'}}
    assert(expected == actual)


def test_set_up_dirs_Diffex():
    actual = watermelon_init._set_up_dirs(type='diffex', project_id='ABC')
    expected = [
        ('diffex_results', 'analysis_ABC/diffex_results'),
        ('deliverables', 'analysis_ABC/deliverables'),
        ('report', 'analysis_ABC/report')
        ]
    assert(expected == [(k,v) for k,v in actual.items()])

def test_set_up_dirs_AlignQc():
    actual = watermelon_init._set_up_dirs(type='align_qc', project_id='ABC')
    expected = [
        ('alignment_results', 'analysis_ABC/alignment_results'),
        ('deliverables', 'analysis_ABC/deliverables'),
        ('report', 'analysis_ABC/report')
        ]
    assert(expected == [(k,v) for k,v in actual.items()])

def test_set_up_dirs_UnknownTypeReturnsEmptyDict():
    actual = watermelon_init._set_up_dirs(type='foo', project_id='ABC')
    assert(0 == len(actual))


def test_validate_genomes():
    ref_dict = watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, "GRCh38", "align_qc").get("references")
    # Base case - no exceptions should be raised (see the assert False in the except)
    raised = False
    try:
        watermelon_init.validate_genomes(ref_dict)
    except Exception:
        raised = True
    assert(raised == False)
    # Invalid case - fake filepath
    ref_dict["fasta"] = "foo"
    raised = False
    try:
        watermelon_init.validate_genomes(ref_dict) # This time it will raise an exception
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
