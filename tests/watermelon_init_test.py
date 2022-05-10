from argparse import Namespace
import os
import pytest

import pandas as pd

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import watermelon_init

def test_generate_samplesheet_Basecase(tmp_path):
    fq_parent_dir_path = tmp_path / "fq_parent_dir"
    fq_parent_dir_path.mkdir()
    (fq_parent_dir_path / "sample_01_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_01_R2.fastq.gz").touch()
    (fq_parent_dir_path / "sample_02_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_02_R2.fastq.gz").touch()
    (fq_parent_dir_path / "sample_03_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_03_R2.fastq.gz").touch()

    sample_fq_regex = watermelon_init._DEFAULT_SAMPLE_FASTQ_REGEX
    autoglob_ext = watermelon_init._DEFAULT_AUTOGLOB_EXT
    sample_sheet_df = watermelon_init.generate_samplesheet([fq_parent_dir_path], sample_fq_regex, autoglob_ext)
    assert(sample_sheet_df.shape == (3,2))
    assert(list(sample_sheet_df.columns) == ['sample', 'input_glob'])
    assert(list(sample_sheet_df['sample']) == ['sample_01', 'sample_02', 'sample_03'])
    #assert(sum(sample_sheet_df['input_glob'].str.contains(str(fq_parent_dir_path))) == 3)
    #assert(list(sample_sheet_df['input_glob'].apply(os.path.basename)) == list(sample_sheet_df['sample']))

def test_generate_samplesheet_UsesNaturalSorting(tmp_path):
    fq_parent_dir_path = tmp_path / "fq_parent_dir"
    fq_parent_dir_path.mkdir()
    (fq_parent_dir_path / "sample_1_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_1_R2.fastq.gz").touch()
    (fq_parent_dir_path / "sample_10_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_10_R2.fastq.gz").touch()
    (fq_parent_dir_path / "sample_2_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_2_R2.fastq.gz").touch()
    (fq_parent_dir_path / "sample_20_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_20_R2.fastq.gz").touch()
    (fq_parent_dir_path / "sample_3_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_3_R2.fastq.gz").touch()
    (fq_parent_dir_path / "sample_30_R1.fastq.gz").touch()
    (fq_parent_dir_path / "sample_30_R2.fastq.gz").touch()

    sample_fq_regex = watermelon_init._DEFAULT_SAMPLE_FASTQ_REGEX
    autoglob_ext = watermelon_init._DEFAULT_AUTOGLOB_EXT
    sample_sheet_df = watermelon_init.generate_samplesheet([fq_parent_dir_path], sample_fq_regex, autoglob_ext)
    assert(list(sample_sheet_df['sample']) == ['sample_1', 'sample_2', 'sample_3', 'sample_10', 'sample_20', 'sample_30'])


def test_generate_samplesheet_MissingParentDirRaisesRuntimeError(tmp_path):
    fq_parent_dir_path = tmp_path / "fq_parent_dir"
    fq_parent_dir = str(fq_parent_dir_path)
    sample_fq_regex = watermelon_init._DEFAULT_SAMPLE_FASTQ_REGEX
    autoglob_ext = watermelon_init._DEFAULT_AUTOGLOB_EXT
    with pytest.raises(RuntimeError, match=fq_parent_dir + " does not exist") as e_info:
        watermelon_init.generate_samplesheet([fq_parent_dir], sample_fq_regex, autoglob_ext)

def test_generate_samplesheet_MissingOneParentDirRaisesRuntimeError(tmp_path):
    fq_parent_dir_path1 = tmp_path / "fq_parent_dir1"
    fq_parent_dir_path1.mkdir()
    (fq_parent_dir_path1 / "subdir").mkdir()
    fq_parent_dir_path2 = tmp_path / "fq_parent_dir2"
    sample_fq_regex = watermelon_init._DEFAULT_SAMPLE_FASTQ_REGEX
    autoglob_ext = watermelon_init._DEFAULT_AUTOGLOB_EXT
    with pytest.raises(RuntimeError, match=str(fq_parent_dir_path2) + " does not exist") as e_info:
        watermelon_init.generate_samplesheet([fq_parent_dir_path1, fq_parent_dir_path2], sample_fq_regex, autoglob_ext)


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

def test_validate_count_matrix_validMatrix():
    # Base case - no exceptions should be raised (see the assert False in the except)
    counts_str = '''gene_id	sample_01	sample_02
gene_foo	55	0
gene_bar	22	15
gene_baz	0	338'''
    counts_file_like = StringIO(counts_str)
    raised = False
    try:
        watermelon_init.validate_count_matrix(counts_file_like)
    except Exception as foo:
        raised = True
    assert(raised == False)

def test_validate_count_matrix_invalidMatrix():
    # Invalid case - has non-numeric columns
    counts_str = '''gene_id	description	sample_01	sample_02
gene_foo	The amazing foo gene	55	0
gene_bar	The incredible bar gene	22	15
gene_baz	The indisputably awesome baz gene	0	338'''
    counts_file_like = StringIO(counts_str)
    with pytest.raises(RuntimeError, match="The count matrix must only contain"):
        watermelon_init.validate_count_matrix(counts_file_like)

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

def test_validate_fastq_inputs():
    # Assemble an in-memory samplesheet for data contained in this repo
    fastq_basedir = os.path.join(watermelon_init._WATERMELON_ROOT, "data", "sim_reads_human_chr22")
    header = ["sample,input_glob"]
    content = ["sample_0{},{}/sample_0{}_*.fastq.gz".format(x, fastq_basedir, x) for x in range(1,7)]
    samplesheet_str = "\n".join(header + content)
    samplesheet_file_like = StringIO(samplesheet_str)
    ss_df = pd.read_csv(samplesheet_file_like)
    # Base case
    watermelon_init.validate_fastq_inputs(
        ss_df=ss_df,
        sample_col="sample",
        fq_col="input_glob"
    )
    # Test broken glob
    ss_df.at[0, "input_glob"] = "samplefoo.xyz"
    with pytest.raises(RuntimeError, match = "No fastq files found"):
        watermelon_init.validate_fastq_inputs(
            ss_df=ss_df,
            sample_col="sample",
            fq_col="input_glob"
        )


def test_set_up_refs_validKey():
    ref_dict = watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, "GRCh38", "align_qc")
    assert(sorted(ref_dict.keys()) == ["ensembl_version", "fastq_screen", "genome", "references"])


def test_set_up_refs_invalidKey():
    with pytest.raises(RuntimeError) as e_info:
        watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, "foobar", "align_qc")
