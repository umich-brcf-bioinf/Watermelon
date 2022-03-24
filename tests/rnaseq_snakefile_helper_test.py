#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import io
import mock
import os
import pandas as pd
import pytest
import sys
import time
import unittest
import yaml

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures.tempdirectory import TempDirectory

from scripts import rnaseq_snakefile_helper # local module that we're testing
from tests import testing_utils # local module

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
WATERMELON_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))
EXAMPLE_CONFIGFILE_PATH = os.path.join(WATERMELON_BASE_DIR, 'config', 'example_config.yaml')

#Allow suppression of stdout
class DevNull(object):
    def write(self, data):
        pass

#class PhenotypeManagerTest(unittest.TestCase):
#    def setUp(self):
#        self.original_wd = os.getcwd()
#
#    def tearDown(self):
#        os.chdir(self.original_wd)
#
#    def test_phenotype_sample_list(self):
#        #Set up CSV lines to write
#        lines = [
#            'sample,A,B,C',
#            's3,a1,b1,c1',
#            's4,a2,b2,c2',
#            's1,a1,b1,c1',
#            's2,a2,b2,c2'
#        ]
#        lines = "\n".join(lines)
#        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
#        test_samplesheet = io.StringIO(lines)
#        self.sample_phenotype_value_dict = self.samplesheet.to_dict(orient='index') #TWS TODO use this instead of PhenotypeManager
#        #Feed it to the PhenotypeManager
#        config = {'samplesheet' : test_samplesheet}
#        manager = rnaseq_snakefile_helper.PhenotypeManager(config)
#
#        actual_dict = manager.phenotype_sample_list
#        actual_dict = testing_utils.ddict2dict(actual_dict) # phenotype_sample_list is built with nested defaultdict. Convert this to builtin dict type
#        expected_dict = {'A' : {'a1': ['s3', 's1'], 'a2': ['s4','s2']},
#                         'B' : {'b1': ['s3', 's1'], 'b2': ['s4','s2']},
#                         'C' : {'c1': ['s3', 's1'], 'c2': ['s4','s2']},}
#        self.assertEqual(expected_dict, actual_dict)
#
#    def test_phenotype_sample_list_missingValues(self):
#        #Set up CSV lines to write
#        lines = [
#            'sample,A,B,C',
#            's3,a1,b1,c1',
#            's4,,b2,c2',
#            's1,a1,,c1',
#            's2,a2,,'
#        ]
#        lines = "\n".join(lines)
#        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
#        test_samplesheet = io.StringIO(lines)
#        #Feed it to the PhenotypeManager
#        config = {'samplesheet' : test_samplesheet}
#        manager = rnaseq_snakefile_helper.PhenotypeManager(config)
#
#        actual_dict = manager.phenotype_sample_list
#        # phenotype_sample_list is built with nested defaultdict. Convert this to builtin dict type
#        actual_dict = testing_utils.ddict2dict(actual_dict)
#        expected_dict = {'A' : {'a1': ['s3', 's1'], 'a2': ['s2']},
#                         'B' : {'b1': ['s3'], 'b2': ['s4']},
#                         'C' : {'c1': ['s3', 's1'], 'c2': ['s4']},}
#        self.assertEqual(expected_dict, actual_dict)
#
#
#    def test_phenotype_sample_list_extraPhenotypeValues(self):
#        #Set up CSV lines to write
#        lines = [
#            'sample,A,B',
#            's1,a1,b1',
#            's2,a2,b2,c2'
#        ]
#        lines = "\n".join(lines)
#        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
#        test_samplesheet = io.StringIO(lines)
#        #Feed it to the PhenotypeManager, pandas should raise ParserError
#        config = {'samplesheet' : test_samplesheet}
#        self.assertRaises(pd.errors.ParserError,
#                    rnaseq_snakefile_helper.PhenotypeManager,
#                    config)
#
#    def test_phenotype_with_replicates(self):
#        #Set up CSV lines to write
#        lines = [
#            'sample,A,D,C,B',
#            's1,a1,d1,c1,b1',
#            's2,a2,d1,c2,b2',
#            's3,a3,d1,,b2',
#            's4,a4,d2,,'
#        ]
#        lines = "\n".join(lines)
#        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
#        test_samplesheet = io.StringIO(lines)
#        #Feed it to the PhenotypeManager
#        config = {'samplesheet' : test_samplesheet}
#        manager = rnaseq_snakefile_helper.PhenotypeManager(config)
#
#        actual = manager.phenotypes_with_replicates
#        self.assertEqual(['B', 'D'], actual)

# Use this for some pytests
def write_foo_file(fname):
    with open(fname, "w") as fh:
        fh.write('foo')


def test_get_sample_fastq_paths_no_fastqs_returns_none(tmp_path):
    s1_glob = str(tmp_path / "sample_1_*.fastq.gz")
    actual = rnaseq_snakefile_helper.get_sample_fastq_paths(s1_glob)
    assert(actual == None)


def test_get_sample_fastq_paths_cant_mix_gz_plaintext(tmp_path):
    # Create the fake inputs
    s1_glob = str(tmp_path / "sample_1_*.fastq*")
    write_foo_file(str(tmp_path / "sample_1_R1.fastq.gz"))
    write_foo_file(str(tmp_path / "sample_1_R2.fastq"))

    with pytest.raises(RuntimeError, match="contains a mixture of fastq and fastq.gz files"):
        rnaseq_snakefile_helper.get_sample_fastq_paths(s1_glob)


def test_fastqs_to_concat(tmp_path):
    # Create the fake inputs
    s1_glob = tmp_path / "sample_1_*.fastq.gz"
    #s1_dir.mkdir()
    write_foo_file(str(tmp_path / "sample_1_R1.fastq.gz")) # sample_1 only has R1
    s2_glob = tmp_path / "sample_2_*.fastq.gz"
    #s2_dir.mkdir()
    write_foo_file(str(tmp_path / "sample_2_R1.fastq.gz"))
    write_foo_file(str(tmp_path / "sample_2_R2.fastq.gz"))
    s3_glob = tmp_path / "sample_3_*.fastq.gz"
    #s3_dir.mkdir()
    write_foo_file(str(tmp_path / "sample_3_R1_L001.fastq.gz")) # sample_3_R1 has 3 lanes
    write_foo_file(str(tmp_path / "sample_3_R1_L002.fastq.gz"))
    write_foo_file(str(tmp_path / "sample_3_R1_L003.fastq.gz"))
    write_foo_file(str(tmp_path / "sample_3_R2.fastq.gz"))

    sample_sheet_str = """sample,input_glob
sample_1,{}
sample_2,{}
sample_3,{}""".format(s1_glob,s2_glob,s3_glob)
    sample_sheet_df = pd.read_csv(io.StringIO(sample_sheet_str), comment='#', dtype='object') \
        .set_index("sample", drop=True)

    expected = {
        'sample_1': {'1': [str(tmp_path / "sample_1_R1.fastq.gz")]},
        'sample_2': {'1': [str(tmp_path / "sample_2_R1.fastq.gz")], '2': [str(tmp_path / "sample_2_R2.fastq.gz")]},
        'sample_3': {'1': [str(tmp_path / "sample_3_R1_L001.fastq.gz"), str(tmp_path / "sample_3_R1_L002.fastq.gz"), str(tmp_path / "sample_3_R1_L003.fastq.gz")], '2': [str(tmp_path / "sample_3_R2.fastq.gz")]}
    }

    actual = rnaseq_snakefile_helper.fastqs_to_concat(
        samplesheet=sample_sheet_df,
        capture_regex=r'.*_R(\d+)[_L0-9]*\.fastq.*',
        )
    actual = testing_utils.ddict2dict(actual)
    assert(expected == actual)


def test_sample_bnames_from_filenames(tmp_path):
    # Create the fake inputs
    s1_glob = tmp_path / "sample_1_*.fastq.gz"
    #s1_dir.mkdir()
    write_foo_file(str(tmp_path / "sample_1_R1.fastq.gz")) # sample_1 only has R1
    s2_glob = tmp_path / "sample_2_*.fastq.gz"
    #s2_dir.mkdir()
    write_foo_file(str(tmp_path / "sample_2_R1.fastq.gz"))
    write_foo_file(str(tmp_path / "sample_2_R2.fastq.gz"))
    s3_glob = tmp_path / "sample_3_*.fastq.gz"
    #s3_dir.mkdir()
    write_foo_file(str(tmp_path / "sample_3_R1.fastq.gz"))
    write_foo_file(str(tmp_path / "sample_3_R2.fastq.gz"))
    s4_glob = tmp_path / "sample_4_*.fastq.gz"
    #s4_dir.mkdir()
    write_foo_file(str(tmp_path / "sample_4_R1.fastq.gz"))
    write_foo_file(str(tmp_path / "sample_4_R2.fastq.gz"))

    sample_sheet_str = """sample,input_glob
sample_1,{}
sample_2,{}
sample_3,{}
sample_4,{}""".format(s1_glob,s2_glob,s3_glob,s4_glob)
    sample_sheet_df = pd.read_csv(io.StringIO(sample_sheet_str), comment='#', dtype='object') \
        .set_index("sample", drop=True)

    expected = {
        'sample_1': ['sample_1_R1'],
        'sample_2': ['sample_2_R1', 'sample_2_R2'],
        'sample_3': ['sample_3_R1', 'sample_3_R2'],
        'sample_4': ['sample_4_R1', 'sample_4_R2']
    }

    actual = rnaseq_snakefile_helper.sample_bnames_from_filenames(
        samplesheet=sample_sheet_df,
        capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*',
        bname_fmt='{}_R{}'
        )
    assert(expected == actual)


def test_gather_basenames():
    samplename_list = ['sample_1', 'sample_2', 'sample_3']
    sample_bnames_dict = {
        'sample_1': ['sample_1_R1'],
        'sample_2': ['sample_2_R1', 'sample_2_R2'],
        'sample_3': ['sample_3_R1', 'sample_3_R2']
    }
    expected = ['sample_1_R1', 'sample_2_R1', 'sample_2_R2', 'sample_3_R1', 'sample_3_R2']

    actual_bnames = rnaseq_snakefile_helper.gather_basenames(sample_bnames_dict, samplename_list)
    assert (expected == actual_bnames)


def test_gather_basenames_error_if_missing_input():
    sample_bnames_dict = {
        'sample_1': ['sample_1_R1.fastq.gz']
    }
    samplename_list = ['sample_1', 'sample_2', 'sample_3']
    with pytest.raises(RuntimeError):
        rnaseq_snakefile_helper.gather_basenames(sample_bnames_dict, samplename_list)
