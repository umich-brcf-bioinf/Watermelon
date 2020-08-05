#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import io
import mock
import os
import sys
import time
import unittest
import yaml
from pandas.errors import ParserError

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

class PhenotypeManagerTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_phenotype_sample_list(self):
        #Set up CSV lines to write
        lines = [
            'sample,A,B,C',
            's3,a1,b1,c1',
            's4,a2,b2,c2',
            's1,a1,b1,c1',
            's2,a2,b2,c2'
        ]
        lines = "\n".join(lines)
        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
        test_samplesheet = io.StringIO(lines)
        #Feed it to the PhenotypeManager
        config = {'samplesheet' : test_samplesheet}
        manager = rnaseq_snakefile_helper.PhenotypeManager(config)

        actual_dict = manager.phenotype_sample_list
        actual_dict = testing_utils.ddict2dict(actual_dict) # phenotype_sample_list is built with nested defaultdict. Convert this to builtin dict type
        expected_dict = {'A' : {'a1': ['s3', 's1'], 'a2': ['s4','s2']},
                         'B' : {'b1': ['s3', 's1'], 'b2': ['s4','s2']},
                         'C' : {'c1': ['s3', 's1'], 'c2': ['s4','s2']},}
        self.assertEqual(expected_dict, actual_dict)

    def test_phenotype_sample_list_missingValues(self):
        #Set up CSV lines to write
        lines = [
            'sample,A,B,C',
            's3,a1,b1,c1',
            's4,,b2,c2',
            's1,a1,,c1',
            's2,a2,,'
        ]
        lines = "\n".join(lines)
        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
        test_samplesheet = io.StringIO(lines)
        #Feed it to the PhenotypeManager
        config = {'samplesheet' : test_samplesheet}
        manager = rnaseq_snakefile_helper.PhenotypeManager(config)

        actual_dict = manager.phenotype_sample_list
        # phenotype_sample_list is built with nested defaultdict. Convert this to builtin dict type
        actual_dict = testing_utils.ddict2dict(actual_dict)
        expected_dict = {'A' : {'a1': ['s3', 's1'], 'a2': ['s2']},
                         'B' : {'b1': ['s3'], 'b2': ['s4']},
                         'C' : {'c1': ['s3', 's1'], 'c2': ['s4']},}
        self.assertEqual(expected_dict, actual_dict)


    def test_phenotype_sample_list_extraPhenotypeValues(self):
        #Set up CSV lines to write
        lines = [
            'sample,A,B',
            's1,a1,b1',
            's2,a2,b2,c2'
        ]
        lines = "\n".join(lines)
        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
        test_samplesheet = io.StringIO(lines)
        #Feed it to the PhenotypeManager, pandas should raise ParserError
        config = {'samplesheet' : test_samplesheet}
        self.assertRaises(ParserError,
                    rnaseq_snakefile_helper.PhenotypeManager,
                    config)

    def test_phenotype_with_replicates(self):
        #Set up CSV lines to write
        lines = [
            'sample,A,D,C,B',
            's1,a1,d1,c1,b1',
            's2,a2,d1,c2,b2',
            's3,a3,d1,,b2',
            's4,a4,d2,,'
        ]
        lines = "\n".join(lines)
        #Use stringIO to create in-memory stream from lines - this is file-like-object readable by pandas
        test_samplesheet = io.StringIO(lines)
        #Feed it to the PhenotypeManager
        config = {'samplesheet' : test_samplesheet}
        manager = rnaseq_snakefile_helper.PhenotypeManager(config)

        actual = manager.phenotypes_with_replicates
        self.assertEqual(['B', 'D'], actual)

class InputFastqManagerTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()
        self.original_stderr = sys.stderr
        sys.stderr = DevNull()

    def tearDown(self):
        os.chdir(self.original_wd)
        sys.stderr = self.original_stderr

    def test_input_paths_dict_fastq(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_1_dir = os.path.join(temp_dir_path, 'sample_1')
            os.mkdir(sample_1_dir)
            sample_1_expected = [os.path.join(sample_1_dir, 'sample_1_R1.fastq.gz')]
            for file in sample_1_expected:
                temp_dir.write(file, b'foo')

            sample_2_dir = os.path.join(temp_dir_path,'sample_2')
            os.mkdir(sample_2_dir)
            sample_2_expected = [os.path.join(sample_2_dir, "sample_2_R{}.fastq.gz".format(x)) for x in [1,2] ]
            for file in sample_2_expected:
                temp_dir.write(file, b'foo')

            #sample 3 is not gzipped, other samples are
            sample_3_dir = os.path.join(temp_dir_path, 'sample_3')
            os.mkdir(sample_3_dir)
            sample_3_expected = [os.path.join(sample_3_dir, "sample_3_R{}.fastq".format(x)) for x in [1,2] ]
            for file in sample_3_expected:
                temp_dir.write(file, b'foo')

            manager = rnaseq_snakefile_helper.InputFastqManager(input_dir = temp_dir_path, capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')
            actual_dict = manager.input_paths_dict

            expected_dict = {
                'sample_1' : sample_1_expected,
                'sample_2' : sample_2_expected,
                'sample_3' : sample_3_expected
            }

            self.assertEqual(expected_dict, actual_dict)

    def test_input_paths_dict_fastq_warns_sample_no_fastqs(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_1_dir = os.path.join(temp_dir_path, 'sample_1')
            os.mkdir(sample_1_dir)
            temp_dir.write(os.path.join(temp_dir_path, 'non_fastq.txt'), b'foo')

            sample_2_dir = os.path.join(temp_dir_path,'sample_2')
            os.mkdir(sample_2_dir)
            sample_2_expected = [os.path.join(sample_2_dir, "sample_2_R{}.fastq.gz".format(x)) for x in [1,2] ]
            for file in sample_2_expected:
                temp_dir.write(file, b'foo')

            self.assertWarns(Warning, rnaseq_snakefile_helper.InputFastqManager, input_dir=temp_dir_path, capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')

    def test_input_paths_dict_fastq_errors_no_fastqs_any_sample(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_1_dir = os.path.join(temp_dir_path, 'sample_1')
            os.mkdir(sample_1_dir)
            temp_dir.write(os.path.join(temp_dir_path, 'non_fastq.txt'), b'foo')

            self.assertRaises(RuntimeError, rnaseq_snakefile_helper.InputFastqManager, input_dir=temp_dir_path, capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')

    def test_input_paths_dict_fastq_cant_mix_gz_plaintext(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_1_dir = os.path.join(temp_dir_path, 'sample_1')
            os.mkdir(sample_1_dir)
            sample_1_file1 = os.path.join(sample_1_dir, 'sample_1_R1.fastq.gz')
            sample_1_file2 = os.path.join(sample_1_dir, 'sample_1_R2.fastq')
            for file in [sample_1_file1, sample_1_file2]:
                temp_dir.write(file, b'foo')

            self.assertRaises(RuntimeError, rnaseq_snakefile_helper.InputFastqManager, input_dir=temp_dir_path, capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')

    def test_fastqs_to_concat_dict(self):
        mock_input_paths_dict = {
            'sample_1': [
                '/path/to/sample_1_L001_R1.fastq.gz',
                '/path/to/sample_1_L002_R1.fastq.gz'
            ],
            'sample_2': [
                '/path/to/sample_2_R1.fastq.gz',
                '/path/to/sample_2_R2.fastq.gz'
            ],
            'sample_3': [
                '/path/to/sample_3_L001_R1.fastq',
                '/path/to/sample_3_L002_R1.fastq',
                '/path/to/sample_3_L001_R2.fastq',
                '/path/to/sample_3_L002_R2.fastq'
            ]
        }
        def __init__(self, input_dir, capture_regex):
            self.input_dir = input_dir
            self.capture_regex = capture_regex
            self.input_paths_dict = mock_input_paths_dict

            # Specifically testing the _fastqs_to_concat_from_filenames method here
            self.fastqs_to_concat_dict = self._fastqs_to_concat_from_filenames(capture_regex=r'.*_R(\d+)\.fastq.*')

        expected = {
            'sample_1': {'1': ['/path/to/sample_1_L001_R1.fastq.gz', '/path/to/sample_1_L002_R1.fastq.gz']},
            'sample_2': {'1': ['/path/to/sample_2_R1.fastq.gz'], '2': ['/path/to/sample_2_R2.fastq.gz']},
            'sample_3': {'1': ['/path/to/sample_3_L001_R1.fastq', '/path/to/sample_3_L002_R1.fastq'], '2': ['/path/to/sample_3_L001_R2.fastq', '/path/to/sample_3_L002_R2.fastq']}
        }

        with mock.patch.object(rnaseq_snakefile_helper.InputFastqManager, '__init__', __init__):
            manager = rnaseq_snakefile_helper.InputFastqManager(input_dir='foo', capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')
            actual = manager.fastqs_to_concat_dict
            self.assertEqual(expected, actual)

    def test_sample_bnames_dict(self):
        mock_input_paths_dict = {
            'sample_1': [
                '/path/to/sample_1_L001_R1.fastq.gz',
                '/path/to/sample_1_L002_R1.fastq.gz'
            ],
            'sample_2': [
                '/path/to/sample_2_R1.fastq.gz',
                '/path/to/sample_2_R2.fastq.gz'
            ],
            'sample_3': [
                '/path/to/sample_3_L001_R1.fastq',
                '/path/to/sample_3_L002_R1.fastq',
                '/path/to/sample_3_L001_R2.fastq',
                '/path/to/sample_3_L002_R2.fastq'
            ],
            'sample_4': [
                '/path/to/foo_blahblah_R1.fastq.gz',
                '/path/to/foo_barbazbi_R1.fastq.gz',
                '/path/to/whateverever_R2.fastq.gz',
                '/path/to/doesntmatter_R2.fastq.gz'
            ]
        }
        def __init__(self, input_dir, capture_regex):
            self.input_dir = input_dir
            self.capture_regex = capture_regex
            self.input_paths_dict = mock_input_paths_dict
            self.sample_bnames_dict = self._sample_bnames_from_filenames(
                                            capture_regex=r'.*_R(\d+)\.fastq.*',
                                            bname_fmt='{}_R{}')

        expected = {
            'sample_1': ['sample_1_R1'],
            'sample_2': ['sample_2_R1', 'sample_2_R2'],
            'sample_3': ['sample_3_R1', 'sample_3_R2'],
            'sample_4': ['sample_4_R1', 'sample_4_R2']
        }

        with mock.patch.object(rnaseq_snakefile_helper.InputFastqManager, '__init__', __init__):
            manager = rnaseq_snakefile_helper.InputFastqManager(input_dir='foo', capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')
            actual = manager.sample_bnames_dict
            self.assertEqual(expected, actual)


    def test_gather_basenames(self):
        mock_sample_bnames_dict = {
            'sample_1': ['sample_1_R1'],
            'sample_2': ['sample_2_R1', 'sample_2_R2'],
            'sample_3': ['sample_3_R1', 'sample_3_R2']
        }
        def __init__(self, input_dir, capture_regex):
            self.input_dir = input_dir
            self.capture_regex = capture_regex
            self.sample_bnames_dict = mock_sample_bnames_dict

        expected = ['sample_1_R1', 'sample_2_R1', 'sample_2_R2', 'sample_3_R1', 'sample_3_R2']

        with mock.patch.object(rnaseq_snakefile_helper.InputFastqManager, '__init__', __init__):
            manager = rnaseq_snakefile_helper.InputFastqManager(input_dir='foo', capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')
            sample_list = ['sample_1', 'sample_2', 'sample_3']
            actual = manager.gather_basenames(sample_list)
            self.assertEqual(expected, actual)

    def test_gather_basenames_error_if_missing_input(self):
        mock_sample_bnames_dict = {
            'sample_1': ['sample_1_R1.fastq.gz']
        }
        def __init__(self, input_dir, capture_regex):
            self.input_dir = input_dir
            self.capture_regex = capture_regex
            self.sample_bnames_dict = mock_sample_bnames_dict

        with mock.patch.object(rnaseq_snakefile_helper.InputFastqManager, '__init__', __init__):
            manager = rnaseq_snakefile_helper.InputFastqManager(input_dir='foo', capture_regex=r'.*_R(\d+)[_0-9]*\.fastq.*')
            sample_list = ['sample_1', 'sample_2', 'sample_3']
            self.assertRaises(RuntimeError, manager.gather_basenames, sample_list)


# FIXME: Revive any of these?
# class RnaseqSnakefileHelperTest(unittest.TestCase):
#
#     def test_get_sample_reads(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
#             rnaseq_snakefile_helper._mkdir(sample_0_dir)
#             temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
#             sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
#             rnaseq_snakefile_helper._mkdir(sample_1_dir)
#             temp_dir.write(sample_1_dir + '02_R1.fastq.gz', b'foo')
#             temp_dir.write(sample_1_dir + '02_R2.fastq.gz', b'foo')
#             samples = ['sample_0', 'sample_1']
#
#             actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
#         self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R1_PE', 'R2_PE']},
#                          actual_sample_reads)
#
#     def test_get_sample_reads_onlyConsiderR1R2(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
#             rnaseq_snakefile_helper._mkdir(sample_0_dir)
#             temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
#             temp_dir.write(sample_0_dir + '01_R5.fastq.gz', b'foo')
#             sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
#             rnaseq_snakefile_helper._mkdir(sample_1_dir)
#             temp_dir.write(sample_1_dir + '02_R1.fastq.gz', b'foo')
#             temp_dir.write(sample_1_dir + '02_R9.fastq.gz', b'foo')
#             samples = ['sample_0', 'sample_1']
#
#             actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
#         self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R1_SE']},
#                          actual_sample_reads)
#
#     def test_get_sample_reads_okIfR2PresentButR1Missing(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
#             rnaseq_snakefile_helper._mkdir(sample_0_dir)
#             temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
#             sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
#             rnaseq_snakefile_helper._mkdir(sample_1_dir)
#             temp_dir.write(sample_1_dir + '02_R2.fastq.gz', b'foo')
#             samples = ['sample_0', 'sample_1']
#
#             actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
#         self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R2_SE']},
#                          actual_sample_reads)
#
#     def test_get_sample_reads_considersFastqOrFastqGz(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
#             rnaseq_snakefile_helper._mkdir(sample_0_dir)
#             temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
#             sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
#             rnaseq_snakefile_helper._mkdir(sample_1_dir)
#             temp_dir.write(sample_1_dir + '02_R1.fastq', b'foo')
#             samples = ['sample_0', 'sample_1']
#
#             actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
#         self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R1_SE']},
#                          actual_sample_reads)
#
#     def test_get_sample_reads_missingReads(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
#             rnaseq_snakefile_helper._mkdir(sample_0_dir)
#             temp_dir.write(sample_0_dir + '01_RX.fastq.gz', b'foo')
#             samples = ['sample_0']
#
#             actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
#         self.assertEqual({'sample_0': []},
#                          actual_sample_reads)
#
#     def test_get_sample_reads_ignoresFilesThatDoNotEndWithFastq(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
#             rnaseq_snakefile_helper._mkdir(sample_0_dir)
#             temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
#             temp_dir.write(sample_0_dir + '01_R2.foo', b'foo')
#             temp_dir.write(sample_0_dir + '01_R2.fasta', b'foo')
#             temp_dir.write(sample_0_dir + '01_R2', b'foo')
#             temp_dir.write(sample_0_dir + 'fastq.01_R2', b'foo')
#             samples = ['sample_0']
#
#             actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
#         self.assertEqual({'sample_0': ['R1_SE']},
#                          actual_sample_reads)
#
#     def test_flattened_sample_reads(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
#             rnaseq_snakefile_helper._mkdir(sample_0_dir)
#             temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
#             sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
#             rnaseq_snakefile_helper._mkdir(sample_1_dir)
#             temp_dir.write(sample_1_dir + '02_R1.fastq.gz', b'foo')
#             temp_dir.write(sample_1_dir + '02_R2.fastq.gz', b'foo')
#             samples = ['sample_0', 'sample_1']
#
#             actual_sample_reads = rnaseq_snakefile_helper.flattened_sample_reads(temp_dir_path, samples)
#         self.assertEqual([('sample_0', 'R1_SE'), ('sample_1', 'R1_PE'), ('sample_1', 'R2_PE')],
#                          actual_sample_reads)
#
#     def test_expand_sample_read_endedness_all(self):
#         sample_read_endedness_format = 's={sample}|re={read_endedness}'
#         all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
#                                       ('s2', 'R1_SE'),
#                                       ('s3', 'R1_PE'), ('s3', 'R2_PE')]
#         actual = rnaseq_snakefile_helper.expand_sample_read_endedness(sample_read_endedness_format,
#                                                                       all_flattened_sample_reads)
#         expected = ['s=s1|re=R1_PE', 's=s1|re=R2_PE',
#                     's=s2|re=R1_SE',
#                     's=s3|re=R1_PE', 's=s3|re=R2_PE']
#         self.assertEqual(expected, actual)
#
#     def test_expand_sample_read_endedness_singleSample(self):
#         sample_read_endedness_format = 's={sample}|re={read_endedness}'
#         all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
#                                       ('s2', 'R1_SE'),
#                                       ('s3', 'R1_PE'), ('s3', 'R2_PE')]
#         actual = rnaseq_snakefile_helper.expand_sample_read_endedness(sample_read_endedness_format,
#                                                                       all_flattened_sample_reads,
#                                                                       's3')
#         expected = ['s=s3|re=R1_PE', 's=s3|re=R2_PE']
#         self.assertEqual(expected, actual)
#
#     def test_expand_sample_read_endedness_none(self):
#         sample_read_endedness_format = 's={sample}|re={read_endedness}'
#         all_flattened_sample_reads = []
#         actual = rnaseq_snakefile_helper.expand_sample_read_endedness(sample_read_endedness_format,
#                                                                       all_flattened_sample_reads)
#         self.assertEqual([], actual)
#
#     def test_expand_read_stats_if_paired(self):
#         filename_format = 's={sample}'
#         all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
#                                       ('s2', 'R1_SE')]
#         sample = 's1'
#         actual = rnaseq_snakefile_helper.expand_read_stats_if_paired(\
#                 filename_format,
#                 all_flattened_sample_reads,
#                 sample)
#         self.assertEqual(['s=s1'], actual)
#
#     def test_expand_read_stats_if_paired_emptyIfSingle(self):
#         filename_format = 's={sample}'
#         all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
#                                       ('s2', 'R1_SE')]
#         sample = 's2'
#         actual = rnaseq_snakefile_helper.expand_read_stats_if_paired(\
#                 filename_format,
#                 all_flattened_sample_reads,
#                 sample)
#         self.assertEqual([], actual)
#
#     def test_expand_read_stats_if_paired_emptyIfMissing(self):
#         filename_format = 's={sample}'
#         all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
#                                       ('s2', 'R1_SE')]
#         sample = 's3'
#         actual = rnaseq_snakefile_helper.expand_read_stats_if_paired(\
#                 filename_format,
#                 all_flattened_sample_reads,
#                 sample)
#         self.assertEqual([], actual)
