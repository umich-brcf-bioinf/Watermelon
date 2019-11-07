#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import os
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

class PhenotypeManagerTest(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_phenotype_sample_list(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            #Set up CSV lines to write
            lines = [
                'sample,A,B,C',
                's3,a1,b1,c1',
                's4,a2,b2,c2',
                's1,a1,b1,c1',
                's2,a2,b2,c2'
            ]
            lines = "\n".join(lines)
            #Create test samplesheet in this tempdir
            test_samplesheet_file = os.path.join(temp_dir_path, 'samplesheet.csv')
            with open(test_samplesheet_file, 'w') as sample_sheet:
                sample_sheet.writelines(lines)
            #Feed it to the PhenotypeManager
            config = {'sample_description_file' : test_samplesheet_file}
            manager = rnaseq_snakefile_helper.PhenotypeManager(config)

        actual_dict = manager.phenotype_sample_list
        actual_dict = testing_utils.ddict2dict(actual_dict) # phenotype_sample_list is built with nested defaultdict. Convert this to builtin dict type
        expected_dict = {'A' : {'a1': ['s3', 's1'], 'a2': ['s4','s2']},
                         'B' : {'b1': ['s3', 's1'], 'b2': ['s4','s2']},
                         'C' : {'c1': ['s3', 's1'], 'c2': ['s4','s2']},}
        self.assertEqual(expected_dict, actual_dict)

    def test_phenotype_sample_list_missingValues(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            #Set up CSV lines to write
            lines = [
                'sample,A,B,C',
                's3,a1,b1,c1',
                's4,,b2,c2',
                's1,a1,,c1',
                's2,a2,,'
            ]
            lines = "\n".join(lines)
            #Create test samplesheet in this tempdir
            test_samplesheet_file = os.path.join(temp_dir_path, 'samplesheet.csv')
            with open(test_samplesheet_file, 'w') as sample_sheet:
                sample_sheet.writelines(lines)
            #Feed it to the PhenotypeManager
            config = {'sample_description_file' : test_samplesheet_file}
            manager = rnaseq_snakefile_helper.PhenotypeManager(config)

        actual_dict = manager.phenotype_sample_list
        actual_dict = testing_utils.ddict2dict(actual_dict) # phenotype_sample_list is built with nested defaultdict. Convert this to builtin dict type
        expected_dict = {'A' : {'a1': ['s3', 's1'], 'a2': ['s2']},
                         'B' : {'b1': ['s3'], 'b2': ['s4']},
                         'C' : {'c1': ['s3', 's1'], 'c2': ['s4']},}
        self.assertEqual(expected_dict, actual_dict)


    def test_phenotype_sample_list_extraPhenotypeValues(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            #Set up CSV lines to write
            lines = [
                'sample,A,B',
                's1,a1,b1',
                's2,a2,b2,c2'
            ]
            lines = "\n".join(lines)
            #Create test samplesheet in this tempdir
            test_samplesheet_file = os.path.join(temp_dir_path, 'samplesheet.csv')
            with open(test_samplesheet_file, 'w') as sample_sheet:
                sample_sheet.writelines(lines)
            #Feed it to the PhenotypeManager
            config = {'sample_description_file' : test_samplesheet_file}

            self.assertRaises(ParserError,
                        rnaseq_snakefile_helper.PhenotypeManager,
                        config)

    def test_phenotype_with_replicates(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            #Set up CSV lines to write
            lines = [
                'sample,A,D,C,B',
                's1,a1,d1,c1,b1',
                's2,a2,d1,c2,b2',
                's3,a3,d1,,b2',
                's4,a4,d2,,'
            ]
            lines = "\n".join(lines)
            #Create test samplesheet in this tempdir
            test_samplesheet_file = os.path.join(temp_dir_path, 'samplesheet.csv')
            with open(test_samplesheet_file, 'w') as sample_sheet:
                sample_sheet.writelines(lines)
            #Feed it to the PhenotypeManager
            config = {'sample_description_file' : test_samplesheet_file}
            manager = rnaseq_snakefile_helper.PhenotypeManager(config)

        actual = manager.phenotypes_with_replicates
        self.assertEqual(['B', 'D'], actual)

class RnaseqSnakefileHelperTest(unittest.TestCase):

    def test_get_sample_reads(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
            rnaseq_snakefile_helper._mkdir(sample_0_dir)
            temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
            sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
            rnaseq_snakefile_helper._mkdir(sample_1_dir)
            temp_dir.write(sample_1_dir + '02_R1.fastq.gz', b'foo')
            temp_dir.write(sample_1_dir + '02_R2.fastq.gz', b'foo')
            samples = ['sample_0', 'sample_1']

            actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
        self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R1_PE', 'R2_PE']},
                         actual_sample_reads)

    def test_get_sample_reads_onlyConsiderR1R2(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
            rnaseq_snakefile_helper._mkdir(sample_0_dir)
            temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
            temp_dir.write(sample_0_dir + '01_R5.fastq.gz', b'foo')
            sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
            rnaseq_snakefile_helper._mkdir(sample_1_dir)
            temp_dir.write(sample_1_dir + '02_R1.fastq.gz', b'foo')
            temp_dir.write(sample_1_dir + '02_R9.fastq.gz', b'foo')
            samples = ['sample_0', 'sample_1']

            actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
        self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R1_SE']},
                         actual_sample_reads)

    def test_get_sample_reads_okIfR2PresentButR1Missing(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
            rnaseq_snakefile_helper._mkdir(sample_0_dir)
            temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
            sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
            rnaseq_snakefile_helper._mkdir(sample_1_dir)
            temp_dir.write(sample_1_dir + '02_R2.fastq.gz', b'foo')
            samples = ['sample_0', 'sample_1']

            actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
        self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R2_SE']},
                         actual_sample_reads)

    def test_get_sample_reads_considersFastqOrFastqGz(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
            rnaseq_snakefile_helper._mkdir(sample_0_dir)
            temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
            sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
            rnaseq_snakefile_helper._mkdir(sample_1_dir)
            temp_dir.write(sample_1_dir + '02_R1.fastq', b'foo')
            samples = ['sample_0', 'sample_1']

            actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
        self.assertEqual({'sample_0': ['R1_SE'], 'sample_1': ['R1_SE']},
                         actual_sample_reads)

    def test_get_sample_reads_missingReads(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
            rnaseq_snakefile_helper._mkdir(sample_0_dir)
            temp_dir.write(sample_0_dir + '01_RX.fastq.gz', b'foo')
            samples = ['sample_0']

            actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
        self.assertEqual({'sample_0': []},
                         actual_sample_reads)

    def test_get_sample_reads_ignoresFilesThatDoNotEndWithFastq(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
            rnaseq_snakefile_helper._mkdir(sample_0_dir)
            temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
            temp_dir.write(sample_0_dir + '01_R2.foo', b'foo')
            temp_dir.write(sample_0_dir + '01_R2.fasta', b'foo')
            temp_dir.write(sample_0_dir + '01_R2', b'foo')
            temp_dir.write(sample_0_dir + 'fastq.01_R2', b'foo')
            samples = ['sample_0']

            actual_sample_reads = rnaseq_snakefile_helper._get_sample_reads(temp_dir_path, samples)
        self.assertEqual({'sample_0': ['R1_SE']},
                         actual_sample_reads)

    def test_flattened_sample_reads(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_0_dir = os.path.join(temp_dir_path, 'sample_0', '')
            rnaseq_snakefile_helper._mkdir(sample_0_dir)
            temp_dir.write(sample_0_dir + '01_R1.fastq.gz', b'foo')
            sample_1_dir = os.path.join(temp_dir_path,'sample_1', '')
            rnaseq_snakefile_helper._mkdir(sample_1_dir)
            temp_dir.write(sample_1_dir + '02_R1.fastq.gz', b'foo')
            temp_dir.write(sample_1_dir + '02_R2.fastq.gz', b'foo')
            samples = ['sample_0', 'sample_1']

            actual_sample_reads = rnaseq_snakefile_helper.flattened_sample_reads(temp_dir_path, samples)
        self.assertEqual([('sample_0', 'R1_SE'), ('sample_1', 'R1_PE'), ('sample_1', 'R2_PE')],
                         actual_sample_reads)

    def test_expand_sample_read_endedness_all(self):
        sample_read_endedness_format = 's={sample}|re={read_endedness}'
        all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
                                      ('s2', 'R1_SE'),
                                      ('s3', 'R1_PE'), ('s3', 'R2_PE')]
        actual = rnaseq_snakefile_helper.expand_sample_read_endedness(sample_read_endedness_format,
                                                                      all_flattened_sample_reads)
        expected = ['s=s1|re=R1_PE', 's=s1|re=R2_PE',
                    's=s2|re=R1_SE',
                    's=s3|re=R1_PE', 's=s3|re=R2_PE']
        self.assertEqual(expected, actual)

    def test_expand_sample_read_endedness_singleSample(self):
        sample_read_endedness_format = 's={sample}|re={read_endedness}'
        all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
                                      ('s2', 'R1_SE'),
                                      ('s3', 'R1_PE'), ('s3', 'R2_PE')]
        actual = rnaseq_snakefile_helper.expand_sample_read_endedness(sample_read_endedness_format,
                                                                      all_flattened_sample_reads,
                                                                      's3')
        expected = ['s=s3|re=R1_PE', 's=s3|re=R2_PE']
        self.assertEqual(expected, actual)

    def test_expand_sample_read_endedness_none(self):
        sample_read_endedness_format = 's={sample}|re={read_endedness}'
        all_flattened_sample_reads = []
        actual = rnaseq_snakefile_helper.expand_sample_read_endedness(sample_read_endedness_format,
                                                                      all_flattened_sample_reads)
        self.assertEqual([], actual)

    def test_expand_read_stats_if_paired(self):
        filename_format = 's={sample}'
        all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
                                      ('s2', 'R1_SE')]
        sample = 's1'
        actual = rnaseq_snakefile_helper.expand_read_stats_if_paired(\
                filename_format,
                all_flattened_sample_reads,
                sample)
        self.assertEqual(['s=s1'], actual)

    def test_expand_read_stats_if_paired_emptyIfSingle(self):
        filename_format = 's={sample}'
        all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
                                      ('s2', 'R1_SE')]
        sample = 's2'
        actual = rnaseq_snakefile_helper.expand_read_stats_if_paired(\
                filename_format,
                all_flattened_sample_reads,
                sample)
        self.assertEqual([], actual)

    def test_expand_read_stats_if_paired_emptyIfMissing(self):
        filename_format = 's={sample}'
        all_flattened_sample_reads = [('s1', 'R1_PE'), ('s1', 'R2_PE'),
                                      ('s2', 'R1_SE')]
        sample = 's3'
        actual = rnaseq_snakefile_helper.expand_read_stats_if_paired(\
                filename_format,
                all_flattened_sample_reads,
                sample)
        self.assertEqual([], actual)
