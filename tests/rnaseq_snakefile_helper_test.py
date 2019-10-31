#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import os
import time
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures.tempdirectory import TempDirectory
import yaml

import scripts.rnaseq_snakefile_helper as rnaseq_snakefile_helper
from scripts.rnaseq_snakefile_helper import  PhenotypeManager

class PhenotypeManagerTest(unittest.TestCase):
    def test_phenotype_sample_list(self):
        phenotype_labels = 'A | B | C'
        sample_phenotype_value_dict = {'s3' : ' a1 | b1 | c1 ',
                                       's4' : ' a2 | b2 | c2 ',
                                       's1' : ' a1 | b1 | c1 ',
                                       's2' : ' a2 | b2 | c2 ',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                     'samples' : sample_phenotype_value_dict},
                                   delimiter)
        actual_dict = manager.phenotype_sample_list
        expected_dict = {'A' : {'a1': ['s1', 's3'], 'a2': ['s2','s4']},
                         'B' : {'b1': ['s1', 's3'], 'b2': ['s2','s4']},
                         'C' : {'c1': ['s1', 's3'], 'c2': ['s2','s4']},}

        self.assertEqual(expected_dict, actual_dict)

    def test_phenotype_sample_list_missingValues(self):
        phenotype_labels = 'A | B | C'
        sample_phenotype_value_dict = {'s3' : ' a1 | b1 | c1 ',
                                       's4' : '    | b2 | c2 ',
                                       's1' : ' a1 |    | c1 ',
                                       's2' : ' a2 |    | ',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                     'samples' : sample_phenotype_value_dict},
                                   delimiter)
        actual_dict = manager.phenotype_sample_list
        expected_dict = {'A' : {'a1': ['s1', 's3'], 'a2': ['s2']},
                         'B' : {'b1': ['s3'], 'b2': ['s4']},
                         'C' : {'c1': ['s1', 's3'], 'c2': ['s4']},}

        self.assertEqual(expected_dict, actual_dict)

    def test_phenotype_sample_list_missingPhenotypeValues(self):
        phenotype_labels = 'A|B'
        sample_phenotype_value_dict = {'s1' : ' a1 | b1',
                                       's2' : ' a2 ',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                     'samples' : sample_phenotype_value_dict},
                                   delimiter)
        self.assertRaisesRegexp(ValueError,
                                r'expected 2 .* but sample s2 .*1',
                                getattr,
                                manager,
                                'phenotype_sample_list')

    def test_phenotype_sample_list_extraPhenotypeValues(self):
        phenotype_labels = 'A|B'
        sample_phenotype_value_dict = {'s1' : ' a1 | b1',
                                       's2' : ' a2 | b2 | c2',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                     'samples' : sample_phenotype_value_dict},
                                   delimiter)

        self.assertRaisesRegexp(ValueError,
                                r'expected 2 .* but sample s2 .*3',
                                getattr,
                                manager,
                                'phenotype_sample_list')

    def test_phenotype_sample_list_trailingDelimIsExtraPhenotypeValue(self):
        phenotype_labels = 'A|B'
        sample_phenotype_value_dict = {'s1' : ' a1 | b1 |',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                    'samples' : sample_phenotype_value_dict},
                                   delimiter)

        self.assertRaisesRegexp(ValueError,
                                r'expected 2 .* but sample s1 .*3',
                                getattr,
                                manager,
                                'phenotype_sample_list')

    def test_phenotype_sample_list_missingPhenotypeLabel(self):
        phenotype_labels = 'A||C'
        sample_phenotype_value_dict = {'s1' : ' a1 | b1 | c1',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                     'samples' : sample_phenotype_value_dict},
                                   delimiter)
        self.assertRaisesRegexp(ValueError,
                                r'label of phenotype 2 is empty',
                                getattr,
                                manager,
                                'phenotype_sample_list')

    def test_phenotype_sample_list_trailingDelimIsMissingPhenotypeLabel(self):
        phenotype_labels = 'A|B|C|'
        sample_phenotype_value_dict = {'s1' : ' a1 | b1 | c1',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                    'samples' : sample_phenotype_value_dict},
                                   delimiter)

        self.assertRaisesRegexp(ValueError,
                                r'label of phenotype 4 is empty',
                                getattr,
                                manager,
                                'phenotype_sample_list')

    def test_phenotype_with_replicates(self):
        phenotype_labels = 'A | D | C | B'
        sample_dict = {'s1' : ' a1 | d1 | c1 | b1',
                       's2' : ' a2 | d1 | c2 | b2',
                       's3' : ' a3 | d1 |    | b2',
                       's4' : ' a4 | d2 |    |   ',}
        delimiter = '|'
        manager = PhenotypeManager({'phenotypes' : phenotype_labels,
                                     'samples' : sample_dict},
                                   delimiter)
        actual = manager.phenotypes_with_replicates

        self.assertEqual(['B', 'D'], actual)

    def test_separated_comparisons(self):
        comparison_dict = {'phenoLabel1' : ['1A_v_1B','1C_v_1D'],
                           'phenoLabel2' : ['2A_v_2B','2C_v_2D']}
        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.separated_comparisons('#')

        expected_comparisons = {'phenoLabel1' : '1A_v_1B#1C_v_1D',
                                'phenoLabel2' : '2A_v_2B#2C_v_2D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_separated_comparisons_stripsWhitespace(self):
        comparison_dict = {'phenoLabel1' : ['\t1A_v_1B\t','\t1C_v_1D\t'],
                           'phenoLabel2' : ['  2A_v_2B  ','  2C_v_2D  ']}
        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.separated_comparisons('$')

        expected_comparisons = {'phenoLabel1' : '1A_v_1B$1C_v_1D',
                                'phenoLabel2' : '2A_v_2B$2C_v_2D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_separated_comparisons_emptyInputReturnsEmptyDict(self):
        manager = PhenotypeManager({'comparisons' : {}})
        actual_comparisons = manager.separated_comparisons(',')
        self.assertEqual({}, actual_comparisons)

    def test_concatenated_comparison_values(self):
        comparison_dict = {'phenoLabel1' : ['1A_v_1B','1C_v_1D'],
                           'phenoLabel2' : ['2A_v_2B','2C_v_2D']}
        comparison_infix = '_v_'
        delimiter = '|^|'
        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)
        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D',
                                'phenoLabel2' : '2A|^|2B|^|2C|^|2D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_concatenated_comparison_values_stripsWhitespace(self):
        comparison_dict = {'phenoLabel1' : [' \t1A_v_1B \t', '\t 1C_v_1D \t']}
        comparison_infix = '_v_'
        delimiter = '|^|'

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)

        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_concatenated_comparison_values_duplicateValuesCoalesced(self):
        comparison_dict = {'phenoLabel1' : ['1A_v_1B','1A_v_1C','1A_v_1D']}
        comparison_infix = '_v_'
        delimiter = '|^|'

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)

        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D'}
        self.assertEqual(expected_comparisons, actual_comparisons)


    def test_concatenated_comparison_values_reorders(self):
        comparison_dict = {'phenoLabel1' : ['1D_v_1B','1C_v_1A']}
        comparison_infix = '_v_'
        delimiter = '|^|'

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)

        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_phenotypes_comparisons_all_tuple(self):
        comparison_dict = {'phenoLabel1' : ['1A_v_1B','1C_v_1D'],
                           'phenoLabel2' : ['2A_v_2B','2C_v_2D']}

        manager = PhenotypeManager({'comparisons' : comparison_dict})

        actual = manager.phenotypes_comparisons_all_tuple

        self.assertEqual(['phenoLabel1', 'phenoLabel1', 'phenoLabel2', 'phenoLabel2'],
                         actual.phenotypes)
        self.assertEqual(['1A_v_1B', '1C_v_1D', '2A_v_2B', '2C_v_2D'],
                         actual.comparisons)

    def test_phenotypes_comparisons_replicate_tuple(self):
        phenotype_labels =    'pheno1 | pheno2 | pheno3 '
        sample_dict = {'s1' : 'a1     | b1        | c1 ',
                       's2' : 'a1     | b2        | c2 ',
                       's3' : 'a2     | b3        | c2 ',}
        comparison_dict = {'pheno1' : ['1A_v_1B','1C_v_1D'],
                           'pheno2' : ['2A_v_2B','2C_v_2D'],
                           'pheno3' : ['3A_v_3B','3C_v_3D'],}

        delimiter = '|'
        manager = PhenotypeManager({'phenotypes': phenotype_labels,
                                    'samples': sample_dict,
                                    'comparisons' : comparison_dict},
                                    delimiter)

        actual = manager.phenotypes_comparisons_replicates_tuple

        self.assertEqual(['pheno1', 'pheno1', 'pheno3', 'pheno3'],
                         actual.phenotypes)
        self.assertEqual(['1A_v_1B', '1C_v_1D', '3A_v_3B', '3C_v_3D'],
                         actual.comparisons)


    def test_phenotypes_comparisons_all_tuple_oneLabelOneComparison(self):
        comparison_dict = {'phenoLabel1' : ['1A_v_1B']}

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual = manager.phenotypes_comparisons_all_tuple

        self.assertEqual(['phenoLabel1'],
                         actual.phenotypes)
        self.assertEqual(['1A_v_1B'],
                         actual.comparisons)

    def test_phenotypes_comparisons_all_tuple_stripsWhitespace(self):
        comparison_dict = {'phenoLabel1' : [' 1A_v_1B\t',' 1C_v_1D '],
                           'phenoLabel2' : [' 2A_v_2B\t',' 2C_v_2D ']}
        manager = PhenotypeManager({'comparisons' : comparison_dict})

        actual = manager.phenotypes_comparisons_all_tuple

        self.assertEqual(['phenoLabel1', 'phenoLabel1', 'phenoLabel2', 'phenoLabel2'],
                         actual.phenotypes)
        self.assertEqual(['1A_v_1B', '1C_v_1D', '2A_v_2B', '2C_v_2D'],
                         actual.comparisons)


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
