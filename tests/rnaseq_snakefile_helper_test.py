#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import os
import time
import unittest

from testfixtures.tempdirectory import TempDirectory

import scripts.rnaseq_snakefile_helper as rnaseq_snakefile_helper
from scripts.rnaseq_snakefile_helper import PhenotypeManager

class ChecksumManagerTest(unittest.TestCase):
    def assertChecksumFile(self, config_dir, config_filename, expected_checksum=None):
        config_path = os.path.join(config_dir, config_filename)
        self.assertTrue(os.path.exists(config_path))
        with open(config_path, 'r') as file:
            lines = file.readlines()
        self.assertEqual(1, len(lines))
        if expected_checksum:
            self.assertEqual(expected_checksum, lines[0])

    def test_checksum_reset_all_initializesConfig(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : 'foo',
                      'B' : 'bar',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, configThing=config)

            self.assertChecksumFile(config_dir, 'configThing-A.watermelon.md5', '.watermelon.md5:acbd18db4cc2f85cedef654fccc4a4d8')
            self.assertChecksumFile(config_dir, 'configThing-B.watermelon.md5', '.watermelon.md5:37b51d194a7513e45b56f6524f2d51f2')
            self.assertChecksumFile(config_dir, 'configThing-C.watermelon.md5', '.watermelon.md5:73feffa4b7f6bb68e44cf984c85f6e88')

    def test_checksum_reset_all_doesNotOverwriteUnchangedKeysForSimpleValues(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : 'foo',
                      'B' : 'bar',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            def mod_time(filename):
                return os.path.getmtime(os.path.join(config_dir, filename))

            checksum_A1_timestamp = mod_time('config-A.watermelon.md5')
            checksum_B1_timestamp = mod_time('config-B.watermelon.md5')
            checksum_C1_timestamp = mod_time('config-C.watermelon.md5')

            time.sleep(1)

            config = {'A' : 'foo2',
                      'B' : 'bar',
                      'C' : 'baz2'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            checksum_A2_timestamp = mod_time('config-A.watermelon.md5')
            checksum_B2_timestamp = mod_time('config-B.watermelon.md5')
            checksum_C2_timestamp = mod_time('config-C.watermelon.md5')

    def test_checksum_reset_all_doesNotOverwriteUnchangedKeysForComplexValues(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : {'1': 'foo', '2': 'bar'},
                      'B' : {'1': 'foo', '2': 'bar'},
                      'C' : {'1': 'foo', '2': 'bar'}}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            def mod_time(filename):
                return os.path.getmtime(os.path.join(config_dir, filename))

            checksum_A1_timestamp = mod_time('config-A.watermelon.md5')
            checksum_B1_timestamp = mod_time('config-B.watermelon.md5')
            checksum_C1_timestamp = mod_time('config-C.watermelon.md5')

            time.sleep(1)

            config = {'A' : {'2': 'bar', '1': 'foo'},
                      'B' : {'2': 'bar', '1': 'foo'},
                      'C' : {'2': 'bar', '1': 'foo'}}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            checksum_A2_timestamp = mod_time('config-A.watermelon.md5')
            checksum_B2_timestamp = mod_time('config-B.watermelon.md5')
            checksum_C2_timestamp = mod_time('config-C.watermelon.md5')


    def test_checksum_reset_all_addsNewFilesForNewKeys(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : 'foo',
                      'B' : 'bar'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            config = {'A' : 'foo',
                      'B' : 'bar',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            self.assertChecksumFile(config_dir, 'config-A.watermelon.md5', '.watermelon.md5:acbd18db4cc2f85cedef654fccc4a4d8')
            self.assertChecksumFile(config_dir, 'config-B.watermelon.md5', '.watermelon.md5:37b51d194a7513e45b56f6524f2d51f2')
            self.assertChecksumFile(config_dir, 'config-C.watermelon.md5', '.watermelon.md5:73feffa4b7f6bb68e44cf984c85f6e88')

    def test_checksum_reset_all_removesChecksumFilesForMissingKeys(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : 'foo',
                      'B' : 'bar',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            config = {'A' : 'foo',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            actual_files = sorted([f for f in os.listdir(config_dir) if os.path.isfile(os.path.join(config_dir, f))])
            self.assertEqual(['config-A.watermelon.md5', 'config-C.watermelon.md5'],
                             actual_files)

    def test_checksum_reset_all_equivalentObjectsDoNotResetChecksum(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'samples' : {'s1':'1',
                                   's2':'2',
                                   's3':'3'},
                      'comparisons' : {'a_b':'a vs b',
                                       'c_d':'c vs d'}
                      }

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            self.assertChecksumFile(config_dir, 'config-samples.watermelon.md5', '.watermelon.md5:94a5e110e7244bc21d27fc4e563346f4')
            self.assertChecksumFile(config_dir, 'config-comparisons.watermelon.md5', '.watermelon.md5:540157f53ef390438a880625b6d5beb4')

            config = {'samples' : {'s1':'1',
                                   's2':'2',
                                   's3':'3'},
                      'comparisons' : {'a_b':'a vs b',
                                       'c_d':'c vs d'}
                      }
            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            self.assertChecksumFile(config_dir, 'config-samples.watermelon.md5', '.watermelon.md5:94a5e110e7244bc21d27fc4e563346f4')
            self.assertChecksumFile(config_dir, 'config-comparisons.watermelon.md5', '.watermelon.md5:540157f53ef390438a880625b6d5beb4')

    def test_checksum_reset_all_distinctObjectsResetsChecksum(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'samples' : {'s1':'1',
                                   's2':'2',
                                   's3':'3'},
                      'comparisons' : {'a_b':'a vs b',
                                       'c_d':'c vs d'}
                      }

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            self.assertChecksumFile(config_dir, 'config-samples.watermelon.md5', '.watermelon.md5:94a5e110e7244bc21d27fc4e563346f4')
            self.assertChecksumFile(config_dir, 'config-comparisons.watermelon.md5', '.watermelon.md5:540157f53ef390438a880625b6d5beb4')

            config = {'samples' : {'s1':'1',
                                   's2':'2'},
                      'comparisons' : {'a_b':'a vs b'}
                      }
            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            self.assertChecksumFile(config_dir, 'config-samples.watermelon.md5', '.watermelon.md5:1c9f24c825613d5fd1cdc2df726d5bb8')
            self.assertChecksumFile(config_dir, 'config-comparisons.watermelon.md5', '.watermelon.md5:18a615d74c5db4f9edc688955a8c7516')

    def test_checksum_reset_all_createsChecksumDir(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = os.path.join(temp_dir_path, 'parent', 'child', 'config_checksums')
            config = {'samples' : {'s1':'1',
                                   's2':'2'},
                      'comparisons' : {'a_b':'a vs b',
                                       'c_d':'c vs d'}
                      }

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            self.assertChecksumFile(config_dir, 'config-samples.watermelon.md5')
            self.assertChecksumFile(config_dir, 'config-comparisons.watermelon.md5')

            config = {'samples' : {'s1':'1',
                                   's2':'2'},
                      'comparisons' : {'a_b':'a vs b'},
                      'references' : 'A'}
            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config=config)

            self.assertChecksumFile(config_dir, 'config-samples.watermelon.md5')
            self.assertChecksumFile(config_dir, 'config-comparisons.watermelon.md5')
            self.assertChecksumFile(config_dir, 'config-references.watermelon.md5')

    def test_checksum_reset_all_multiple_dicts(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            checksum_dir = temp_dir_path
            config = {'samples' : {'s1':'1',
                                   's2':'2',
                                   's3':'3'}}
            phenotype_samples = {'gender': {'M':['s1','s2'], 'F': ['s3']},
                                 'diet' : {'H' :['s1','s2'], 'L': ['s3']}}
            phenotype_comparisons = {'gender': 'M_v_F',
                                     'diet' : 'L_v_H'}

            rnaseq_snakefile_helper.checksum_reset_all(checksum_dir,
                                                       config=config,
                                                       phenotype_samples=phenotype_samples,
                                                       phenotype_comparisons=phenotype_comparisons)

            self.assertChecksumFile(checksum_dir,
                                    'config-samples.watermelon.md5',
                                    '.watermelon.md5:94a5e110e7244bc21d27fc4e563346f4')
            self.assertChecksumFile(checksum_dir,
                                    'phenotype_samples-gender.watermelon.md5',
                                    '.watermelon.md5:2ceb07f4677c8a6255e4913ce26c0c95')
            self.assertChecksumFile(checksum_dir,
                                    'phenotype_samples-diet.watermelon.md5',
                                    '.watermelon.md5:6e615144105118241c3bf19e47d36580')
            self.assertChecksumFile(checksum_dir,
                                    'phenotype_comparisons-gender.watermelon.md5',
                                    '.watermelon.md5:8253a6220bedb78b4b4380b152d518e0')
            self.assertChecksumFile(checksum_dir,
                                    'phenotype_comparisons-diet.watermelon.md5',
                                    '.watermelon.md5:0ef89ab1c75b2de0359c4dff6363105f')

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

    def test_separated_comparisons(self):
        comparison_dict = {'phenoLabel1' : '1A_v_1B\n1C_v_1D',
                           'phenoLabel2' : '2A_v_2B\n2C_v_2D'}
        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.separated_comparisons('#')

        expected_comparisons = {'phenoLabel1' : '1A_v_1B#1C_v_1D',
                                'phenoLabel2' : '2A_v_2B#2C_v_2D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_separated_comparisons_stripsWhitespace(self):
        comparison_dict = {'phenoLabel1' : '\t1A_v_1B\t\n\t1C_v_1D\t',
                           'phenoLabel2' : '  2A_v_2B  \n  2C_v_2D  '}
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
        comparison_dict = {'phenoLabel1' : '1A_v_1B\n1C_v_1D',
                           'phenoLabel2' : '2A_v_2B\n2C_v_2D'}
        comparison_infix = '_v_'
        delimiter = '|^|'
        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)
        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D',
                                'phenoLabel2' : '2A|^|2B|^|2C|^|2D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_concatenated_comparison_values_stripsWhitespace(self):
        comparison_dict = {'phenoLabel1' : ' \t1A_v_1B \t\n\t 1C_v_1D \t'}
        comparison_infix = '_v_'
        delimiter = '|^|'

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)

        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_concatenated_comparison_values_duplicateValuesCoalesced(self):
        comparison_dict = {'phenoLabel1' : '1A_v_1B\n1A_v_1C\n1A_v_1D'}
        comparison_infix = '_v_'
        delimiter = '|^|'

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)

        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D'}
        self.assertEqual(expected_comparisons, actual_comparisons)


    def test_concatenated_comparison_values_reorders(self):
        comparison_dict = {'phenoLabel1' : '1D_v_1B\n1C_v_1A'}
        comparison_infix = '_v_'
        delimiter = '|^|'

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual_comparisons = manager.concatenated_comparison_values(delimiter)

        expected_comparisons = {'phenoLabel1' : '1A|^|1B|^|1C|^|1D'}
        self.assertEqual(expected_comparisons, actual_comparisons)

    def test_phenotypes_comparisons_tuple(self):
        comparison_dict = {'phenoLabel1' : '1A_v_1B\n1C_v_1D',
                           'phenoLabel2' : '2A_v_2B\n2C_v_2D'}

        manager = PhenotypeManager({'comparisons' : comparison_dict})

        actual = manager.phenotypes_comparisons_tuple

        self.assertEqual(['phenoLabel1', 'phenoLabel1', 'phenoLabel2', 'phenoLabel2'],
                         actual.phenotypes)
        self.assertEqual(['1A_v_1B', '1C_v_1D', '2A_v_2B', '2C_v_2D'],
                         actual.comparisons)

    def test_phenotypes_comparisons_tuple_oneLabelOneComparison(self):
        comparison_dict = {'phenoLabel1' : '1A_v_1B'}

        manager = PhenotypeManager({'comparisons' : comparison_dict})
        actual = manager.phenotypes_comparisons_tuple

        self.assertEqual(['phenoLabel1'],
                         actual.phenotypes)
        self.assertEqual(['1A_v_1B'],
                         actual.comparisons)

    def test_phenotypes_comparisons_tuple_stripsWhitespace(self):
        comparison_dict = {'phenoLabel1' : ' 1A_v_1B\t\n 1C_v_1D ',
                           'phenoLabel2' : ' 2A_v_2B\t\n 2C_v_2D '}
        manager = PhenotypeManager({'comparisons' : comparison_dict})

        actual = manager.phenotypes_comparisons_tuple

        self.assertEqual(['phenoLabel1', 'phenoLabel1', 'phenoLabel2', 'phenoLabel2'],
                         actual.phenotypes)
        self.assertEqual(['1A_v_1B', '1C_v_1D', '2A_v_2B', '2C_v_2D'],
                         actual.comparisons)

