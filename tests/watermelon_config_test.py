#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import unittest

from scripts import watermelon_config

class WatermelonConfigTest(unittest.TestCase):

    def test_split_config_boolean_true(self):
        self.assertEqual(['yes'], watermelon_config.split_config_list(True))

    def test_split_config_boolean_false(self):
        self.assertEqual(['no'], watermelon_config.split_config_list(False))

    def test_split_config_list(self):
        self.assertEqual(['A','B','C'], watermelon_config.split_config_list('A^B^C'))

    def test_split_config_list_customDelim(self):
        self.assertEqual(['A','B','C'],
                         watermelon_config.split_config_list('A|B|C', '|'))

    def test_split_config_list_stripsWhitespace(self):
        self.assertEqual(['A','B1 B2','C'],
                         watermelon_config.split_config_list('  A  | B1 B2 |\t\tC\n', '|'))

    def test_transform_config_passthrough(self):
        config = {'a': 1, 'b': 2}
        expected_config = dict(config)

        watermelon_config.transform_config(config)
        self.assertEqual(expected_config, config)

    def test_transform_config_obsolete_dir_keys(self):
        config = {'a': 1,
                  'b': 2,
                  'input_dir': 'Tuttle_RS1_tests/00-multiplexed_reads',
                  'alignment_output_dir': 'alignment1',
                  'diffex_output_dir': 'diffex1',
                  'deliverables_output_dir': 'deliverables1',
                  }
        expected_config = {\
                  'a': 1,
                  'b': 2,
                  'dirs': {'input': 'Tuttle_RS1_tests/00-multiplexed_reads',
                           'alignment_output': 'alignment1',
                           'diffex_output': 'diffex1',
                           'deliverables_output': 'deliverables1',
                           }
                  }
        watermelon_config.transform_config(config)
        self.assertEqual(expected_config, config)

    def test_transform_config_mergesStrayKeysIntoExistingDirKey(self):
        config = {'a': 1,
                  'b': 2,
                  'dirs': {'input': 'Tuttle_RS1_tests/00-multiplexed_reads',
                           'alignment_output': 'alignment1',
                           },
                  'diffex_output_dir': 'diffex1',
                  'deliverables_output_dir': 'deliverables1',
                  }
        expected_config = {\
                  'a': 1,
                  'b': 2,
                  'dirs': {'input': 'Tuttle_RS1_tests/00-multiplexed_reads',
                           'alignment_output': 'alignment1',
                           'diffex_output': 'diffex1',
                           'deliverables_output': 'deliverables1',
                           }
                  }
        watermelon_config.transform_config(config)
        self.assertEqual(expected_config, config)

    def test_transform_config_raisesOnDuplicateKeys(self):
        config = {'a': 1,
                  'b': 2,
                  'dirs': {'input': 'new_input',
                           'alignment_output': 'new_alignment',
                           },
                  'input_dir': 'old_input',
                  'alignment_output_dir': 'old_alignment',
                  }
        self.assertRaisesRegexp(ValueError,
                                (r'found both old and new style.*'
                                 r'alignment_output_dir,dirs:alignment_output; '
                                 r'input_dir,dirs:input'),
                                watermelon_config.transform_config,
                                config)
