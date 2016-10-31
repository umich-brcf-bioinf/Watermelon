#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import os
import time
import unittest

from testfixtures.tempdirectory import TempDirectory

import scripts.rnaseq_snakefile_helper as rnaseq_snakefile_helper

class ConfigurationTest(unittest.TestCase):
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

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'A.watermelon.md5', '.watermelon.md5:acbd18db4cc2f85cedef654fccc4a4d8')
            self.assertChecksumFile(config_dir, 'B.watermelon.md5', '.watermelon.md5:37b51d194a7513e45b56f6524f2d51f2')
            self.assertChecksumFile(config_dir, 'C.watermelon.md5', '.watermelon.md5:73feffa4b7f6bb68e44cf984c85f6e88')

    def test_checksum_reset_all_doesNotOverwriteUnchangedKeys(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : 'foo',
                      'B' : 'bar',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            def mod_time(filename):
                return os.path.getmtime(os.path.join(config_dir, filename))

            checksum_A1_timestamp = mod_time('A.watermelon.md5')
            checksum_B1_timestamp = mod_time('B.watermelon.md5')
            checksum_C1_timestamp = mod_time('C.watermelon.md5')

            time.sleep(2)

            config = {'A' : 'foo2',
                      'B' : 'bar',
                      'C' : 'baz2'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            checksum_A2_timestamp = mod_time('A.watermelon.md5')
            checksum_B2_timestamp = mod_time('B.watermelon.md5')
            checksum_C2_timestamp = mod_time('C.watermelon.md5')

    def test_checksum_reset_all_addsNewFilesForNewKeys(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : 'foo',
                      'B' : 'bar'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            config = {'A' : 'foo',
                      'B' : 'bar',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'A.watermelon.md5', '.watermelon.md5:acbd18db4cc2f85cedef654fccc4a4d8')
            self.assertChecksumFile(config_dir, 'B.watermelon.md5', '.watermelon.md5:37b51d194a7513e45b56f6524f2d51f2')
            self.assertChecksumFile(config_dir, 'C.watermelon.md5', '.watermelon.md5:73feffa4b7f6bb68e44cf984c85f6e88')

    def test_checksum_reset_all_removesChecksumFilesForMissingKeys(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = temp_dir_path
            config = {'A' : 'foo',
                      'B' : 'bar',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            config = {'A' : 'foo',
                      'C' : 'baz'}

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            actual_files = sorted([f for f in os.listdir(config_dir) if os.path.isfile(os.path.join(config_dir, f))])
            self.assertEqual(['A.watermelon.md5', 'C.watermelon.md5'], actual_files)

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

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'samples.watermelon.md5', '.watermelon.md5:94a5e110e7244bc21d27fc4e563346f4')
            self.assertChecksumFile(config_dir, 'comparisons.watermelon.md5', '.watermelon.md5:540157f53ef390438a880625b6d5beb4')

            config = {'samples' : {'s1':'1',
                                   's2':'2',
                                   's3':'3'},
                      'comparisons' : {'a_b':'a vs b',
                                       'c_d':'c vs d'}
                      }
            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'samples.watermelon.md5', '.watermelon.md5:94a5e110e7244bc21d27fc4e563346f4')
            self.assertChecksumFile(config_dir, 'comparisons.watermelon.md5', '.watermelon.md5:540157f53ef390438a880625b6d5beb4')

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

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'samples.watermelon.md5', '.watermelon.md5:94a5e110e7244bc21d27fc4e563346f4')
            self.assertChecksumFile(config_dir, 'comparisons.watermelon.md5', '.watermelon.md5:540157f53ef390438a880625b6d5beb4')

            config = {'samples' : {'s1':'1',
                                   's2':'2'},
                      'comparisons' : {'a_b':'a vs b'}
                      }
            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'samples.watermelon.md5', '.watermelon.md5:1c9f24c825613d5fd1cdc2df726d5bb8')
            self.assertChecksumFile(config_dir, 'comparisons.watermelon.md5', '.watermelon.md5:18a615d74c5db4f9edc688955a8c7516')

    def test_checksum_reset_all_createsChecksumDir(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_dir = os.path.join(temp_dir_path, 'parent', 'child', 'config_checksums')
            config = {'samples' : {'s1':'1',
                                   's2':'2'},
                      'comparisons' : {'a_b':'a vs b',
                                       'c_d':'c vs d'}
                      }

            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'samples.watermelon.md5')
            self.assertChecksumFile(config_dir, 'comparisons.watermelon.md5')

            config = {'samples' : {'s1':'1',
                                   's2':'2'},
                      'comparisons' : {'a_b':'a vs b'},
                      'references' : 'A'}
            rnaseq_snakefile_helper.checksum_reset_all(config_dir, config)

            self.assertChecksumFile(config_dir, 'samples.watermelon.md5')
            self.assertChecksumFile(config_dir, 'comparisons.watermelon.md5')
            self.assertChecksumFile(config_dir, 'references.watermelon.md5')
