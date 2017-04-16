#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import unittest

from scripts import watermelon_config

class WatermelonConfigTest(unittest.TestCase):

    def test_split_config_list(self):
        self.assertEqual(['A','B','C'], watermelon_config.split_config_list('A^B^C'))

    def test_split_config_list_customDelim(self):
        self.assertEqual(['A','B','C'],
                         watermelon_config.split_config_list('A|B|C', '|'))

    def test_split_config_list_stripsWhitespace(self):
        self.assertEqual(['A','B1 B2','C'],
                         watermelon_config.split_config_list('  A  | B1 B2 |\t\tC\n', '|'))
