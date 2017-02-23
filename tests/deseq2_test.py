#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import unittest

import scripts.deseq2 as deseq2

class Deseq2Test(unittest.TestCase):
    def test_build_sample_metadata_list_basecase(self):
        config = {'phenotypes': ' diet ^ gender ',
                  'samples' : {'sampleA' : 'HF ^ male',
                               'sampleB' : 'LF ^ female',
                               'sampleC' : '   ^ male',
                               'sampleD' : 'HF ^     ',
                               'sampleE' : '   ^     '}}

        lines = deseq2._build_sample_metadata_list(config)

        self.assertEqual(6, len(lines))
        line = iter(lines)
        self.assertEqual(['phenotype', 'diet', 'gender'], next(line))
        self.assertEqual(['sampleA', 'HF', 'male'], next(line))
        self.assertEqual(['sampleB','LF','female'], next(line))
        self.assertEqual(['sampleC','','male'], next(line))
        self.assertEqual(['sampleD','HF',''], next(line))
        self.assertEqual(['sampleE','',''], next(line))