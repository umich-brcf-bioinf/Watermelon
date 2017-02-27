#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division
import os
import unittest
from testfixtures.tempdirectory import TempDirectory
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

    def test_build_sample_metadata_basecase(self):
        config = {'phenotypes': ' diet ^ gender ',
                  'samples' : {'sampleA' : 'HF ^ male'}}

        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_metadata_file_name = os.path.join(temp_dir_path, 'sample_metadata.txt')
            deseq2.build_sample_metadata(config, sample_metadata_file_name)
            with open(sample_metadata_file_name, 'r') as actual_file:
                lines = actual_file.readlines()

        self.assertEqual(2, len(lines))
        line = iter(lines)
        self.assertEqual('phenotype\tdiet\tgender\n', next(line))
        self.assertEqual('sampleA\tHF\tmale\n', next(line))