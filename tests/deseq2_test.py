#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division
import os
import unittest
from testfixtures.tempdirectory import TempDirectory
import scripts.deseq2 as deseq2


class Deseq2Test(unittest.TestCase):
    def test_build_sample_metadata_list_basecase(self):
        config = {'main_factors': 'yes  ^ yes    ^ no',
                  'phenotypes'  : 'diet ^ gender ^ male.diet',
                  'samples'     : 
                     {'sampleA' : 'HF   ^ male   ^ HF',
                      'sampleB' : 'LF   ^ female ^',
                      'sampleC' : '     ^ male   ^ ',
                      'sampleD' : 'HF   ^        ^ ',
                      'sampleE' : '     ^        ^ '}}

        lines = deseq2._build_sample_metadata_list(config)

        self.assertEqual(6, len(lines))
        line = iter(lines)
        self.assertEqual(['phenotype', 'diet', 'gender', 'male.diet', 'combinatoric_group'], next(line))
        self.assertEqual(['sampleA', 'HF', 'male', 'HF', 'HF^male'], next(line))
        self.assertEqual(['sampleB','LF','female', '', 'LF^female'], next(line))
        self.assertEqual(['sampleC','','male', '', '^male'], next(line))
        self.assertEqual(['sampleD','HF','', '', 'HF^'], next(line))
        self.assertEqual(['sampleE','','', '', '^'], next(line))

    def test_build_sample_metadata_basecase(self):
        config = {'main_factors': 'yes   ^ yes',
                  'phenotypes'  : 'diet  ^ gender ',
                  'samples'     :
                      {'sampleA' : 'HF   ^ male'}}

        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_metadata_file_name = os.path.join(temp_dir_path, 'sample_metadata.txt')
            deseq2.build_sample_metadata(config, sample_metadata_file_name)
            with open(sample_metadata_file_name, 'r') as actual_file:
                lines = actual_file.readlines()

        self.assertEqual(2, len(lines))
        line = iter(lines)
        self.assertEqual('phenotype\tdiet\tgender\tcombinatoric_group\n', next(line))
        self.assertEqual('sampleA\tHF\tmale\tHF^male\n', next(line))