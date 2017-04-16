#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division
import os
import unittest
from testfixtures.tempdirectory import TempDirectory
import scripts.deseq2_helper as deseq2_helper


class Deseq2Test(unittest.TestCase):
    def test_build_sample_metadata_list_basecase(self):
        config = {'main_factors': 'yes  ^ yes    ^ no',
                  'phenotypes'  : 'diet ^ gender ^ male.diet',
                  'samples'     : 
                     {'sampleA' : 'HF   ^ male   ^ HF',
                      'sampleB' : 'LF   ^ female ^',
                      'sampleC' : '     ^ male   ^ ',
                      'sampleD' : 'HF   ^        ^ ',
                      'sampleE' : '     ^        ^ ',
                      'sampleF' : '     ^        ^ '}}

        lines = deseq2_helper._build_sample_metadata_list(config)

        line = iter(lines)
        self.assertEqual(['sample_name', 'diet', 'diet^replicate', 'gender', 'gender^replicate', 'male.diet', 'male.diet^replicate', 'combinatoric_group', 'combinatoric_group^replicate'], next(line))
        self.assertEqual(['sampleA',   'HF', '1', 'male',   '1',  'HF', '1', 'HF^male',  '1'], next(line))
        self.assertEqual(['sampleB',   'LF', '1', 'female', '1',  '',    '', 'LF^female', '1'], next(line))
        self.assertEqual(['sampleC',   '',   '',  'male',   '2',  '',    '', '^male',     '1'], next(line))
        self.assertEqual(['sampleD',   'HF', '2', '', '',         '',    '', 'HF^',       '1'], next(line))
        self.assertEqual(['sampleE',   '', '',    '', '',         '',    '', '^' ,        '1'], next(line))
        self.assertEqual(['sampleF',   '', '',    '', '',         '',    '', '^' ,        '2'], next(line))
        self.assertRaises(StopIteration, next, line)

    def test_build_sample_metadata_basecase(self):
        config = {'main_factors': 'yes   ^ yes',
                  'phenotypes'  : 'diet  ^ gender ',
                  'samples'     :
                      {'sampleA' : 'HF   ^ male'}}

        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            sample_metadata_file_name = os.path.join(temp_dir_path, 'sample_metadata.txt')
            deseq2_helper.build_sample_metadata(config, sample_metadata_file_name)
            with open(sample_metadata_file_name, 'r') as actual_file:
                lines = actual_file.readlines()

        line = iter(lines)
        self.assertEqual('sample_name\tdiet\tdiet^replicate\tgender\tgender^replicate\tcombinatoric_group\tcombinatoric_group^replicate\n', next(line))
        self.assertEqual('sampleA\tHF\t1\tmale\t1\tHF^male\t1\n', next(line))
        self.assertRaises(StopIteration, next, line)

    def test_build_contrasts_list_basecase(self):
        config = {'comparisons': 
                                {'gender'    : ['male_v_female'],
                                'diet'       : ['HF_v_ND', 'HF_v_DRG'],
                                'male.diet'  : ['HF_v_ND', 'DRG_v_ND'],
                                'female.diet': ['HF_v_ND']}}
        

        lines = deseq2_helper._build_contrasts_list(config)

        line = iter(lines)
        self.assertEqual(['factor', 'test_level', 'reference_level', 'base_file_name'], next(line))
        self.assertEqual(['diet', 'HF', 'DRG',  'HF_v_DRG'], next(line))
        self.assertEqual(['diet', 'HF', 'ND',  'HF_v_ND'], next(line))
        self.assertEqual(['female.diet', 'HF', 'ND', 'HF_v_ND'], next(line))
        self.assertEqual(['gender', 'male', 'female', 'male_v_female'], next(line))
        self.assertEqual(['male.diet', 'DRG', 'ND', 'DRG_v_ND'], next(line))
        self.assertEqual(['male.diet', 'HF', 'ND', 'HF_v_ND'], next(line))
        self.assertRaises(StopIteration, next, line)

    def test_build_contrasts_basecase(self):
        config = {'comparisons': {'gender' : ['male_v_female']}}
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            contrast_file_name = os.path.join(temp_dir_path, 'deseq2_contrast.txt')

            deseq2_helper.build_contrasts(config, contrast_file_name)
            with open(contrast_file_name, 'r') as actual_file:
                lines = actual_file.readlines()

        line = iter(lines)
        self.assertEqual('factor\ttest_level\treference_level\tbase_file_name\n', next(line))
        self.assertEqual('gender\tmale\tfemale\tmale_v_female\n', next(line))
        self.assertRaises(StopIteration, next, line)
