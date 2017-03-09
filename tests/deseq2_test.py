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

        line = iter(lines)
        self.assertEqual(['phenotype', 'diet', 'gender', 'male.diet', 'combinatoric_group'], next(line))
        self.assertEqual(['sampleA', 'HF', 'male', 'HF', 'HF^male'], next(line))
        self.assertEqual(['sampleB','LF','female', '', 'LF^female'], next(line))
        self.assertEqual(['sampleC','','male', '', '^male'], next(line))
        self.assertEqual(['sampleD','HF','', '', 'HF^'], next(line))
        self.assertEqual(['sampleE','','', '', '^'], next(line))
        self.assertRaises(StopIteration, next, line)

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

        line = iter(lines)
        self.assertEqual('phenotype\tdiet\tgender\tcombinatoric_group\n', next(line))
        self.assertEqual('sampleA\tHF\tmale\tHF^male\n', next(line))
        self.assertRaises(StopIteration, next, line)

    def test_build_contrasts_list_basecase(self):
        config = {'comparisons': 
                                {'gender'    : ['male_v_female'],
                                'diet'       : ['HF_v_ND', 'HF_v_DRG'],
                                'male.diet'  : ['HF_v_ND', 'DRG_v_ND'],
                                'female.diet': ['HF_v_ND']}}
        

        lines = deseq2._build_contrasts_list(config, '/tmp')

        line = iter(lines)
        self.assertEqual(['factor', 'test_level', 'reference_level', 'directory_name', 'base_file_name'], next(line))
        self.assertEqual(['diet', 'HF', 'DRG', '/tmp/diet', 'HF_v_DRG'], next(line))
        self.assertEqual(['diet', 'HF', 'ND', '/tmp/diet', 'HF_v_ND'], next(line))
        self.assertEqual(['female.diet', 'HF', 'ND', '/tmp/female.diet', 'HF_v_ND'], next(line))
        self.assertEqual(['gender', 'male', 'female', '/tmp/gender', 'male_v_female'], next(line))
        self.assertEqual(['male.diet', 'DRG', 'ND', '/tmp/male.diet', 'DRG_v_ND'], next(line))
        self.assertEqual(['male.diet', 'HF', 'ND', '/tmp/male.diet', 'HF_v_ND'], next(line))
        self.assertRaises(StopIteration, next, line)

    def test_build_contrasts_basecase(self):
        config = {'comparisons': {'gender' : ['male_v_female']}}
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            contrast_file_name = os.path.join(temp_dir_path, 'deseq2_contrast.txt')
            comparison_file_prefix = 'analysis_02_27/diffex_results/16-deseq2'

            deseq2.build_contrasts(config, comparison_file_prefix, contrast_file_name)
            with open(contrast_file_name, 'r') as actual_file:
                lines = actual_file.readlines()

        line = iter(lines)
        self.assertEqual('factor\ttest_level\treference_level\tdirectory_name\tbase_file_name\n', next(line))
        self.assertEqual('gender\tmale\tfemale\tanalysis_02_27/diffex_results/16-deseq2/gender\tmale_v_female\n', next(line))
        self.assertRaises(StopIteration, next, line)
