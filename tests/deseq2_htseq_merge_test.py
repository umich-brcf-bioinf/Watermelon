#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division
import os
import unittest
from testfixtures.tempdirectory import TempDirectory


class Deseq2HtseqMergeTest(unittest.TestCase):
    def test_1(self):
        self.assertEqual(1,2)


    # def test_build_sample_metadata_basecase(self):
    #     with TempDirectory() as temp_dir:
    #         temp_dir_path = temp_dir.path
    #         sample_metadata_file_name = os.path.join(temp_dir_path, 'sample_metadata.txt')
    #         deseq2_helper.build_sample_metadata(config, with_replicates, sample_metadata_file_name)
    #         with open(sample_metadata_file_name, 'r') as actual_file:
    #             lines = actual_file.readlines()
    #
    #     line = iter(lines)
    #     self.assertEqual('sample_name\tdiet\tdiet^replicate\tgender\tgender^replicate\tcombinatoric_group\tcombinatoric_group^replicate\n', next(line))
    #     self.assertEqual('sampleA\tHF\t1\tmale\t1\tHF^male\t1\n', next(line))
    #     self.assertEqual('sampleB\tHF\t2\tfemale\t1\tHF^female\t1\n', next(line))
    #     self.assertEqual('sampleC\tLF\t1\tmale\t2\tLF^male\t1\n', next(line))
    #     self.assertEqual('sampleD\tLF\t2\tmale\t3\tLF^male\t2\n', next(line))
    #     self.assertRaises(StopIteration, next, line)
