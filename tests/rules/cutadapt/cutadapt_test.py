from __future__ import print_function, absolute_import
import filecmp
import gzip
import os
import shutil
import subprocess
import sys
import unittest

from testfixtures import TempDirectory

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
SNAKEFILE_PATH = os.path.join(TEST_DIR, '..', '..', '..', 'rnaseq.snakefile')
CONFIGFILE_PATH = os.path.join(TEST_DIR, 'end_trim_test', 'end_trim.yaml')

def gunzip(source_file, dest_file):
    with gzip.GzipFile(source_file, 'rb') as inF, \
         open(dest_file, 'wb') as outF:
        data = inF.read()
        outF.write(data)

class CutadaptTest(unittest.TestCase):
    def test_end_trim(self):
        source_dir = os.path.join(TEST_DIR, 'end_trim_test', 'output')
        expected_dir = os.path.join(TEST_DIR, 'end_trim_test', 'expected', '02-cutadapt')
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            dest_dir = os.path.join(temp_dir_path, 'output')
            shutil.copytree(source_dir, dest_dir) 
            os.chdir(dest_dir)
            command = '''snakemake --cores 2 \
     --snakefile {0} \
     --configfile {1} \
     --force 02-cutadapt/Sample_0_trimmed_R1.fastq.gz 02-cutadapt/Sample_1_trimmed_R1.fastq.gz
'''.format(SNAKEFILE_PATH, CONFIGFILE_PATH)
            subprocess.check_output(command, shell=True)
            gunzip('02-cutadapt/Sample_0_trimmed_R1.fastq.gz', '02-cutadapt/Sample_0_trimmed_R1.fastq')
            gunzip('02-cutadapt/Sample_1_trimmed_R1.fastq.gz', '02-cutadapt/Sample_1_trimmed_R1.fastq')
            s0_matched = filecmp.cmp(expected_dir + '/Sample_0_trimmed_R1.fastq',
                        '02-cutadapt/Sample_0_trimmed_R1.fastq',
                        shallow=True)
            s1_matched = filecmp.cmp(expected_dir + '/Sample_1_trimmed_R1.fastq',
                        '02-cutadapt/Sample_1_trimmed_R1.fastq',
                        shallow=True)

        self.assertTrue(s0_matched, "s0")
        self.assertTrue(s1_matched, "s1")



# HOME_DIR=Watermelon/
# TEST_DIR=Watermelon/tests/rules/cutadapt/end_trim_test/
# TEMP_DIR=/tmp/output/
# cd $TEST_DIR
# cp -a $TEST_DIR/output/ $TEMP_DIR
# cd $TEMP_DIR
# # python can do everything above
# snakemake --cores 2 \
#     --snakefile $HOME_DIR/Snakefile \
#     --configfile $TEST_DIR/config.yaml \
#     --force 02-cutadapt/Sample_0_trimmed_R1.fastq.gz 02-cutadapt/Sample_1_trimmed_R1.fastq.gz
# # python can do everything below
# gunzip 02-cutadapt/*.gz
# diff -r $TEST_DIR/expected/02-cutadapt 02-cutadapt