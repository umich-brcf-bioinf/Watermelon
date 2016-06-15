#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import os
import shutil
import subprocess
import sys
import unittest

from testfixtures.tempdirectory import TempDirectory

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class SnakeFileTest(unittest.TestCase):
    def test_concat_reads(self):
        with TempDirectory() as tmp_dir:
            os.chdir(tmp_dir.path)
            os.mkdir("00-multiplexed_reads")
            os.mkdir("00-multiplexed_reads/Sample_01")
            os.mkdir("00-multiplexed_reads/Sample_02")

            subprocess.check_output("echo 01-A > 00-multiplexed_reads/Sample_01/01_stuff_R1_001.fastq.gz", shell=True)
            subprocess.check_output("echo 01-B > 00-multiplexed_reads/Sample_01/01_stuff_R1_002.fastq.gz", shell=True)
            subprocess.check_output("echo 01-C > 00-multiplexed_reads/Sample_01/01_stuff_R1_003.fastq.gz", shell=True)
            subprocess.check_output("echo 02-A > 00-multiplexed_reads/Sample_02/02_stuff_R1_001.fastq.gz", shell=True)
            subprocess.check_output("echo 02-B > 00-multiplexed_reads/Sample_02/02_stuff_R1_002.fastq.gz", shell=True)
            subprocess.check_output("echo 02-C > 00-multiplexed_reads/Sample_02/02_stuff_R1_003.fastq.gz", shell=True)

            shutil.copy("/ccmb/BioinfCore/SoftwareDev/projects/Watermelon/Snakefile",
                    tmp_dir.path)

            command = "snakemake 2>&1 2>/dev/null"
            output = subprocess.check_output(command, shell=True)
            sample_A_filename = os.path.join(tmp_dir.path,
                                             '01-raw_reads',
                                             'Sample_01_R1.fastq.gz')
            with open(sample_A_filename, 'r') as sample_A:
                sample_A_contents = sample_A.read()
        self.assertEquals('01-A\n01-B\n01-C\n', sample_A_contents)