#pylint: disable=too-many-public-methods, invalid-name, no-self-use
import subprocess
import sys
import unittest

class VersionTest(unittest.TestCase):
    def test_snakemake_version(self):
        output = str(subprocess.check_output("snakemake --version",
                    shell=True))
        self.assertRegexpMatches(output, r"5\..*", "wrong snakemake version")

    def test_mutt_version(self):
        output = str(subprocess.check_output("mutt -v | head -1",
                    shell=True))
        self.assertRegexpMatches(output, r"Mutt", "mutt not installed")
