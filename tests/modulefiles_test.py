#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import os
import shlex
import subprocess
import sys
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

_TEST_DIR = os.path.realpath(os.path.dirname(__file__))
_MODULES_DIR = os.path.realpath(os.path.join(_TEST_DIR, '..', 'modulefiles'))

class BfxCoreBaseTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.dev_null = open("/tmp/rnaseq_stdout.tmp","w")

    def tearDown(self):
        self.dev_null.close()
        unittest.TestCase.tearDown(self)

    @staticmethod
    def build_command(command_text):
        command = ("module purge; "
                   "module use {}; "
                   "module load watermelon_dependencies 2>/dev/null; "
                   "{}").format(_MODULES_DIR, command_text)
        return command

    def check_command(self, command, expected_regexp, message):
        #print(command)
        output = str(subprocess.check_output(command, shell=True))
        self.assertRegexpMatches(output, expected_regexp, message)

class WatermelonModuleTest(BfxCoreBaseTestCase):
    def test_module_switches_python_version(self):
        command = ("module purge; "
                   "module use {}; "
                   "module load watermelon 2>/dev/null; "
                   "python --version 2>&1").format(_MODULES_DIR)
        self.check_command(command, "Python 3.6.1", "wrong python version")

    def test_python3_modules_present(self):
        missing_modules = []
        modules = ["numpy",
                   "pandas",
                   "testfixtures",
                   "xlsxwriter",
                   "yaml"]
        for module_name in modules:
            command = ("module purge; "
                       "module load python/3.6.1; "
                       "python -c 'import {}'".format(module_name))
            try:
                subprocess.check_output(command, stderr=self.dev_null, shell=True)
            except subprocess.CalledProcessError:
                missing_modules.append(module_name)
        self.assertEquals([], missing_modules)

    def test_snakemake_version(self):
        command = ("module purge; "
                   "module use {}; "
                   "module load watermelon 2>/dev/null; "
                   "python -c 'import snakemake as s; print(s.__version__)'").format(_MODULES_DIR)
        output = str(subprocess.check_output(command, stderr=self.dev_null, shell=True))
        self.assertRegexpMatches(output, r"3\..*", "wrong snakemake version")


class WatermelonRnaseqModuleTest(BfxCoreBaseTestCase):
    def test_module_switches_python_version(self):
        command = ("module purge; "
                   "module use {}; "
                   "module load python/3.6.1; "
                   "module load watermelon_dependencies 2>/dev/null; "
                   "python --version 2>&1").format(_MODULES_DIR)
        self.check_command(command, "Python 2.7.9", "wrong python version")

    def test_bbmap_versions(self):
        command = self.build_command("bbmap.sh --version 2>&1 | grep 'BBMap version' || echo 'not installed'")
        self.check_command(command, "37.02", "BBMap wrong version")

    def test_bowtie_version(self):
        command = self.build_command("bowtie2 --version | awk 'NR==1 {print $NF}'")
        self.check_command(command, "2.2.1", "wrong bowtie2 version")

    def test_cufflinks_version(self):
        command = self.build_command("cufflinks 2>&1 | awk 'NR==1'")
        self.check_command(command, "cufflinks v2.2.1", "wrong cufflinks version")

    def test_cutadapt_version(self):
        command = self.build_command("cutadapt --version 2>&1")
        self.check_command(command, "1.8.1", "wrong cutadapt version")

    def test_fastq_version(self):
        command = self.build_command("fastqc --version 2>&1")
        self.check_command(command, "FastQC v0.11.3", "wrong fastqc version")

    def test_fastq_screen_versions(self):
        command = self.build_command("fastq_screen --version 2>&1 || echo 'not installed'")
        self.check_command(command, "v0.11.1", "fastq screen wrong version")

    def test_multiqc_installed(self):
        command = self.build_command("multiqc --version")
        self.check_command(command, "multiqc, version", "multiqc not installed")

    def test_mutt_version(self):
        command = self.build_command("mutt -v | head -1")
        self.check_command(command, "\d*\.\d+", "mutt not installed")

    def test_picard_version(self):
        command = self.build_command("java -jar $PICARD_JARS/SortSam.jar --version 2>&1 | cut -d'(' -f1")
        self.check_command(command, "1.77", "wrong picard version")

    def test_python_version(self):
        command = self.build_command("python --version 2>&1")
        self.check_command(command, "Python 2.7.9", "wrong python version")

    def test_python2_modules_present(self):
        missing_modules = []
        modules = ["HTSeq",
                   "numpy",
                   "xlsxwriter",
                   "yaml"]
        for module_name in modules:
            command = self.build_command("python -c 'import {}'".format(module_name))
            try:
                subprocess.check_output(command, stderr=self.dev_null, shell=True)
            except subprocess.CalledProcessError:
                missing_modules.append(module_name)
        self.assertEquals([], missing_modules)

    def test_R_version(self):
        command = self.build_command("Rscript --version 2>&1")
        self.check_command(command, "version 3.3", "wrong R version")

    def test_R_deseq2_libraries_present(self):
        missing_libs = []
        libs = ['BiocParallel',
                'calibrate',
                'data.table',
                'DESeq2',
                'genefilter',
                'geneplotter',
                'GGally',
                'ggfortify',
                'ggplot2',
                'ggrepel',
                'optparse',
                'pheatmap',
                'plotly',
                'RColorBrewer',
                'reshape2',
                'xlsx',
                ]
        for lib_name in libs:
            rscript = ("Rscript --vanilla -e "
                       "  'result<-1-require({}); "
                       "   quit(status=result)'").format(lib_name)
            command = self.build_command(rscript)
            try:
                subprocess.check_output(command, stderr=self.dev_null, shell=True)
            except subprocess.CalledProcessError:
                missing_libs.append(lib_name)
        self.assertEquals([], missing_libs)

    def test_R_tuxedo_libraries_present(self):
        missing_libs = []
        libs = ['cummeRbund']
        for lib_name in libs:
            rscript = ("Rscript --vanilla -e "
                       "  'result<-1-require({}); "
                       "   quit(status=result)'").format(lib_name)
            command = self.build_command(rscript)
            try:
                subprocess.check_output(command, stderr=self.dev_null, shell=True)
            except subprocess.CalledProcessError:
                missing_libs.append(lib_name)
        self.assertEquals([], missing_libs)

    def test_samtools_version(self):
        command = self.build_command("(samtools 2>&1 | grep 'Version') || echo -e samtools not loaded")
        self.check_command(command, "Version: 1.5", "wrong samtools version")

    def test_tophat_version(self):
        command = self.build_command("tophat --version 2>&1")
        self.check_command(command, "TopHat v2.0.13", "wrong tophat version")
