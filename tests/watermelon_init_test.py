#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import functools
import os
import os.path
import sys
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures.tempdirectory import TempDirectory
import yaml

import scripts.watermelon_init as watermelon_init

_TESTS_DIR = os.path.realpath(os.path.dirname(__file__))
_WATERMELON_ROOT = os.path.dirname(_TESTS_DIR)
_CONFIG_DIR = os.path.join(_WATERMELON_ROOT, 'config')

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

class WatermelonInitTest(unittest.TestCase):
    def test_build_postlude(self):
        args = Namespace(inputs_dir='INPUTS_DIR',
                         source_fastq_dir='SOURCE_FASTQ_DIR',
                         analysis_dir='ANALYSIS_DIR',
                         config_file='CONFIG_FILE',
                         job_suffix='_JOB_SUFFIX')
        sample_count = 42
        file_count = 168
        actual_postlude = watermelon_init._build_postlude(args, sample_count, file_count)
        self.assertRegexpMatches(actual_postlude,
                                 r'INPUTS_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'SOURCE_FASTQ_DIR \| 42 samples \| 168 files')
        self.assertRegexpMatches(actual_postlude,
                                 r'ANALYSIS_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'screen -S watermelon_JOB_SUFFIX')
        self.assertRegexpMatches(actual_postlude,
                                 r'watermelon --dry-run CONFIG_FILE')

    def test_parse_args(self):
        command_line_args = ('--genome_build mm10 '
                             '--job_suffix _11_01_A '
                             '--x_working_dir /tmp.foo '
                             '--x_template_config /tmp/watermelon/template_config.yaml '
                             '--x_genome_references /tmp/watermelon/genome_references.yaml '
                             'DNASeqCore/Run_1286/rhim/Run_1286').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        expected_args = ['analysis_dir',
                         'config_file',
                         'genome_build',
                         'x_genome_references',
                         'inputs_dir',
                         'job_suffix',
                         'source_fastq_dir',
                         'x_template_config',
                         'x_working_dir']
        self.assertEqual(sorted(expected_args), sorted(vars(args)))
        self.assertEqual('mm10', args.genome_build)
        self.assertEqual('_11_01_A', args.job_suffix)
        self.assertEqual('DNASeqCore/Run_1286/rhim/Run_1286', args.source_fastq_dir)
        join = os.path.join
        basedir = '/tmp.foo'
        expected_analysis_filepath = join(basedir, 'analysis_11_01_A')
        self.assertEqual(expected_analysis_filepath, args.analysis_dir)
        self.assertEqual(os.path.join(expected_analysis_filepath, 'config_11_01_A.yaml'),
                         args.config_file)
        self.assertEqual(os.path.join(basedir, 'inputs', '00-multiplexed_reads'),
                         args.inputs_dir)
        self.assertEqual('/tmp/watermelon/template_config.yaml', args.x_template_config)
        self.assertEqual('/tmp/watermelon/genome_references.yaml', args.x_genome_references)

    def test_parse_args_workingDirDefaultsToCwd(self):
        command_line_args = ('--genome_build mm10 '
                             '--job_suffix _11_01_A '
                             'DNASeqCore/Run_1286/rhim/Run_1286').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        self.assertEqual(os.getcwd(), args.x_working_dir) 

    def test_parse_args_templateConfigDefaultsToWatermelonConfigDir(self):
        command_line_args = ('--genome_build mm10 '
                             '--job_suffix _11_01_A '
                             'DNASeqCore/Run_1286/rhim/Run_1286').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        self.assertEqual(os.path.join(_CONFIG_DIR, 'template_config.yaml'),
                         args.x_template_config) 

    def test_parse_args_genomeReferencesDefaultsToWatermelonConfigDir(self):
        command_line_args = ('--genome_build mm10 '
                             '--job_suffix _11_01_A '
                             'DNASeqCore/Run_1286/rhim/Run_1286').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        self.assertEqual(os.path.join(_CONFIG_DIR, 'genome_references.yaml'),
                         args.x_genome_references) 

    def test_GENOME_BUILD_OPTIONS_matchGenomeReferenceKeys(self):
        with open(os.path.join(_CONFIG_DIR, 'genome_references.yaml'), 'r') as yaml_file:
            genome_references = yaml.load(yaml_file)
        self.assertEqual(sorted(watermelon_init.GENOME_BUILD_OPTIONS),
                         sorted(genome_references.keys()))

    def test_make_top_level_dirs(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            real = functools.partial(os.path.join, temp_dir_path)
            args = Namespace(analysis_dir=real('ANALYSIS_11_10_A'),
                             inputs_dir=real('INPUTS'))
            watermelon_init._make_top_level_dirs(args)

            def is_dir(o):
                return os.path.isdir(os.path.join(temp_dir_path,o))
            actual_dirs = sorted([o for o in os.listdir(temp_dir_path) if is_dir(o)])
            self.assertEqual(['ANALYSIS_11_10_A','INPUTS'],
                             actual_dirs)

    def test_populate_input_dir(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_fastq_dir = os.path.join(temp_dir_path, 'source_fastq_dir')
            os.mkdir(source_fastq_dir)
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_A'))
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_B'))
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_C'))

            os.mkdir(os.path.join(temp_dir_path, 'working_dir'))
            inputs_dir = os.path.join(temp_dir_path, 'working_dir', 'inputs', '00-multiplexed')
            os.mkdir(os.path.join(temp_dir_path, 'working_dir', 'inputs'))
            os.mkdir(inputs_dir)
            samples = {'sA' : os.path.join(source_fastq_dir, 'Sample_A'),
                       'sB' : os.path.join(source_fastq_dir, 'Sample_B'),
                       'sC' : os.path.join(source_fastq_dir, 'Sample_C'),}

            watermelon_init._populate_inputs_dir(inputs_dir, samples)

            def is_link(o):
                return os.path.islink(os.path.join(temp_dir_path, inputs_dir, o))
            actual_links = sorted([o for o in os.listdir(inputs_dir) if is_link(o)])
            self.assertEqual(['sA', 'sB' ,'sC'], actual_links)
            self.assertEqual(os.path.join(source_fastq_dir, 'Sample_A'),
                             os.readlink(os.path.join(inputs_dir, 'sA')))
            self.assertEqual(os.path.join(source_fastq_dir, 'Sample_B'),
                             os.readlink(os.path.join(inputs_dir, 'sB')))
            self.assertEqual(os.path.join(source_fastq_dir, 'Sample_C'),
                             os.readlink(os.path.join(inputs_dir, 'sC')))

    def _assertSameFile(self, file1, file2):
        self.assertEqual(os.path.realpath(file1), os.path.realpath(file2))

    def test_initialize_samples(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_fastq_dir = os.path.join(temp_dir_path, 'source_fastq_dir')
            os.mkdir(source_fastq_dir)
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_A'))
            touch(os.path.join(source_fastq_dir, 'Sample_A', 'a1.fastq.gz'))
            touch(os.path.join(source_fastq_dir, 'Sample_A', 'a2.fastq.gz'))
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_B'))
            touch(os.path.join(source_fastq_dir, 'Sample_B', 'b.fastq.gz'))
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_C'))
            touch(os.path.join(source_fastq_dir, 'Sample_C', 'c.fastq.gz'))
            (actual_samples,
             file_count) = watermelon_init._initialize_samples(source_fastq_dir)
        self.assertEqual(['Sample_A', 'Sample_B', 'Sample_C'],
                         sorted(actual_samples.keys()))
        self._assertSameFile(os.path.join(source_fastq_dir, 'Sample_A'),
                             actual_samples['Sample_A'])
        self._assertSameFile(os.path.join(source_fastq_dir, 'Sample_B'),
                             actual_samples['Sample_B'])
        self._assertSameFile(os.path.join(source_fastq_dir, 'Sample_C'),
                             actual_samples['Sample_C'])
        self.assertEqual(4, file_count)

    def test_initialize_samples_missingFastqExcluded(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_fastq_dir = os.path.join(temp_dir_path, 'source_fastq_dir')
            os.mkdir(source_fastq_dir)
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_A'))
            touch(os.path.join(source_fastq_dir, 'Sample_A', 'a.fastq.gz'))
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_B'))
            touch(os.path.join(source_fastq_dir, 'Sample_B', 'b.fastq'))
            os.mkdir(os.path.join(source_fastq_dir, 'Sample_C'))
            touch(os.path.join(source_fastq_dir, 'Sample_C', 'c.foo'))
            (actual_samples,
             file_count) = watermelon_init._initialize_samples(source_fastq_dir)
        self.assertEqual(['Sample_A', 'Sample_B'],
                         sorted(actual_samples.keys()))
        self._assertSameFile(os.path.join(source_fastq_dir, 'Sample_A'),
                             actual_samples['Sample_A'])
        self._assertSameFile(os.path.join(source_fastq_dir, 'Sample_B'),
                             actual_samples['Sample_B'])
        self.assertEqual(2, file_count)

    def test_initialize_samples_emptyIfNoFastqFound(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_fastq_dir = os.path.join(temp_dir_path, 'source_fastq_dir')
            os.mkdir(source_fastq_dir)
            (actual_samples,
             file_count) = watermelon_init._initialize_samples(source_fastq_dir)
        self.assertEqual({}, actual_samples)
        self.assertEqual(0, file_count)

    def test_make_config(self):
        template_config = yaml.load(\
'''foo1:
    bar1: baz1
    hoopy1: frood1
foo2:
    bar2: baz2
''')
        genome_references = yaml.load(\
'''genome:
    hg19

references:
    gtf: /ccmb/BioinfCore/noRibo_human.gtf
    bowtie2_index: /ccmb/BioinfCore/iGenomes/hg19/Sequence/Bowtie2Index
    entrez_gene_info: /ccmb/BioinfCore/entrez_gene_info/2016_09_02/gene_info
''')
        input_dir = '/my/input/dir'
        samples = {'sA':'fooA', 'sB': 'fooB', 'sC': 'fooC'}
        actual_config = watermelon_init._make_config_dict(template_config,
                                                          genome_references,
                                                          input_dir,
                                                          samples)
        expected_keys = ['input_dir',
                         'foo1', 'foo2',
                         'genome', 'references',
                         'samples', 'comparisons']
        self.assertEquals(sorted(expected_keys), sorted(actual_config.keys()))
        self.assertEqual('"' + input_dir + '"', actual_config['input_dir'])
        self.assertEqual(template_config['foo1'], actual_config['foo1'])
        self.assertEqual(template_config['foo2'], actual_config['foo2'])
        self.assertEqual(genome_references['genome'], actual_config['genome'])
        self.assertEqual(genome_references['references'], actual_config['references'])
        self.assertEqual({'sA':'g1', 'sB':'g1', 'sC':'g1'},
                         actual_config['samples'])
        self.assertEqual({'g1_v_g2':'g1_v_g2'}, actual_config['comparisons'])

    def test_write_config_file(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_filename = os.path.join(temp_dir_path, 'config.yaml')
            config_dict = {'input_dir' : 'INPUT_DIR',
                           'samples' : 'SAMPLES',
                           'comparisons' : 'COMPARISONS',
                           'genome' : 'GENOME',
                           'references' : 'REFERENCES',
                           'templateA' : {'TEMPLATE_A1': 'A1', 'TEMPLATE_A2': 'A2'},
                           'templateB' : 'TEMPLATE_B',
                           'templateC' : 'TEMPLATE_C'}
            watermelon_init._write_config_file(config_filename, config_dict)

            with open(config_filename, 'r') as config_file:
                config_lines = [line.strip('\n') for line in config_file.readlines()]
        print(config_lines)
        self.assertEqual(10, len(config_lines))
        line_iter = iter(config_lines)
        self.assertEqual('input_dir: INPUT_DIR', next(line_iter))
        self.assertEqual('samples: SAMPLES', next(line_iter))
        self.assertEqual('comparisons: COMPARISONS', next(line_iter))
        self.assertEqual('genome: GENOME', next(line_iter))
        self.assertEqual('references: REFERENCES', next(line_iter))
        self.assertEqual('templateA:', next(line_iter))
        self.assertEqual('    TEMPLATE_A1: A1', next(line_iter))
        self.assertEqual('    TEMPLATE_A2: A2', next(line_iter))
        self.assertEqual('templateB: TEMPLATE_B', next(line_iter))
        self.assertEqual('templateC: TEMPLATE_C', next(line_iter))

class CommandValidatorTest(unittest.TestCase):
    def ok(self):
        self.assertTrue(True)

    def test_validate_source_fastq_dir_raisesIfInvalid(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_fastq_dir = os.path.join(temp_dir_path, 'source')
            args = Namespace(source_fastq_dir=source_fastq_dir)
            self.assertRaisesRegex(watermelon_init.UsageError,
                                   'Specified source_fastq_dir \[.*\] is not a dir or cannot be read. ',
                                   validator._validate_source_fastq_dir,
                                   args)

    def test_validate_source_fastq_dir_ok(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            os.mkdir(os.path.join(temp_dir_path, 'source'))
            source_fastq_dir = os.path.join(temp_dir_path, 'source')
            args = Namespace(source_fastq_dir=source_fastq_dir)
            validator._validate_source_fastq_dir(args)
            self.ok()

    def test_validate_overwrite_check_ok(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            analysis_dir = os.path.join(temp_dir_path, 'analysis')
            inputs_dir = os.path.join(temp_dir_path, 'inputs')
            args = Namespace(analysis_dir=analysis_dir,
                             inputs_dir=inputs_dir)
            validator._validate_overwrite_check(args)
            self.ok()

    def test_validate_overwrite_check_raisesIfAnalysisExists(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            os.mkdir(os.path.join(temp_dir_path, 'analysis'))
            analysis_dir = os.path.join(temp_dir_path, 'analysis')
            inputs_dir = os.path.join(temp_dir_path, 'inputs')
            args = Namespace(analysis_dir=analysis_dir, inputs_dir=inputs_dir)
            self.assertRaisesRegex(watermelon_init.UsageError,
                                   r'analysis_dir \[.*\] exist\(s\)',
                                   validator._validate_overwrite_check,
                                   args)

    def test_validate_overwrite_check_raisesIfInputsExists(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            os.mkdir(os.path.join(temp_dir_path, 'inputs'))
            inputs_dir = os.path.join(temp_dir_path, 'inputs')
            analysis_dir = os.path.join(temp_dir_path, 'analysis')
            args = Namespace(analysis_dir=analysis_dir, inputs_dir=inputs_dir)
            self.assertRaisesRegex(watermelon_init.UsageError,
                                   r'inputs_dir \[.*\] exist\(s\)',
                                   validator._validate_overwrite_check,
                                   args)


# class WatermelonInitFunctoinalTest(unittest.TestCase):
#     def execute(self, command):
#         exit_code = 0
#         try:
#             actual_output = subprocess.check_output(command, shell=True)
#         except subprocess.CalledProcessError as e:
#             exit_code = e.returncode
#             actual_output = e.output
#         return exit_code, str(actual_output)
# 
#     def test_commandReturnsCorrectRowAndColumnCount(self):
#         with TempDirectory() as temp_dir:
#             temp_dir_path = temp_dir.path
#             os.
#             script_name = os.path.join(SCRIPTS_DIR, 'watermelon_init.py')
# 
#             redirect_output = '2>/dev/null'
#             command = 'python {} {}'.format(script_name, redirect_output)
#             exit_code, command_output = self.execute(command)
# 
#         self.assertEqual(0, exit_code, command)
#         self.assertEqual('foo', command_output)
