#pylint: disable=locally-disabled,too-many-public-methods, invalid-name
#pylint: disable=locally-disabled,no-self-use,missing-docstring,protected-access
#pylint: disable=locally-disabled,deprecated-method
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import errno
import functools
import os
import subprocess
import sys
import unittest
try:
    #pylint: disable=locally-disabled,unused-import
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures.tempdirectory import TempDirectory
import pandas as pd
from pandas.util.testing import assert_frame_equal

import yaml

import scripts.watermelon_init as watermelon_init

_TESTS_DIR = os.path.realpath(os.path.dirname(__file__))
_WATERMELON_ROOT = os.path.dirname(_TESTS_DIR)
_CONFIG_DIR = os.path.join(_WATERMELON_ROOT, 'config')
_BIN_DIR = os.path.join(_WATERMELON_ROOT, 'bin')

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)

def _create_files(base_dir, source_files):
    for parent, child in source_files.items():
        if isinstance(child, dict):
            parent_dir = os.path.join(base_dir, parent)
            _mkdir(parent_dir)
            _create_files(parent_dir, child)
        else:
            with open(os.path.join(base_dir, parent), 'wt') as f:
                print(child, file=f)

def _listfiles(directory):
    for root, _, files in os.walk(directory):
        for filename in files:
            yield os.path.join(root, filename)

class WatermelonInitTest(unittest.TestCase):
    def setUp(self):
        self.original_cwd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_cwd)

    def _assertSameFile(self, file1, file2):
        self.assertEqual(os.path.realpath(file1), os.path.realpath(file2))

    def test_build_input_summary_text(self):
        runs_df_contents = StringIO(\
        '''run_index|run
        A|/foo/bar/A
        B|/foo/bar/B
        C|/foo/baz/C
        D|/foo/baz/D''')
        runs_df = pd.read_csv(runs_df_contents, sep='|')
        sample_run_counts_df_contents = StringIO(\
        '''sample/run_index|A|B|C|D
        sample_1|2|4| |
        sample_2| | |8| 16''')
        sample_run_counts_df = pd.read_csv(sample_run_counts_df_contents, sep='|')
        input_summary = Namespace(total_file_count=732,
                                  total_sample_count=78,
                                  total_run_count=4,
                                  runs_df=runs_df,
                                  sample_run_counts_df=sample_run_counts_df)

        actual_text = watermelon_init._build_input_summary_text(input_summary)
        self.assertRegexpMatches(actual_text,
                                 r'Found 732 files for 78 samples across 4 run')
        self.assertRegexpMatches(actual_text,
                                 r'A.*/foo/bar/A')
        self.assertRegexpMatches(actual_text,
                                 r'D.*/foo/baz/D')
        self.assertRegexpMatches(actual_text,
                                 r'sample_1.*2.*4')
        self.assertRegexpMatches(actual_text,
                                 r'sample_2.*8.*16')

    def test_build_postlude(self):
        args = Namespace(input_dir='INPUT_DIR',
                         input_runs_dir='INPUT_RUNS_DIR',
                         input_samples_dir='INPUT_SAMPLES_DIR',
                         analysis_dir='ANALYSIS_DIR',
                         config_file='path/to/CONFIG_FILE',
                         job_suffix='_JOB_SUFFIX',
                         x_working_dir='WORKING_DIR')
        linker_results = 'LINKER_RESULTS\n'
        summary_text = 'SUMMARY_TEXT\n'

        actual_postlude = watermelon_init._build_postlude(args,
                                                          linker_results,
                                                          summary_text)
        self.assertRegexpMatches(actual_postlude,
                                 r'LINKER_RESULTS')
        self.assertRegexpMatches(actual_postlude,
                                 r'SUMMARY_TEXT')
        self.assertRegexpMatches(actual_postlude,
                                 r'ANALYSIS_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'CONFIG_FILE')
        self.assertRegexpMatches(actual_postlude,
                                 r'INPUT_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'INPUT_RUNS_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'INPUT_SAMPLES_DIR')
        self.assertRegexpMatches(actual_postlude,
                                 r'screen -S watermelon_JOB_SUFFIX')
        self.assertRegexpMatches(actual_postlude,
                                 r'watermelon --dry-run -c CONFIG_FILE')

    def test_parse_args(self):
        command_line_args = ('--genome_build mm10 '
                             '--job_suffix _11_01_A '
                             '--x_working_dir /tmp.foo '
                             '--x_template_config /tmp/watermelon/template_config.yaml '
                             '--x_genome_references /tmp/watermelon/genome_references.yaml '
                             'DNASeqCore/Run_1/rhim/Run_1 DNASeqCore/Run_2/rhim/Run_2').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        expected_args = ['analysis_dir',
                         'config_file',
                         'genome_build',
                         'x_genome_references',
                         'tmp_input_dir',
                         'input_dir',
                         'input_runs_dir',
                         'input_samples_dir',
                         'job_suffix',
                         'source_fastq_dirs',
                         'x_template_config',
                         'x_working_dir']
        self.assertEqual(sorted(expected_args), sorted(vars(args)))
        self.assertEqual('mm10', args.genome_build)
        self.assertEqual('_11_01_A', args.job_suffix)
        self.assertEqual(['DNASeqCore/Run_1/rhim/Run_1',
                          'DNASeqCore/Run_2/rhim/Run_2'],
                         args.source_fastq_dirs)
        join = os.path.join
        basedir = '/tmp.foo'
        expected_analysis_filepath = join(basedir, 'analysis_11_01_A')
        self.assertEqual(expected_analysis_filepath, args.analysis_dir)
        self.assertEqual(os.path.join(expected_analysis_filepath,
                                      'config_11_01_A.yaml'),
                         args.config_file)
        self.assertEqual(os.path.join(basedir, '.inputs'),
                         args.tmp_input_dir)
        self.assertEqual(os.path.join(basedir, 'inputs'),
                         args.input_dir)
        self.assertEqual('00-source_runs', args.input_runs_dir)
        self.assertEqual('01-source_samples', args.input_samples_dir)
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

    def test_link_run_dirs(self):
        source_files = {
            'Run_1': {
                'tuttle': {
                    'project': {
                        'sampleA': {
                            '1.fastq':'A',
                            '2.fastq': 'AT'},
                        'sampleB': {
                            '3.fastq':'ATC',
                            '4.fastq': 'ATCG'}
                    }
                }
            },
            'Run_2': {
                'tuttle': {
                    'project': {
                        'sampleA': {
                            '1.fastq':'AA',
                            '2.fastq': 'AATT'},
                        'sampleB': {
                            '3.fastq':'AATTCC',
                            '4.fastq': 'AATTCCGG'}
                    }
                }
            }}
        j = os.path.join
        linker = watermelon_init._Linker()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            orig_wd = os.getcwd()
            try:
                os.chdir(temp_dir_path)
                watermelon_runs_dir = j('00-source_runs', '')
                source_run_basedir = j('DNASeqCore', '')

                _create_files(source_run_basedir, source_files)
                source_run_dirs = [j(source_run_basedir, 'Run_1', 'tuttle', 'project'),
                                   j(source_run_basedir, 'Run_2', 'tuttle', 'project')]

                args = Namespace(source_fastq_dirs=source_run_dirs,
                                 tmp_input_dir=temp_dir_path,
                                 input_runs_dir=watermelon_runs_dir)
                watermelon_init._link_run_dirs(args, linker)

                actual_files = list(_listfiles(j(temp_dir_path, watermelon_runs_dir)))
                actual_link_count = sum([os.stat(f).st_nlink == 2 for f in actual_files])
            finally:
                os.chdir(orig_wd)

        expected_files = ['DNASeqCore^Run_1^tuttle^project/sampleA/1.fastq',
                          'DNASeqCore^Run_1^tuttle^project/sampleA/2.fastq',
                          'DNASeqCore^Run_1^tuttle^project/sampleB/3.fastq',
                          'DNASeqCore^Run_1^tuttle^project/sampleB/4.fastq',
                          'DNASeqCore^Run_2^tuttle^project/sampleA/1.fastq',
                          'DNASeqCore^Run_2^tuttle^project/sampleA/2.fastq',
                          'DNASeqCore^Run_2^tuttle^project/sampleB/3.fastq',
                          'DNASeqCore^Run_2^tuttle^project/sampleB/4.fastq']
        prefix_path = lambda x: temp_dir_path.replace(os.path.sep, '^') + '^' + x
        expected_files = [prefix_path(f) for f in expected_files]
        expected_files = [j(temp_dir_path, watermelon_runs_dir, f) for f in expected_files]
        self.assertEqual(set(expected_files), set(actual_files))
        self.assertEqual(len(actual_files), actual_link_count)

    def test_merge_samples_dirs(self):
        source_runs = {
            'DNASeqCore^Run_1^tuttle^project': {
                'sampleA': {
                    '1.fastq':'A',
                    '2.fastq': 'AT'},
                'sampleB': {
                    '3.fastq':'ATC',
                    '4.fastq': 'ATCG'}
            },
            'DNASeqCore^Run_2^tuttle^project': {
                'sampleA': {
                    '1.fastq':'AA',
                    '2.fastq': 'AATT'},
                'sampleB': {
                    '3.fastq':'AATTCC',
                    '4.fastq': 'AATTCCGG'}
            },
        }
        j = os.path.join
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            orig_wd = os.getcwd()
            try:
                os.chdir(temp_dir_path)
                args = Namespace(tmp_input_dir=j(temp_dir_path, '.input', ''),
                                 input_runs_dir=j('00-source_runs', ''),
                                 input_samples_dir=j('01-source_samples', ''))

                input_runs_dir = j(args.tmp_input_dir, args.input_runs_dir)
                _mkdir(input_runs_dir)
                _create_files(input_runs_dir, source_runs)
                watermelon_init._merge_sample_dirs(args)

                actual_files = list(_listfiles(j(args.tmp_input_dir,
                                                 args.input_samples_dir)))
                actual_link_count = sum([os.stat(f).st_nlink == 2 for f in actual_files])
            finally:
                os.chdir(orig_wd)

        expected_files = ['sampleA/DNASeqCore^Run_1^tuttle^project^1.fastq',
                          'sampleA/DNASeqCore^Run_1^tuttle^project^2.fastq',
                          'sampleA/DNASeqCore^Run_2^tuttle^project^1.fastq',
                          'sampleA/DNASeqCore^Run_2^tuttle^project^2.fastq',
                          'sampleB/DNASeqCore^Run_1^tuttle^project^3.fastq',
                          'sampleB/DNASeqCore^Run_1^tuttle^project^4.fastq',
                          'sampleB/DNASeqCore^Run_2^tuttle^project^3.fastq',
                          'sampleB/DNASeqCore^Run_2^tuttle^project^4.fastq']
        prefix_path = lambda x: j(args.tmp_input_dir, args.input_samples_dir, x)
        expected_files = [prefix_path(f) for f in expected_files]
        self.assertEqual(set(expected_files), set(actual_files))
        self.assertEqual(len(actual_files), actual_link_count)

    def test_build_input_summary(self):
        source_files = {
            '01-source_samples': {
                'sample_1': {
                    'foo^run1^1.fastq.gz':'A',
                    'foo^run2^2.fastq.gz': 'AT',},
                'sample_2': {
                    'foo^run2^2.fastq.gz': 'ATCG',},
                'sample_3': {
                    'foo^run1^1.fastq.gz':'ATC',}
                },
            }
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            _create_files(temp_dir_path, source_files)
            args = Namespace(tmp_input_dir=temp_dir_path,
                             input_samples_dir='01-source_samples')
            actual_summary = watermelon_init._build_input_summary(args)

        actual_keys = vars(actual_summary).keys()
        expected_keys = set(['runs_df',
                             'sample_run_counts_df',
                             'total_sample_count',
                             'total_run_count',
                             'total_file_count',
                             'samples',
                             'samples_missing_fastq_files',
                             'runs_missing_fastq_files'])
        self.assertEqual(expected_keys, actual_keys)
        self.assertEqual(['sample_1', 'sample_2', 'sample_3'],
                         actual_summary.samples)
        self.assertEqual(3, actual_summary.total_sample_count)
        self.assertEqual(2, actual_summary.total_run_count)
        self.assertEqual(4, actual_summary.total_file_count)
        self.assertEqual([], actual_summary.samples_missing_fastq_files)
        self.assertEqual([], actual_summary.runs_missing_fastq_files)

        runs_df_contents = StringIO(\
'''run_index|run
A|foo/run1
B|foo/run2''')
        expected_runs_df = pd.read_csv(runs_df_contents, sep='|', index_col=0)
        assert_frame_equal(expected_runs_df,
                           actual_summary.runs_df)

        sample_run_counts_df_contents = StringIO(\
'''sample|A|B|sample_total
sample_1|1|1|2
sample_2|0|1|1
sample_3|1|0|1
run_total|2|2|4''')
        expected_df = pd.read_csv(sample_run_counts_df_contents, sep='|', index_col=0)
        expected_df.columns.names=['run_index']

        assert_frame_equal(expected_df, actual_summary.sample_run_counts_df)

    def test_build_input_summary_samplesMissingFastqFiles(self):
        source_files = {
            '01-source_samples': {
                'sample_1': {
                    'foo^run1^1.fastq.gz':'A',},
                'sample_2': {
                    'foo^run2^2.fastq': 'ATCG',},
                'sample_3': {
                    'foo^run1^1.foo':'ATC',},
                'sample_4': {
                    'foo^run1^2.':'',},
                },
            }
        j = os.path.join
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            _create_files(temp_dir_path, source_files)

            args = Namespace(tmp_input_dir=temp_dir_path,
                             input_samples_dir='01-source_samples')
            actual_summary = watermelon_init._build_input_summary(args)

        sample_run_counts_df_contents = StringIO(\
'''sample|A|B|sample_total
sample_1|1|0|1
sample_2|0|1|1
sample_3|0|0|0
sample_4|0|0|0
run_total|1|1|2''')
        expected_df = pd.read_csv(sample_run_counts_df_contents, sep='|', index_col=0)
        expected_df.columns.names=['run_index']

        assert_frame_equal(expected_df, actual_summary.sample_run_counts_df)
        self.assertEqual(['sample_3', 'sample_4'],
                         actual_summary.samples_missing_fastq_files)

    def test_build_input_summary_runsMissingFastqFiles(self):
        source_files = {
            '01-source_samples': {
                'sample_1': {
                    'foo^run1^1.fastq.gz':'A',
                    'foo^run2^1.foo':'A',
                    'foo^run3^1.foo':'A',},
                'sample_2': {
                    'foo^run1^1.fastq': 'ATCG',
                    'foo^run2^1.foo': 'ATCG',
                    'foo^run3^1.foo': 'ATCG',},
                },
            }
        j = os.path.join
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            _create_files(temp_dir_path, source_files)
            args = Namespace(tmp_input_dir=temp_dir_path,
                             input_samples_dir='01-source_samples')
            actual_summary = watermelon_init._build_input_summary(args)

        sample_run_counts_df_contents = StringIO(\
'''sample|A|B|C|sample_total
sample_1|1|0|0|1
sample_2|1|0|0|1
run_total|2|0|0|2''')
        expected_df = pd.read_csv(sample_run_counts_df_contents, sep='|', index_col=0)
        expected_df.columns.names=['run_index']

        assert_frame_equal(expected_df, actual_summary.sample_run_counts_df)

        self.assertEqual(['B','C'], actual_summary.runs_missing_fastq_files)

    def test_make_config_dict(self):
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
    gtf: /ccmb/BioinfCore/hg19_noRibo.gtf
    bowtie2_index: /ccmb/BioinfCore/iGenomes/hg19/Sequence/Bowtie2Index
    entrez_gene_info: /ccmb/BioinfCore/entrez_gene_info/2016_09_02/gene_info
''')
        args = Namespace(input_dir = '/my/input/dir',
                         input_samples_dir='samples')
        samples = ['sA', 'sB', 'sC', 'sD', 'sE', 'sF']
        actual_config = watermelon_init._make_config_dict(template_config,
                                                          genome_references,
                                                          args,
                                                          samples)
        expected_keys = ['dirs',
                         'foo1', 'foo2',
                         'genome', 'references',
                         'main_factors', 'phenotypes', 'samples', 'comparisons']
        self.assertEquals(sorted(expected_keys), sorted(actual_config.keys()))
        self.assertEqual('/my/input/dir/samples', actual_config['dirs']['input'])
        self.assertEqual(template_config['foo1'], actual_config['foo1'])
        self.assertEqual(template_config['foo2'], actual_config['foo2'])
        self.assertEqual(genome_references['genome'], actual_config['genome'])
        self.assertEqual(genome_references['references'], actual_config['references'])
        self.assertEqual('yes    ^ yes      ^ no', actual_config['main_factors'])
        self.assertEqual('gender ^ genotype ^ gender.genotype', actual_config['phenotypes'])
        self.maxDiff = None
        self.assertEqual({'sA':'female ^ MutA     ^ female.MutA',
                          'sB':'female ^ MutB     ^ female.MutB',
                          'sC':'female ^ WT       ^ female.WT',
                          'sD':'male   ^ MutA     ^ male.MutA',
                          'sE':'male   ^ MutB     ^ male.MutB',
                          'sF':'male   ^ WT       ^ male.WT',
                         },
                         actual_config['samples'])
        self.assertEqual(['male_v_female'],
                         actual_config['comparisons']['gender'])
        self.assertEqual(['MutA_v_WT', 'MutB_v_WT'],
                         actual_config['comparisons']['genotype'])

    def test_make_config_dict_intelligentlyMergesGenomeReferences(self):
        template_config = yaml.load(\
'''fastq_screen:
    aligner: bowtie2
''')
        genome_references = yaml.load(\
'''fastq_screen:
    species: human
''')
        args = Namespace(input_dir='/my/input/dir',
                         input_samples_dir='samples')
        samples = []
        actual_config = watermelon_init._make_config_dict(template_config,
                                                          genome_references,
                                                          args,
                                                          samples)
        expected_config = {'aligner': 'bowtie2', 'species': 'human'}
        self.assertEqual(expected_config, actual_config['fastq_screen'])


    def test_write_config_file(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_filename = os.path.join(temp_dir_path, 'config.yaml')
            config_dict = {'dirs': {'input' : 'INPUT_DIR'},
                           'samples' : 'SAMPLES',
                           'comparisons' : 'COMPARISONS',
                           'genome' : 'GENOME',
                           'references' : 'REFERENCES',
                           'templateA' : {'TEMPLATE_A1': 'A1', 'TEMPLATE_A2': 'A2'},
                           'templateB' : 'TEMPLATE_B',
                           'templateC' : 'TEMPLATE_C',
                           'phenotypes' :   'PHENOLABEL1 ^ PHENOLABEL2',
                           'main_factors' : 'yes ^ no'}
            watermelon_init._write_config_file(config_filename, config_dict)

            with open(config_filename, 'r') as config_file:
                config_lines = [line.strip('\n') for line in config_file.readlines()]
        config_prelude = watermelon_init._CONFIG_PRELUDE.split('\n')
        expected_config_lines = 13
        self.assertEqual(expected_config_lines + len(config_prelude),
                         len(config_lines))
        line_iter = iter(config_lines)
        for _ in config_prelude:
            next(line_iter)
        self.assertEqual('dirs:', next(line_iter))
        self.assertEqual('    input: INPUT_DIR', next(line_iter))
        self.assertEqual('main_factors: yes ^ no', next(line_iter))
        self.assertEqual('phenotypes: PHENOLABEL1 ^ PHENOLABEL2', next(line_iter))
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
    def setUp(self):
        self.stderr_saved = sys.stderr
        self.stderr = StringIO()
        sys.stderr = self.stderr

    def tearDown(self):
        sys.stderr = self.stderr_saved

    def ok(self):
        #pylint: disable=locally-disabled,redundant-unittest-assert
        self.assertTrue(True)

    def test_validate_source_fastq_dirs_raisesIfInvalid(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_fastq_dir_a = os.path.join(temp_dir_path, 'dir_a')
            source_fastq_dir_b = os.path.join(temp_dir_path, 'dir_b')
            source_fastq_dir_c = os.path.join(temp_dir_path, 'dir_c')
            os.mkdir(source_fastq_dir_b)
            args = Namespace(source_fastq_dirs=[source_fastq_dir_a,
                                                source_fastq_dir_b,
                                                source_fastq_dir_c])
            msg = (r'Specified source_fastq_dir\(s\) \[.*dir_a,.*dir_c\] not a '
                   r'dir or cannot be read.')
            self.assertRaisesRegexp(watermelon_init._UsageError,
                                    msg,
                                    validator._validate_source_fastq_dirs,
                                    args)

    def test_validate_source_fastq_dir_ok(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            os.mkdir(os.path.join(temp_dir_path, 'source'))
            source_fastq_dir = os.path.join(temp_dir_path, 'source')
            args = Namespace(source_fastq_dirs=[source_fastq_dir])
            validator._validate_source_fastq_dirs(args)
            self.ok()

    def test_validate_overwrite_check_ok(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            analysis_dir = os.path.join(temp_dir_path, 'analysis')
            input_dir = os.path.join(temp_dir_path, 'inputs')
            args = Namespace(analysis_dir=analysis_dir,
                             input_dir=input_dir)
            validator._validate_overwrite_check(args)
            self.ok()

    def test_validate_overwrite_check_raisesIfAnalysisExists(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            analysis_dir = os.path.join(temp_dir_path, 'analysis')
            _mkdir(analysis_dir)
            input_dir = os.path.join(temp_dir_path, 'inputs')
            args = Namespace(analysis_dir=analysis_dir,
                             input_dir=input_dir)
            self.assertRaisesRegexp(watermelon_init._UsageError,
                                    r'dirs: analysis \[.*\] exists',
                                    validator._validate_overwrite_check,
                                    args)

    def test_validate_overwrite_check_raisesIfInputDirExists(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            input_dir = os.path.join(temp_dir_path, 'inputs')
            _mkdir(os.path.join(input_dir))
            analysis_dir = os.path.join(temp_dir_path, 'analysis')
            args = Namespace(analysis_dir=analysis_dir,
                             input_dir=input_dir)
            self.assertRaisesRegexp(watermelon_init._UsageError,
                                    r'dirs: input \[.*\] exists',
                                    validator._validate_overwrite_check,
                                    args)

    def test_validate_runs_have_fastq_files_raisesIfInvalid(self):
        validator = watermelon_init._CommandValidator()
        runs_df_contents = StringIO(\
'''run_index|run
A|/foo/bar/A
B|/foo/bar/B
C|/foo/bar/C''')
        runs_df = pd.read_csv(runs_df_contents, sep='|', index_col='run_index')
        sample_run_counts_df_contents = StringIO(\
'''sample|A|B|C
sample_1|1|0|0
sample_2|1|0|0
run_total|2|0|0''')
        sample_run_counts_df = pd.read_csv(sample_run_counts_df_contents,
                                           sep='|',
                                           index_col='sample')
        input_summary = Namespace(total_file_count=732,
                                  total_sample_count=78,
                                  total_run_count=4,
                                  runs_df=runs_df,
                                  sample_run_counts_df=sample_run_counts_df,
                                  runs_missing_fastq_files=['B', 'C'])

        msg = (r'Some runs missing fastq files: \[B, C\].')
        self.assertRaisesRegexp(watermelon_init._InputValidationError,
                                msg,
                                validator._validate_runs_have_fastq_files,
                                input_summary)
        self.assertRegexpMatches(self.stderr.getvalue(), r'Input summary')

    def test_validate_runs_have_fastq_files_ok(self):
        validator = watermelon_init._CommandValidator()
        input_summary = Namespace(runs_missing_fastq_files=[])
        validator._validate_runs_have_fastq_files(input_summary)
        self.assertEqual(self.stderr.getvalue(), '')

    def test_validate_samples_have_fastq_files_raisesIfInvalid(self):
        validator = watermelon_init._CommandValidator()
        runs_df_contents = StringIO(\
'''run_index|run
A|/foo/bar/A''')
        runs_df = pd.read_csv(runs_df_contents, sep='|', index_col='run_index')
        sample_run_counts_df_contents = StringIO(\
'''sample|A|sample_total
sample_1|1|1
sample_2|0|0
sample_3|0|0
run_total|1|1''')
        sample_run_counts_df = pd.read_csv(sample_run_counts_df_contents,
                                           sep='|',
                                           index_col='sample')
        input_summary = Namespace(total_file_count=732,
                                  total_sample_count=78,
                                  total_run_count=4,
                                  runs_df=runs_df,
                                  sample_run_counts_df=sample_run_counts_df,
                                  samples_missing_fastq_files=['sample_2', 'sample_3'])

        msg = (r'Some samples missing fastq files: \[sample_2, sample_3\].')
        self.assertRaisesRegexp(watermelon_init._InputValidationError,
                                msg,
                                validator._validate_samples_have_fastq_files,
                                input_summary)
        self.assertRegexpMatches(self.stderr.getvalue(), r'Input summary')

    def test_validate_samples_have_fastq_files_ok(self):
        validator = watermelon_init._CommandValidator()
        input_summary = Namespace(samples_missing_fastq_files=[])
        validator._validate_samples_have_fastq_files(input_summary)
        self.assertEqual(self.stderr.getvalue(), '')


class LinkerTest(unittest.TestCase):
    def test_link_allHardlinks(self):
        linker = watermelon_init._Linker()
        j = os.path.join
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_dir = j(temp_dir_path, 'source')
            _mkdir(source_dir)
            dest_dir = j(temp_dir_path, 'dest')
            _mkdir(dest_dir)
            source_filename1 = os.path.join(source_dir, '1.tmp')
            source_filename2 = os.path.join(source_dir, '2.tmp')
            touch(source_filename1)
            touch(source_filename2)
            dest_filename1 = os.path.join(dest_dir, '1.lnk')
            dest_filename2 = os.path.join(dest_dir, '2.lnk')
            linker._link(source_filename1, dest_filename1)
            linker._link(source_filename2, dest_filename2)

            dest_filename1_stat = os.stat(dest_filename1)
            dest_filename2_stat = os.stat(dest_filename2)

            self.assertEqual(os.stat(source_filename1).st_ino,
                              dest_filename1_stat.st_ino)
            self.assertEqual(os.stat(source_filename2).st_ino,
                              dest_filename2_stat.st_ino)

            self.assertEqual(dest_filename1_stat.st_nlink, 2)
            self.assertEqual(dest_filename2_stat.st_nlink, 2)

        self.assertEquals('linked 2 source files: 2 hardlinked', linker.results)

    def test_link_allSymlinks(self):
        linker = watermelon_init._Linker()
        j = os.path.join
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_dir = j(temp_dir_path, 'source')
            _mkdir(source_dir)
            dest_dir = j(temp_dir_path, 'dest')
            _mkdir(dest_dir)
            source_filename1 = os.path.join(source_dir, '1.tmp')
            source_filename2 = os.path.join(source_dir, '2.tmp')
            touch(source_filename1)
            touch(source_filename2)
            dest_filename1 = os.path.join(dest_dir, '1.lnk')
            dest_filename2 = os.path.join(dest_dir, '2.lnk')
            def angry_link(source, dest):
                raise OSError(errno.EXDEV, 'Cross-device link')
            try:
                os_link = os.link
                os.link = angry_link
                linker._link(source_filename1, dest_filename1)
                linker._link(source_filename2, dest_filename2)
            finally:
                os.link = os_link

            dest_filename1_stat = os.stat(dest_filename1, follow_symlinks=False)
            dest_filename2_stat = os.stat(dest_filename2, follow_symlinks=False)

            self.assertNotEqual(os.stat(source_filename1).st_ino,
                              dest_filename1_stat.st_ino)
            self.assertNotEqual(os.stat(source_filename2).st_ino,
                              dest_filename2_stat.st_ino)
        self.assertEquals('linked 2 source files: 2 symlinked', linker.results)

    def test_link_hardlinksAndSymlinks(self):
        linker = watermelon_init._Linker()
        j = os.path.join
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_dir = j(temp_dir_path, 'source')
            _mkdir(source_dir)
            dest_dir = j(temp_dir_path, 'dest')
            _mkdir(dest_dir)
            source_filename1 = os.path.join(source_dir, '1.tmp')
            source_filename2 = os.path.join(source_dir, '2.tmp')
            touch(source_filename1)
            touch(source_filename2)
            dest_filename1 = os.path.join(dest_dir, '1.lnk')
            dest_filename2 = os.path.join(dest_dir, '2.lnk')
            linker._link(source_filename1, dest_filename1)
            def angry_link(source, dest):
                raise OSError(errno.EXDEV, 'Cross-device link')
            try:
                os_link = os.link
                os.link = angry_link
                linker._link(source_filename2, dest_filename2)
            finally:
                os.link = os_link

            dest_filename1_stat = os.stat(dest_filename1, follow_symlinks=False)
            dest_filename2_stat = os.stat(dest_filename2, follow_symlinks=False)

            self.assertEqual(os.stat(source_filename1).st_ino,
                              dest_filename1_stat.st_ino)
            self.assertNotEqual(os.stat(source_filename2).st_ino,
                              dest_filename2_stat.st_ino)
        self.assertEquals('linked 2 source files: 1 hardlinked, 1 symlinked', linker.results)



class WatermelonInitFunctoinalTest(unittest.TestCase):
    def execute(self, command):
        exit_code = 0
        try:
            actual_output = subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as e:
            exit_code = e.returncode
            actual_output = e.output
        return exit_code, str(actual_output)

    def test_commandPresentUsageIfMissingOptions(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path

            script_name = os.path.join(_BIN_DIR, 'watermelon_init')

            redirect_output = ' 2>&1 '
            command = '{} {}'.format(script_name, redirect_output)
            exit_code, command_output = self.execute(command)

        self.assertEqual(2, exit_code, command)
        self.assertRegexpMatches(command_output, 'usage',command)

    def test_commandErrorUsageIfInvalidInput(self):
        source_files = {
            'Run_1': {
                'tuttle': {
                    'project': {
                        'sampleA': {
                            '1.fastq':'A',
                            '2.fastq': 'AT'},
                        'sampleB': {
                            'foo.txt': 'foo',},
                        },
                    },
                },
            }
        j = os.path.join
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            orig_wd = os.getcwd()
            try:
                os.chdir(temp_dir_path)
                source_run_basedir = j('DNASeqCore', '')
                _mkdir(source_run_basedir)
                _create_files(source_run_basedir, source_files)
                source_run_dirs = [j(source_run_basedir, 'Run_1', 'tuttle', 'project'),]

                script_name = os.path.join(_BIN_DIR, 'watermelon_init')

                redirect_output = ' 2>&1 '
                command = '{} {} {} {}'.format(script_name,
                                               '--genome hg19',
                                               ' '.join(source_run_dirs),
                                               redirect_output)
                print(command)
                exit_code, command_output = self.execute(command)
            finally:
                os.chdir(orig_wd)
        self.assertEqual(1, exit_code, command)
        self.assertRegexpMatches(command_output,
                                 'Some samples missing fastq files')
