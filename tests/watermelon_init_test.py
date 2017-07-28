#pylint: disable=locally-disabled,too-many-public-methods, invalid-name
#pylint: disable=locally-disabled,no-self-use,missing-docstring,protected-access
#pylint: disable=locally-disabled,deprecated-method
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import functools
import os
import subprocess
import unittest
try:
    #pylint: disable=locally-disabled,unused-import
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures.tempdirectory import TempDirectory
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

    def test_is_source_fastq_external_TrueForSibling(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            parent_dir = temp_dir_path
            working_dir = os.path.join(parent_dir, 'working_dir')
            sibling_dir = os.path.join(parent_dir, 'sibling_dir')
            self.assertTrue(watermelon_init._is_source_fastq_external(sibling_dir,
                                                                      working_dir=working_dir))
    def test_is_source_fastq_external_TrueForParent(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            parent_dir = temp_dir_path
            working_dir = os.path.join(parent_dir, 'working_dir')
            self.assertTrue(watermelon_init._is_source_fastq_external(parent_dir,
                                                                      working_dir=working_dir))

    def test_is_source_fastq_external_FalseForChild(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            parent_dir = temp_dir_path
            working_dir = os.path.join(parent_dir, 'working_dir')
            child_dir = os.path.join(working_dir, 'child_dir')
            self.assertFalse(watermelon_init._is_source_fastq_external(child_dir,
                                                                       working_dir=working_dir))

    def test_build_postlude(self):
        args = Namespace(inputs_dir='INPUTS_DIR',
                         source_fastq_dir='SOURCE_FASTQ_DIR',
                         analysis_dir='ANALYSIS_DIR',
                         config_file='path/to/CONFIG_FILE',
                         job_suffix='_JOB_SUFFIX',
                         x_working_dir='WORKING_DIR')
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
                                 r'    CONFIG_FILE')
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
                             'DNASeqCore/Run_1286/rhim/Run_1286').split(' ')
        args = watermelon_init._parse_command_line_args(command_line_args)
        expected_args = ['analysis_dir',
                         'config_file',
                         'genome_build',
                         'x_genome_references',
                         'inputs_dir',
                         'is_source_fastq_external',
                         'job_suffix',
                         'source_fastq_dir',
                         'x_template_config',
                         'x_working_dir']
        self.assertEqual(sorted(expected_args), sorted(vars(args)))
        self.assertEqual('mm10', args.genome_build)
        self.assertEqual('_11_01_A', args.job_suffix)
        self.assertEqual(True, args.is_source_fastq_external)
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

    def test_parse_args_inputWhenLocalSourceFastq(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            os.chdir(temp_dir_path)
            source_fastq_dir = "source_fastq"
            _mkdir(os.path.join(temp_dir_path, "source_fastq"))
            args = ('--genome_build mm10 '
                    '--job_suffix _11_01_A '
                    '--x_working_dir {working_dir} '
                    '--x_template_config /tmp/watermelon/template_config.yaml '
                    '--x_genome_references /tmp/watermelon/genome_references.yaml '
                    '{source_fastq_dir}')
            command_line_args = args.format(working_dir=temp_dir_path,
                                            source_fastq_dir=source_fastq_dir).split(' ')
            args = watermelon_init._parse_command_line_args(command_line_args)

        self.assertEqual(temp_dir_path, args.x_working_dir)
        self.assertEqual(False, args.is_source_fastq_external)
        self.assertEqual(os.path.realpath(os.path.join(temp_dir_path, source_fastq_dir)),
                         args.inputs_dir)


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
                             inputs_dir=real('INPUTS'),
                             is_source_fastq_external=True)
            watermelon_init._make_top_level_dirs(args)

            def is_dir(o):
                return os.path.isdir(os.path.join(temp_dir_path, o))
            actual_dirs = sorted([o for o in os.listdir(temp_dir_path) if is_dir(o)])
            self.assertEqual(['ANALYSIS_11_10_A', 'INPUTS'],
                             actual_dirs)

    def test_make_top_level_dirs_bypassInputsWhenSourceNotExternal(self):
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            real = functools.partial(os.path.join, temp_dir_path)
            args = Namespace(analysis_dir=real('ANALYSIS_11_10_A'),
                             inputs_dir=real('INPUTS'),
                             is_source_fastq_external=False)
            watermelon_init._make_top_level_dirs(args)

            def is_dir(o):
                return os.path.isdir(os.path.join(temp_dir_path, o))
            actual_dirs = sorted([o for o in os.listdir(temp_dir_path) if is_dir(o)])
            self.assertEqual(['ANALYSIS_11_10_A'],
                             actual_dirs)


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
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            orig_wd = os.getcwd()
            try:
                os.chdir(temp_dir_path)
                watermelon_run_dir = j('00-source_runs', '')
                source_run_basedir = j('DNASeqCore', '')

                _create_files(source_run_basedir, source_files)
                source_run_dirs = [j(source_run_basedir, 'Run_1', 'tuttle', 'project'),
                                   j(source_run_basedir, 'Run_2', 'tuttle', 'project')]

                watermelon_init._link_run_dirs(source_run_dirs, watermelon_run_dir)

                actual_files = list(_listfiles(j(temp_dir_path, watermelon_run_dir)))
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
        expected_files = [j(temp_dir_path, watermelon_run_dir, f) for f in expected_files]
        self.assertEqual(set(expected_files), set(actual_files))
        self.assertEqual(len(actual_files), actual_link_count)

    # def test_populate_inputs_source_dir(self):
    #     with TempDirectory() as temp_dir:
    #         temp_dir_path = temp_dir.path
    #         source_fastq_dir = os.path.join(temp_dir_path, 'source_fastq_dir')
    #         os.mkdir(source_fastq_dir)
    #         os.mkdir(os.path.join(source_fastq_dir, 'Sample_A'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_A', 'SA_1.fastq.gz'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_A', 'SA_2.fastq.gz'))
    #         os.mkdir(os.path.join(source_fastq_dir, 'Sample_B'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_B', 'SB_1.fastq.gz'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_B', 'SB_2.fastq.gz'))
    #         os.mkdir(os.path.join(source_fastq_dir, 'Sample_C'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_C', 'SC_1.fastq.gz'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_C', 'SC_2.fastq.gz'))
    #
    #         os.mkdir(os.path.join(temp_dir_path, 'working_dir'))
    #         inputs_source_dir = os.path.join(temp_dir_path, 'working_dir', 'inputs', '00-source')
    #         os.mkdir(os.path.join(temp_dir_path, 'working_dir', 'inputs'))
    #         os.mkdir(inputs_source_dir)
    #         samples = {'Sample_A' : os.path.join(source_fastq_dir, 'Sample_A'),
    #                    'Sample_B' : os.path.join(source_fastq_dir, 'Sample_B'),
    #                    'Sample_C' : os.path.join(source_fastq_dir, 'Sample_C'),}
    #         run_samples = {'a_b_c': samples}
    #
    #         watermelon_init._populate_inputs_source_dir(inputs_source_dir, run_samples)
    #
    #         def assertINodesMatch(dirA, dirB, filename):
    #             self.assertEqual(os.stat(os.path.join(dirA, filename)).st_ino,
    #                              os.stat(os.path.join(dirB, filename)).st_ino)
    #         def is_link(o):
    #             return os.path.islink(os.path.join(temp_dir_path, inputs_source_dir, o))
    #
    #
    #         input_run_dirs = os.listdir(inputs_source_dir)
    #         self.assertEqual(1, len(input_run_dirs))
    #         self.assertEqual('a_b_c', input_run_dirs[0])
    #         inputs_source_run_dir = os.path.join(inputs_source_dir, input_run_dirs[0])
    #         actual_dirs = sorted([o for o in os.listdir(inputs_source_run_dir) if not is_link(o)])
    #         self.assertEqual(['Sample_A', 'Sample_B', 'Sample_C'], actual_dirs)
    #         self.assertEqual(['SA_1.fastq.gz', 'SA_2.fastq.gz'],
    #                          os.listdir(os.path.join(inputs_source_run_dir, 'Sample_A')))
    #         assertINodesMatch(source_fastq_dir, inputs_source_run_dir, 'Sample_A/SA_1.fastq.gz')
    #         self.assertEqual(['SB_1.fastq.gz', 'SB_2.fastq.gz'],
    #                          sorted(os.listdir(os.path.join(inputs_source_run_dir, 'Sample_B'))))
    #         assertINodesMatch(source_fastq_dir, inputs_source_run_dir, 'Sample_B/SB_1.fastq.gz')
    #         self.assertEqual(['SC_1.fastq.gz', 'SC_2.fastq.gz'],
    #                          sorted(os.listdir(os.path.join(inputs_source_run_dir, 'Sample_C'))))
    #         assertINodesMatch(source_fastq_dir, inputs_source_run_dir, 'Sample_C/SC_1.fastq.gz')


    # def Xtest_populate_inputs_merged_dir(self):
    #     with TempDirectory() as temp_dir:
    #         temp_dir_path = temp_dir.path
    #         source_fastq_dir = os.path.join(temp_dir_path, 'source_fastq_dir')
    #         os.mkdir(source_fastq_dir)
    #         os.mkdir(os.path.join(source_fastq_dir, 'Sample_A'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_A','SA_1.fastq.gz'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_A','SA_2.fastq.gz'))
    #         os.mkdir(os.path.join(source_fastq_dir, 'Sample_B'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_B','SB_1.fastq.gz'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_B','SB_2.fastq.gz'))
    #         os.mkdir(os.path.join(source_fastq_dir, 'Sample_C'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_C','SC_1.fastq.gz'))
    #         touch(os.path.join(source_fastq_dir, 'Sample_C','SC_2.fastq.gz'))
    #
    #         os.mkdir(os.path.join(temp_dir_path, 'working_dir'))
    #         inputs_merged_dir = os.path.join(temp_dir_path, 'working_dir', 'inputs', '00-merged')
    #         os.mkdir(os.path.join(temp_dir_path, 'working_dir', 'inputs'))
    #         os.mkdir(inputs_source_dir)
    #         samples = {'Sample_A' : os.path.join(source_fastq_dir, 'Sample_A'),
    #                    'Sample_B' : os.path.join(source_fastq_dir, 'Sample_B'),
    #                    'Sample_C' : os.path.join(source_fastq_dir, 'Sample_C'),}
    #
    #         watermelon_init._populate_inputs_source_dir(inputs_merged_dir, samples)
    #
    #         def assertINodesMatch(dirA, dirB, filename):
    #             self.assertEqual(os.stat(os.path.join(dirA, filename)).st_ino,
    #                              os.stat(os.path.join(dirB, filename)).st_ino)
    #         def is_link(o):
    #             return os.path.islink(os.path.join(temp_dir_path, inputs_source_dir, o))
    #         actual_dirs = sorted([o for o in os.listdir(inputs_source_dir) if not is_link(o)])
    #         self.assertEqual(['Sample_A', 'Sample_B', 'Sample_C'], actual_dirs)
    #         self.assertEqual(['SA_1.fastq.gz', 'SA_2.fastq.gz'],
    #                          os.listdir(os.path.join(inputs_source_dir, 'Sample_A')))
    #         assertINodesMatch(source_fastq_dir, inputs_source_dir, 'Sample_A/SA_1.fastq.gz')
    #         self.assertEqual(['SB_1.fastq.gz', 'SB_2.fastq.gz'],
    #                          os.listdir(os.path.join(inputs_source_dir, 'Sample_B')))
    #         assertINodesMatch(source_fastq_dir, inputs_source_dir, 'Sample_B/SB_1.fastq.gz')
    #         self.assertEqual(['SC_1.fastq.gz', 'SC_2.fastq.gz'],
    #                          os.listdir(os.path.join(inputs_source_dir, 'Sample_C')))
    #         assertINodesMatch(source_fastq_dir, inputs_source_dir, 'Sample_C/SC_1.fastq.gz')



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
    gtf: /ccmb/BioinfCore/noRibo_human.gtf
    bowtie2_index: /ccmb/BioinfCore/iGenomes/hg19/Sequence/Bowtie2Index
    entrez_gene_info: /ccmb/BioinfCore/entrez_gene_info/2016_09_02/gene_info
''')
        input_dir = '/my/input/dir'
        samples = {'sA':'fooA', 'sB': 'fooB', 'sC': 'fooC',
                   'sD':'fooA', 'sE': 'fooB', 'sF': 'fooC',
                  }
        actual_config = watermelon_init._make_config_dict(template_config,
                                                          genome_references,
                                                          input_dir,
                                                          samples)
        expected_keys = ['input_dir',
                         'foo1', 'foo2',
                         'genome', 'references',
                         'main_factors', 'phenotypes', 'samples', 'comparisons']
        self.assertEquals(sorted(expected_keys), sorted(actual_config.keys()))
        self.assertEqual(input_dir, actual_config['input_dir'])
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
                           'templateC' : 'TEMPLATE_C',
                           'phenotypes' :   'PHENOLABEL1 ^ PHENOLABEL2',
                           'main_factors' : 'yes ^ no'}
            watermelon_init._write_config_file(config_filename, config_dict)

            with open(config_filename, 'r') as config_file:
                config_lines = [line.strip('\n') for line in config_file.readlines()]
        config_prelude = watermelon_init._CONFIG_PRELUDE.split('\n')
        self.assertEqual(12 + len(config_prelude), len(config_lines))
        line_iter = iter(config_lines)
        for _ in config_prelude:
            next(line_iter)
        self.assertEqual('input_dir: INPUT_DIR', next(line_iter))
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
    def ok(self):
        #pylint: disable=locally-disabled,redundant-unittest-assert
        self.assertTrue(True)

    def test_validate_source_fastq_dir_raisesIfInvalid(self):
        validator = watermelon_init._CommandValidator()
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            source_fastq_dir = os.path.join(temp_dir_path, 'source')
            args = Namespace(source_fastq_dir=source_fastq_dir)
            msg = (r'Specified source_fastq_dir \[.*\] is not a dir or '
                   r'cannot be read. ')
            self.assertRaisesRegexp(watermelon_init._UsageError,
                                    msg,
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
            self.assertRaisesRegexp(watermelon_init._UsageError,
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
            args = Namespace(analysis_dir=analysis_dir,
                             inputs_dir=inputs_dir,
                             is_source_fastq_external=True)
            self.assertRaisesRegexp(watermelon_init._UsageError,
                                    r'inputs_dir \[.*\] exist\(s\)',
                                    validator._validate_overwrite_check,
                                    args)


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
            orig_wd = os.getcwd()
            try:
                os.chdir(temp_dir_path)
                script_name = os.path.join(_BIN_DIR, 'watermelon_init')

                redirect_output = ' 2>&1 '
                command = '{} {}'.format(script_name, redirect_output)
                exit_code, command_output = self.execute(command)
            finally:
                os.chdir(orig_wd)
        self.assertEqual(2, exit_code, command)
        self.assertRegexpMatches(command_output, 'usage', command)
