#!/usr/bin/env python
'''Creates template config file and directories for a watermelon rnaseq job.

Specifically this does three things:
1) watermelon-init accepts a source fastq dir which would typically contain a set of
   sample_dirs each of which would contain fastq files for that sample. To encapsulate
   the project data flow, while enabling data provenance, init creates a local inputs dir
   and creates symlinks to the original source dir.

   Note that this step is skipped if the source fastq is inside the working dir (i.e.
   already local).

2) watermelon-init creates an analysis directory which contains a template watermelon
   config file. The template config file is created by combining the specified genome
   with the sample names in (from the inputs dir). The template config must be reviewed
   and edited to:
   a) adjust run params
   b) add sample phenotypes/aliases
   c) specify sample comparisons

3) watermelon-init creates a readme file that lists basic info about when/how it was run,
   what it did, and what the user has to do to prepare the template config
'''
#pylint: disable=locally-disabled,no-member
from __future__ import print_function, absolute_import, division
import argparse
from collections import defaultdict, OrderedDict
import datetime
import functools
import glob
import itertools
import os
import shutil
import sys
import time
import traceback

import pandas as pd
import yaml

import scripts.watermelon_config as watermelon_config

DESCRIPTION = \
'''Creates template config file and directories for a watermelon rnaseq job.'''

_PATH_SEP = '^'


class _UsageError(Exception):
    '''Raised for malformed command or invalid arguments.'''
    def __init__(self, msg, *args):
        super(_UsageError, self).__init__(msg, *args)


class _CommandValidator(object):
    #pylint: disable=locally-disabled,too-few-public-methods
    def __init__(self):
        pass

    def validate_args(self, args):
        #pylint: disable=locally-disabled, missing-docstring
        for validation in [self._validate_source_fastq_dirs,
                           self._validate_overwrite_check]:
            validation(args)

    @staticmethod
    def _validate_source_fastq_dirs(args):
        bad_dirs = [x for x in args.source_fastq_dirs if not os.path.isdir(x)]
        if bad_dirs:
            msg = ('Specified source_fastq_dir(s) [{}] not a dir or cannot be '
                   'read. Review inputs and try again.'
                  ).format(','.join(sorted(bad_dirs)))
            raise _UsageError(msg)

    @staticmethod
    def _validate_overwrite_check(args):
        existing_files = {}
        if os.path.exists(args.analysis_dir):
            existing_files['analysis_dir'] = args.analysis_dir
        if os.path.exists(args.input_runs_dir):
            existing_files['input_runs_dir'] = args.input_runs_dir
        if os.path.exists(args.input_samples_dir):
            existing_files['input_samples_dir'] = args.input_samples_dir
        if existing_files:
            file_types = [file_type for file_type in sorted(existing_files.keys())]
            file_names = [existing_files[key] for key in file_types]
            msg_format = ('{file_types} [{file_names}] exist(s) (will not overwrite). '
                          'Remove these dir(s)/files(s) and try again.')
            msg = msg_format.format(file_types=",".join(file_types),
                                    file_names=",".join(file_names))
            raise _UsageError(msg)

def _timestamp():
    return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')

README_FILENAME = 'watermelon.README'

GENOME_BUILD_OPTIONS = ('hg19', 'mm10', 'rn5')

_SCRIPTS_DIR = os.path.realpath(os.path.dirname(__file__))
_WATERMELON_ROOT = os.path.dirname(_SCRIPTS_DIR)
_CONFIG_DIR = os.path.join(_WATERMELON_ROOT, 'config')
_DEFAULT_TEMPLATE_CONFIG = os.path.join(_CONFIG_DIR, 'template_config.yaml')
_DEFAULT_GENOME_REFERENCES = os.path.join(_CONFIG_DIR, 'genome_references.yaml')

_FASTQ_GLOBS = ['*.fastq', '*.fastq.gz']

_DEFAULT_MAIN_FACTORS = 'yes    | yes      | no'.\
    replace('|', watermelon_config.DEFAULT_PHENOTYPE_DELIM)
_DEFAULT_PHENOTYPE_LABELS = 'gender | genotype | gender.genotype'.\
    replace('|', watermelon_config.DEFAULT_PHENOTYPE_DELIM)
_DEFAULT_GENDER_VALUES = ['female', 'male']
_DEFAULT_GENOTYPE_VALUES = ['MutA', 'MutB', 'WT']
_DEFAULT_COMPARISONS = {'gender'   : ['male_v_female'],
                        'genotype' : ['MutA_v_WT', 'MutB_v_WT']
                       }

_TODAY = datetime.date.today()
_DEFAULT_JOB_SUFFIX = '_{:02d}_{:02d}'.format(_TODAY.month, _TODAY.day)
_CONFIG_PRELUDE = '''# config created {timestamp}
# To do:
# ------
# 1) Review/adjust samples names
#     Note: names in config, must match the sample dir names in the input dir
# 2) Add a sample group for each sample
# 3) Add comparisons
# 4) Review genome and references
# 5) Review alignment options
# 6) Review trimming options
'''.format(timestamp=_timestamp())

def _setup_yaml():
    """ http://stackoverflow.com/a/8661021 """
    represent_dict_order = lambda self, data: \
        self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)

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

def _initialize_samples(watermelon_sample_dir):
    def _count_fastq(directory):
        count = 0
        for fastq_glob in _FASTQ_GLOBS:
            count += len(glob.glob(os.path.join(directory, fastq_glob)))
        return count

    samples = {}
    for file_name in os.listdir(watermelon_sample_dir):
        full_path = os.path.realpath(os.path.join(watermelon_sample_dir, file_name))
        if os.path.isdir(full_path) and _count_fastq(full_path):
            samples[file_name] = full_path
    return samples

def _link_run_dirs(source_run_dirs, watermelon_runs_dir):
    _mkdir(watermelon_runs_dir)
    for run_dir in source_run_dirs:
        abs_run_dir = os.path.abspath(run_dir)
        mangled_dirname = abs_run_dir.replace(os.path.sep, _PATH_SEP)
        mangled_dirname.rstrip(os.path.sep)
        watermelon_run_dir = os.path.join(watermelon_runs_dir, mangled_dirname)
        shutil.copytree(run_dir, watermelon_run_dir, copy_function=os.link)

def _merge_sample_dirs(watermelon_runs_dir, watermelon_samples_dir):
    j = os.path.join
    watermelon_runs_dir = os.path.abspath(watermelon_runs_dir)
    watermelon_samples_dir = os.path.abspath(watermelon_samples_dir)
    sample_source_link = defaultdict(dict)
    for run_dir in os.listdir(watermelon_runs_dir):
        for sample_dir in os.listdir(j(watermelon_runs_dir, run_dir)):
            for sample_file in os.listdir(j(watermelon_runs_dir, run_dir, sample_dir)):
                source = j(watermelon_runs_dir, run_dir, sample_dir, sample_file)
                link = j(watermelon_samples_dir,
                         sample_dir,
                         run_dir + _PATH_SEP + sample_file)
                new_sample_dir = j(watermelon_samples_dir, sample_dir)
                sample_source_link[new_sample_dir][source] = link

    _mkdir(watermelon_samples_dir)
    for sample_dir, source_links in sample_source_link.items():
        _mkdir(sample_dir)
        for source, link in source_links.items():
            os.link(source, link)

def _get_samples_and_file_count(watermelon_samples_dir):
    return 'no', 0

def _build_phenotypes_samples_comparisons(samples):
    #pylint: disable=locally-disabled,invalid-name
    config = {}
    config[watermelon_config.CONFIG_KEYS.main_factors] = _DEFAULT_MAIN_FACTORS
    config[watermelon_config.CONFIG_KEYS.phenotypes] = _DEFAULT_PHENOTYPE_LABELS

    gender_values = itertools.cycle(_DEFAULT_GENDER_VALUES)
    genotype_values = itertools.cycle(_DEFAULT_GENOTYPE_VALUES)
    pheno_value_items = []
    for _ in samples:
        gender_value = next(gender_values)
        genotype_value = next(genotype_values)
        gender_genotype_value = '.'.join([gender_value, genotype_value])
        pheno_value_items.append((gender_value, genotype_value, gender_genotype_value))
    pheno_fmt = "{0:7}| {1:8} | {2}".replace('|', watermelon_config.DEFAULT_PHENOTYPE_DELIM)
    pheno_value_strings = [pheno_fmt.format(*values) for values in sorted(pheno_value_items)]
    samples_dict = dict([(s, v) for s, v in zip(sorted(samples), pheno_value_strings)])
    config[watermelon_config.CONFIG_KEYS.samples] = samples_dict

    config[watermelon_config.CONFIG_KEYS.comparisons] = _DEFAULT_COMPARISONS
    return config

def _make_config_dict(template_config, genome_references, input_dir, samples):
    config = dict(template_config)
    config[watermelon_config.CONFIG_KEYS.input_dir] = '{}'.format(input_dir)
    config.update(genome_references)
    config.update(_build_phenotypes_samples_comparisons(samples))
    return config

def _build_input_summary(input_samples_dir):
    pd.set_option('display.max_colwidth', 300)
    rows = []
    for root, dirs, files in os.walk(input_samples_dir):
        for name in files:
            file_name = os.path.basename(name)
            sample_name = os.path.basename(root)
            run_name = os.path.dirname(file_name.replace(_PATH_SEP, os.path.sep))
            rows.append((run_name, sample_name, file_name))

    df = pd.DataFrame(data=rows, columns=['run', 'sample', 'file'])
    df['run_index'] = df['run'].rank(method='dense').apply(lambda x: chr(64+int(x)))

    runs_df = df.drop_duplicates(subset=['run_index', 'run'])[['run_index', 'run']].set_index('run_index')
    sample_run_counts_df = df.pivot_table(values='file',index=['sample'], columns='run_index', aggfunc='count').fillna(0).astype(int, errors='ignore')

    total_sample_count = len(sample_run_counts_df)
    total_file_count = sum(sum(sample_run_counts_df.values))
    sample_run_counts_df.loc['run_total']= sample_run_counts_df.sum()
    sample_run_counts_df['sample_total']= sample_run_counts_df.sum(axis=1)

    input_summary = argparse.Namespace(runs_df = runs_df,
                                       sample_run_counts_df=sample_run_counts_df,
                                       total_sample_count=total_sample_count,
                                       total_file_count=total_file_count,
                                       total_run_count=len(runs_df))
    return input_summary

def _build_postlude(args, input_summary):
    input_dir_relative = os.path.relpath(args.input_dir, args.x_working_dir)
    postlude = \
'''
watermelon_init.README
======================

Input summary:
--------------
Found {total_file_count} files for {total_sample_count} samples across {total_run_count} run(s).

Source fastq dirs details (i.e. source sequencing run dirs):
{run_details}

Sample x run file counts:
{sample_run_counts}

Created files and dirs
----------------------
{working_dir}/
    {input_dir_relative}/
        {input_runs_dir}/
        {input_samples_dir}/
    {analysis_relative}/
        {config_basename}

You need to review config file: {config_relative}:
-------------------------------
1) Review/adjust samples names
    Note: if you change names in config, also change the sample dir names in the input dir
2) Review sample phenotype labels and values for each sample.
   Phenotype labels must be distinct.
   Phenotype labels and values must be valid R column names, so each label/value must
   be a letter followed alphanumerics or [-_.]. (So "A24-5" is ok, but "1hr+" is not.
   Also the literals "T", "F", and "NAN" cannot be used as phenotype labels/values.
3) Adjust the main_factors line to indicate whether a phenotype is main (yes) or derived (no).
4) Add comparisons
5) Review genome and references
6) Review alignment options
7) Review trimming options

When the config file looks good:
--------------------------------
$ cd {analysis_dir}
# start a screen session:
$ screen -S watermelon{job_suffix}
# to validate the config and check the execution plan:
$ watermelon --dry-run -c {config_basename}
# to run:
$ watermelon -c {config_basename}
'''.format(total_file_count=input_summary.total_file_count,
           total_sample_count=input_summary.total_sample_count,
           total_run_count=input_summary.total_run_count,
           run_details=input_summary.runs_df.to_string(justify='left'),
           sample_run_counts=input_summary.sample_run_counts_df.to_string(justify='left'),
           working_dir=args.x_working_dir,
           input_dir_relative=input_dir_relative,
           input_runs_dir=os.path.relpath(args.input_runs_dir, input_dir_relative),
           input_samples_dir=os.path.relpath(args.input_samples_dir, input_dir_relative),
           analysis_dir=args.analysis_dir,
           analysis_relative=os.path.relpath(args.analysis_dir, args.x_working_dir),
           config_file=args.config_file,
           config_relative=os.path.relpath(args.config_file, args.x_working_dir),
           config_basename=os.path.basename(args.config_file),
           job_suffix=args.job_suffix,)
    return postlude

def _write_config_file(config_filename, config_dict):
    tmp_config_dict = dict(config_dict)
    ordered_config_dict = OrderedDict()
    named_keys = [watermelon_config.CONFIG_KEYS.input_dir,
                  watermelon_config.CONFIG_KEYS.main_factors,
                  watermelon_config.CONFIG_KEYS.phenotypes,
                  watermelon_config.CONFIG_KEYS.samples,
                  watermelon_config.CONFIG_KEYS.comparisons,
                  watermelon_config.CONFIG_KEYS.genome,
                  watermelon_config.CONFIG_KEYS.references]
    for key in named_keys:
        if key in tmp_config_dict:
            ordered_config_dict[key] = tmp_config_dict.pop(key)
    for key in sorted(tmp_config_dict):
        ordered_config_dict[key] = tmp_config_dict[key]
    with open(config_filename, 'w') as config_file:
        print(_CONFIG_PRELUDE, file=config_file)
        yaml.dump(ordered_config_dict, config_file, default_flow_style=False, indent=4)

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--genome_build',
        type=str,
        choices=GENOME_BUILD_OPTIONS,
        help='Config template will based on this genome',
        required=True)
    parser.add_argument(
        '--job_suffix',
        type=str,
        default=_DEFAULT_JOB_SUFFIX,
        help=('={} This suffix will be appended to the names of analysis/deliverables '
              'dirs and configfile. '
              '(Useful if rerunning a job with revised input fastq or revised config; '
              'can help differentiate simultaneous watermelon jobs in '
              'top/ps.)').format(_DEFAULT_JOB_SUFFIX),
        required=False)
    parser.add_argument(
        'source_fastq_dirs',
        type=str,
        nargs='+',
        help=('One or more paths to run dirs. Each run dir should contain '
              'samples dirs which should in turn contain one or more fastq.gz '
              'files. The watermelon sample names will be derived from the '
              'sample directories.'),
        )
    parser.add_argument(
        '--x_working_dir',
        type=str,
        default=os.getcwd(),
        help=argparse.SUPPRESS)

    parser.add_argument(
        '--x_template_config',
        type=str,
        default=_DEFAULT_TEMPLATE_CONFIG,
        help=argparse.SUPPRESS)

    parser.add_argument(
        '--x_genome_references',
        type=str,
        default=_DEFAULT_GENOME_REFERENCES,
        help=argparse.SUPPRESS)

    args = parser.parse_args(sys_argv)

    realpath = functools.partial(os.path.join, args.x_working_dir)
    args.analysis_dir = realpath('analysis{}'.format(args.job_suffix))
    args.config_file = os.path.join(args.analysis_dir,
                                    'config{}.yaml'.format(args.job_suffix))
    args.input_dir = realpath('inputs')
    args.input_runs_dir = os.path.join(args.input_dir, '00-source_runs')
    args.input_samples_dir = os.path.join(args.input_dir, '01-source_samples')
    return args

def main(sys_argv):
    '''See DESCRIPTION'''
    try:
        args = _parse_command_line_args(sys_argv)
        _CommandValidator().validate_args(args)

        _link_run_dirs(args.source_fastq_dirs, args.input_runs_dir)
        _merge_sample_dirs(args.input_runs_dir, args.input_samples_dir)
        input_summary = _build_input_summary(args.input_samples_dir)
        #check for samples without fastq files
        samples  = _initialize_samples(args.input_samples_dir)
        _mkdir(args.analysis_dir)

        with open(args.x_template_config, 'r') as template_config_file:
            template_config = yaml.load(template_config_file)
        with open(args.x_genome_references, 'r') as genome_references_file:
            genome_references = yaml.load(genome_references_file)[args.genome_build]
        config_dict = _make_config_dict(template_config,
                                        genome_references,
                                        args.input_samples_dir,
                                        samples)
        _write_config_file(args.config_file, config_dict)

        postlude = _build_postlude(args, input_summary)
        print(postlude)
        with open(README_FILENAME, 'w') as readme:
            print(postlude, file=readme)
    except _UsageError as usage_error:
        message = "watermelon-init usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'watermelon-init --help'.", file=sys.stderr)
        sys.exit(1)
    except Exception: #pylint: disable=locally-disabled,broad-except
        print("An unexpected error occurred", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        exit(1)

_setup_yaml()
if __name__ == '__main__':
    main(sys.argv[1:])
