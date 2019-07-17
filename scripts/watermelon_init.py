#!/usr/bin/env python
'''Creates config file, input directories, and validates a sample description file for a watermelon rnaseq job.

Specifically:
1) watermelon_init accepts a set of source fastq dirs and a sample description file.
   * Each source fastq dir (a "run") contains a set of sample_dirs
   * watermelon_init verifies that the sample names provided in the sample
     description file match with those gathered from the source fastq dirs.
   * Each sample_dir contains one or more fastq (.fastq[.gz]) files.
   * watermelon_init creates local input dirs for the runs and also merges
     fastqs from all runs to create a dir of samples. Each sample dir
     is a flat dir of all fastqs for that sample.
   * Absence of fastq files for a run or an individaul sample raises an error.

2) watermelon-init creates an analysis directory which contains a template
   watermelon config file. The template config file contains directory
   information, reference lists, run parameters, etc.,
   These must be reviewed and edited to:
   a) adjust genome and references to match the experiment
   b) adjust trimming, alignment, and fastq_screen options
   c) specify diffex parameters, and add  DESeq2 calls and contrasts
      based on the example pheno_gender stanza (can add as many similar stanzas
      as required by the experiment)

3) watermelon-init creates a readme file that lists basic info about when/how it was run,
   what it did, and what the user has to do to prepare the template config
'''
#pylint: disable=locally-disabled,no-member
from __future__ import print_function, absolute_import, division
import argparse
from collections import defaultdict, OrderedDict
import datetime
import errno
import functools
import itertools
import os
import shutil
import sys
import time
import traceback

import pandas as pd
import yaml
import pdb

import scripts.watermelon_config as watermelon_config
from scripts.watermelon_config import CONFIG_KEYS

DESCRIPTION = \
'''Creates template config file and directories for a watermelon rnaseq job.'''

_PATH_SEP = '^'


class _InputValidationError(Exception):
    '''Raised for problematic input data.'''
    def __init__(self, msg, *args):
        super(_InputValidationError, self).__init__(msg, *args)


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

    def validate_inputs(self, input_summary):
        for validation in [self._validate_runs_have_fastq_files,
                           self._validate_samples_have_fastq_files]:
            validation(input_summary)

    @staticmethod
    def _validate_samples_have_fastq_files(input_summary):
        if input_summary.samples_missing_fastq_files:
            print(_build_input_summary_text(input_summary),
                  file=sys.stderr)
            msg = ('Some samples missing fastq files: [{}]. '
                   'Review details above against inputs and try again.'
                  ).format(', '.join(input_summary.samples_missing_fastq_files))
            raise _InputValidationError(msg)

    @staticmethod
    def _validate_runs_have_fastq_files(input_summary):
        if input_summary.runs_missing_fastq_files:
            print(_build_input_summary_text(input_summary),
                  file=sys.stderr)
            msg = ('Some runs missing fastq files: [{}]. '
                   'Review details above against inputs and try again.'
                  ).format(', '.join(input_summary.runs_missing_fastq_files))
            raise _InputValidationError(msg)


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
        if os.path.exists(args.input_dir):
            existing_files['dirs: input'] = args.input_dir
        if existing_files:
            file_types = [file_type for file_type in sorted(existing_files.keys())]
            file_names = [existing_files[key] for key in file_types]
            msg_format = ('{file_types} [{file_names}] exists (will not overwrite). '
                          'Remove these dir(s)/files(s) and try again.')
            msg = msg_format.format(file_types=",".join(file_types),
                                    file_names=",".join(file_names))
            raise _UsageError(msg)


class _Linker(object):
    '''Attempts to hard link file falling back to symlink as necessary.
    Only copies files.'''
    def __init__(self):
        self.links = defaultdict(int)

    def _link(self, source, dest):
        try:
            os.link(source, dest)
            self.links['hardlinked'] += 1
        except OSError as e:
            if e.errno != errno.EXDEV:
                raise e
            os.symlink(source, dest)
            self.links['symlinked'] += 1

    def copytree(self, source, dest):
        shutil.copytree(source, dest, copy_function=self._link)

    @property
    def results(self):
        x = ', '.join(['{} {}'.format(count, type) for type, count in sorted(self.links.items())])
        return 'linked {} source files: {}'.format(sum(self.links.values()), x)



def _timestamp():
    return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')

def _umich_email():
    return os.getlogin() + '@umich.edu'

README_FILENAME = 'watermelon.README'

GENOME_BUILD_OPTIONS = ('GRCh37', 'GRCh38', 'hg19', 'hg38',
                        'mm10',
                        'rn5', 'rn6',
                        'ce10', 'ce11', 'WBS235',
                        'GRCz10',
                        'ecoMG1655', 'ecoUTI89',
                        'dm6',
                        )

_SCRIPTS_DIR = os.path.realpath(os.path.dirname(__file__))
_WATERMELON_ROOT = os.path.dirname(_SCRIPTS_DIR)
_CONFIG_DIR = os.path.join(_WATERMELON_ROOT, 'config')
_DEFAULT_TEMPLATE_CONFIG = os.path.join(_CONFIG_DIR, 'template_config.yaml')
_DEFAULT_GENOME_REFERENCES = os.path.join(_CONFIG_DIR, 'genome_references.yaml')

_FASTQ_GLOBS = ['*.fastq', '*.fastq.gz']

_TODAY = datetime.date.today()
_DEFAULT_JOB_SUFFIX = '_{:02d}_{:02d}'.format(_TODAY.month, _TODAY.day)
_CONFIG_PRELUDE = '''# config created {timestamp}
# To do:
# ------
# 1) Review genome and references
# 2) Review alignment, trimming, and fastq_screen options
# 3) Modify the diffex options - use the pheno_gender stanza as an example, and
#    set up the DESeq2 calls and comparisons that match the analysis. For example,
#    the contrasts should contain phenotype labels and values which correspond
#    to the column names and values in the sample sheet, respectively.
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

import collections

def _dict_merge(dct, merge_dct):
    """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    # https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            _dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

def _link_run_dirs(args, linker):
    source_run_dirs = args.source_fastq_dirs
    watermelon_runs_dir = os.path.join(args.tmp_input_dir, args.input_runs_dir)
    _mkdir(watermelon_runs_dir)
    for run_dir in source_run_dirs:
        abs_run_dir = os.path.abspath(run_dir)
        mangled_dirname = abs_run_dir.replace(os.path.sep, _PATH_SEP)
        mangled_dirname.rstrip(os.path.sep)
        watermelon_run_dir = os.path.join(watermelon_runs_dir, mangled_dirname)
        linker.copytree(run_dir, watermelon_run_dir)

def _merge_sample_dirs(args):
    j = os.path.join
    watermelon_runs_dir = os.path.abspath(j(args.tmp_input_dir,
                                            args.input_runs_dir))
    watermelon_samples_dir = os.path.abspath(j(args.tmp_input_dir,
                                               args.input_samples_dir))
    _mkdir(watermelon_samples_dir)
    for run_dir in os.listdir(watermelon_runs_dir):
        for sample_dir in os.listdir(j(watermelon_runs_dir, run_dir)):
            new_sample_dir = j(watermelon_samples_dir, sample_dir)
            _mkdir(new_sample_dir)
            files = os.listdir(j(watermelon_runs_dir, run_dir, sample_dir))
            if files:
                for sample_file in files:
                    source = j(watermelon_runs_dir, run_dir, sample_dir, sample_file)
                    link = j(new_sample_dir,
                             run_dir + _PATH_SEP + sample_file)
                    os.link(source, link)
            else:
                from pathlib import Path
                Path(j(new_sample_dir, run_dir + _PATH_SEP + 'empty_dir')).touch()

def _is_fastq_file(file_name):
    return file_name.endswith('.fastq') or file_name.endswith('.fastq.gz')

def _build_input_summary(args):
    input_samples_dir = os.path.join(args.tmp_input_dir, args.input_samples_dir)
    rows = []
    for root, dirs, files in os.walk(input_samples_dir):
        for name in files:
            file_name = os.path.basename(name)
            sample_name = os.path.basename(root)
            run_name = os.path.dirname(file_name.replace(_PATH_SEP, os.path.sep))
            rows.append((run_name, sample_name, file_name))

    df = pd.DataFrame(data=rows, columns=['run', 'sample', 'file'])
    df['run_index'] = df['run'].rank(method='dense').apply(lambda x: chr(64+int(x)))
    df['is_fastq_file'] = df['file'].apply(_is_fastq_file)

    runs_df = df.drop_duplicates(subset=['run_index', 'run'])[['run_index', 'run']].set_index('run_index').sort_index(axis=0)
    sample_run_counts_df = df.pivot_table(values='is_fastq_file',
                                          index=['sample'],
                                          columns='run_index',
                                          aggfunc='sum')
    sample_run_counts_df = sample_run_counts_df.fillna(0).astype(int, errors='ignore')

    total_sample_count = len(sample_run_counts_df)
    total_file_count = sum(sum(sample_run_counts_df.values))
    samples = sorted(sample_run_counts_df.index)
    sample_run_counts_df.loc['run_total'] = sample_run_counts_df.sum()
    sample_run_counts_df['sample_total'] = sample_run_counts_df.sum(axis=1)
    samples_missing_fastq_files = sorted(sample_run_counts_df[sample_run_counts_df['sample_total'] == 0].index)

    run_counts_df = sample_run_counts_df[sample_run_counts_df.index == 'run_total'].transpose()
    runs_missing_fastq_files = sorted(run_counts_df[run_counts_df['run_total'] == 0].index)

    input_summary = argparse.Namespace(runs_df=runs_df,
                                       sample_run_counts_df=sample_run_counts_df,
                                       total_sample_count=total_sample_count,
                                       total_file_count=total_file_count,
                                       total_run_count=len(runs_df),
                                       samples=samples,
                                       samples_missing_fastq_files=samples_missing_fastq_files,
                                       runs_missing_fastq_files=runs_missing_fastq_files)
    return input_summary

'''
Verify that the samplesheet exists, and that the sample names
in the samplesheet correspond to those in the input directories
'''
def _validate_sample_sheet(samples, sheet_file):
    #pylint: disable=locally-disabled,invalid-name

    if not os.path.isfile(sheet_file):
        msg = ("The specified sample sheet file {} cannot be found.".format(sheet_file))
        raise _InputValidationError(msg)

    samplesheet = pd.read_csv(sheet_file).set_index("sample", drop=True)
    sheet_samples = set(samplesheet.index)
    input_samples = set(samples)

    if not sheet_samples.issubset(input_samples):
        msg = ("The following sample names cannot be found in the input directories: "
        + str(sheet_samples.difference(input_samples)))
        raise _InputValidationError(msg)

def _make_config_dict(template_config, genome_references, args, samples):
    config = dict(template_config)

    dirs = config.get(CONFIG_KEYS.dirs, {})
    dirs[CONFIG_KEYS.dirs_input] = os.path.join(args.input_dir,
                                                args.input_samples_dir)
    config[CONFIG_KEYS.dirs] = dirs

    if 'alignment_options' in config:
        config['alignment_options']['rsem_ref_prefix'] = args.genome_build

    _dict_merge(config, genome_references)

    return config

def _build_input_summary_text(input_summary):
    pd.set_option('display.max_colwidth', 300)
    text = \
'''Input summary:
--------------
Found {total_file_count} files for {total_sample_count} samples across {total_run_count} run(s).

Run dir(s):
-----------
{run_details}

Sample x run fastq file counts:
-------------------------------
{sample_run_counts}
'''.format(total_file_count=input_summary.total_file_count,
           total_sample_count=input_summary.total_sample_count,
           total_run_count=input_summary.total_run_count,
           run_details=input_summary.runs_df.to_string(justify='left'),
           sample_run_counts=input_summary.sample_run_counts_df.to_string(justify='left'))
    return text

def _build_postlude(args, linker_results, input_summary_text):
    input_dir_relative = os.path.relpath(args.input_dir, args.x_working_dir)
    postlude = \
'''
watermelon_init.README
======================

{linker_results}
{input_summary_text}

Created files and dirs
----------------------
{working_dir}/
    {config_basename}
    watermelon
    {input_dir_relative}/
        {input_runs_dir}/
        {input_samples_dir}/

You need to review the sample sheet and config file:
samplesheet:
{samplesheet_relative}
-------------------------------
1) Review sample phenotype labels and values for each sample.
   Phenotype labels must be distinct.
   Phenotype labels and values must be valid R column names, so each label/value must
   be a letter followed alphanumerics or [-_.]. (So "A24-5" is ok, but "1hr+" is not.
   Also the literals "T", "F", and "NAN" cannot be used as phenotype labels/values.

config:
{config_relative}
-------------------------------
1) Review genome and references
2) Review alignment, trimming, and fastq_screen options
3) Modify the diffex options - use the pheno_gender stanza as an example, and
   set up the DESeq2 calls and comparisons that match the analysis. For example,
   the contrasts should contain phenotype labels and values which correspond
   to the column names and values in the sample sheet.

When the config file looks good:
--------------------------------
# start a screen session:
$ screen -S watermelon{job_suffix}
# to validate the config and check the execution plan:
$ watermelon --dry-run -c {config_basename}
# to run:
$ watermelon -c {config_basename}
'''.format(linker_results=linker_results,
           input_summary_text=input_summary_text,
           working_dir=args.x_working_dir,
           input_dir_relative=input_dir_relative,
           input_runs_dir=args.input_runs_dir,
           input_samples_dir=args.input_samples_dir,
           config_file=args.config_file,
           config_relative=os.path.relpath(args.config_file, args.x_working_dir),
           samplesheet_relative=os.path.relpath(args.sample_sheet, args.x_working_dir),
           config_basename=os.path.basename(args.config_file),
           job_suffix=args.job_suffix,)
    return postlude

def _write_config_file(config_filename, config_dict):
    tmp_config_dict = dict(config_dict)
    ordered_config_dict = OrderedDict()
    named_keys = [CONFIG_KEYS.dirs,
                  CONFIG_KEYS.main_factors,
                  CONFIG_KEYS.phenotypes,
                  CONFIG_KEYS.samples,
                  CONFIG_KEYS.comparisons,
                  CONFIG_KEYS.genome,
                  CONFIG_KEYS.references]
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
        '--sample_sheet',
        type=str,
        required=True,
        help=('A CSV file containing sample names and phenotype/characteristic '
        'information which correspond to these samples. Watermelon-init will verify that '
        'sample names in this file match with those found in the input directories. '),
        )
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
    args.config_file = 'config{}.yaml'.format(args.job_suffix)
    args.tmp_input_dir = realpath('.inputs')
    args.input_runs_dir = '00-source_runs'
    args.input_samples_dir = '01-source_samples'
    args.input_dir = realpath('inputs')
    return args

def main(sys_argv):
    '''See DESCRIPTION'''
    try:
        args = _parse_command_line_args(sys_argv)
        _CommandValidator().validate_args(args)

        print("Copying Watermelon source to " + os.getcwd())
        shutil.copytree(_WATERMELON_ROOT, os.path.join(os.getcwd(), "Watermelon"))
        os.symlink("Watermelon/bin/watermelon", "watermelon")

        if os.path.isdir(args.tmp_input_dir):
            shutil.rmtree(args.tmp_input_dir)
        linker = _Linker()
        _link_run_dirs(args, linker)
        _merge_sample_dirs(args)
        input_summary = _build_input_summary(args)

        _CommandValidator().validate_inputs(input_summary)
        os.rename(args.tmp_input_dir, args.input_dir)

        with open(args.x_template_config, 'r') as template_config_file:
            template_config = yaml.load(template_config_file, yaml.FullLoader)
        with open(args.x_genome_references, 'r') as genome_references_file:
            genome_references = yaml.load(genome_references_file, yaml.FullLoader)[args.genome_build]
        config_dict = _make_config_dict(template_config,
                                        genome_references,
                                        args,
                                        input_summary.samples)

        _validate_sample_sheet(input_summary.samples, args.sample_sheet)
        config_dict['sample_description_file'] = os.path.abspath(args.sample_sheet)
        _write_config_file(args.config_file, config_dict)

        postlude = _build_postlude(args,
                                   linker.results,
                                   _build_input_summary_text(input_summary))
        print(postlude)
        with open(README_FILENAME, 'w') as readme:
            print(postlude, file=readme)
    except _UsageError as usage_error:
        message = "watermelon-init usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'watermelon-init --help'.", file=sys.stderr)
        sys.exit(1)
    except _InputValidationError as input_validation_error:
        print(input_validation_error, file=sys.stderr)
        sys.exit(1)
    except Exception: #pylint: disable=locally-disabled,broad-except
        print("An unexpected error occurred", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        exit(1)

_setup_yaml()
if __name__ == '__main__':
    main(sys.argv[1:])
