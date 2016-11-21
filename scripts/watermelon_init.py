#!/usr/bin/env python
from __future__ import print_function, absolute_import, division
import argparse
from collections import OrderedDict
import datetime
import functools
import glob
import os
import sys
import time
import traceback

import yaml

README_FILENAME = 'watermelon.README'

class _UsageError(Exception):
    '''Raised for malformed command or invalid arguments.'''
    def __init__(self, msg, *args):
        super(_UsageError, self).__init__(msg, *args)

class _CommandValidator(object):
    def __init__(self):
        pass

    def validate_args(self, args):
        for validation in [self._validate_source_fastq_dir,
                           self._validate_overwrite_check]:
            validation(args)

    def _validate_source_fastq_dir(self, args):
        if not os.path.isdir(args.source_fastq_dir):
            msg = ('Specified source_fastq_dir [{}] is not a dir or cannot be read. '
                   'Review inputs and try again.').format(args.source_fastq_dir)
            raise _UsageError(msg)

    def _validate_overwrite_check(self, args):
        existing_files = {}
        if os.path.exists(args.analysis_dir):
            existing_files['analysis_dir'] = args.analysis_dir
        if os.path.exists(args.inputs_dir) and args.is_source_fastq_external:
            existing_files['inputs_dir'] = args.inputs_dir
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

DESCRIPTION = \
'''Creates template config file and directories for a watermelon rnaseq job.'''
GENOME_BUILD_OPTIONS = ('hg19', 'mm10', 'rn5')

_SCRIPTS_DIR = os.path.realpath(os.path.dirname(__file__))
_WATERMELON_ROOT = os.path.dirname(_SCRIPTS_DIR)
_CONFIG_DIR = os.path.join(_WATERMELON_ROOT, 'config')
_DEFAULT_TEMPLATE_CONFIG = os.path.join(_CONFIG_DIR, 'template_config.yaml')
_DEFAULT_GENOME_REFERENCES = os.path.join(_CONFIG_DIR, 'genome_references.yaml')

_FASTQ_GLOBS = ['*.fastq', '*.fastq.gz']

_CONFIG_KEYS = argparse.Namespace(comparisons='comparisons',
                                  genome='genome',
                                  input_dir='input_dir',
                                  references='references',
                                  samples='samples')
_DEFAULT_PHENOTYPE = 'g1'
_DEFAULT_COMPARISON = {'g1_v_g2' : 'g1_v_g2'}

_TODAY = datetime.date.today()
_DEFAULT_JOB_SUFFIX = '_{:02d}_{:02d}'.format(_TODAY.month, _TODAY.day)
_CONFIG_PRELUDE=\
'''# config created {timestamp}
# To do:
# ------
# 1) Review/adjust samples names
#     Note: if you change names in config, also change the sample dir names in the input dir
# 2) Add a sample group for each sample
# 3) Add comparisons
# 4) Review genome and references
# 5) Review alignment options
# 6) Review trimming options
'''.format(timestamp=_timestamp())

def _setup_yaml():
  """ http://stackoverflow.com/a/8661021 """
  represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
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

def _initialize_samples(source_fastq_dir):
    def _count_fastq(dir):
        count = 0
        for fastq_glob in _FASTQ_GLOBS:
            count += len(glob.glob(os.path.join(dir, fastq_glob)))
        return count

    samples = {}
    file_count = 0
    for file_name in os.listdir(source_fastq_dir):
        full_path = os.path.realpath(os.path.join(source_fastq_dir, file_name))
        if os.path.isdir(full_path) and _count_fastq(full_path):
            samples[file_name] = full_path
            file_count += _count_fastq(full_path)
    return samples, file_count

def _is_source_fastq_external(source_fastq_dir, working_dir=os.getcwd()):
    a = os.path.realpath(source_fastq_dir)
    b = os.path.realpath(working_dir)
    common_prefix = os.path.commonprefix([a, b])
    return common_prefix != b

def _populate_inputs_dir(inputs_dir, samples):
    for sample_name, sample_dir_target in samples.items():
        os.symlink(sample_dir_target, os.path.join(inputs_dir, sample_name))

def _make_config_dict(template_config, genome_references, input_dir, samples):
    config = dict(template_config)
    config[_CONFIG_KEYS.input_dir] = '"{}"'.format(input_dir)
    config.update(genome_references)
    config[_CONFIG_KEYS.samples] = dict([(name, _DEFAULT_PHENOTYPE)for name in samples])
    config[_CONFIG_KEYS.comparisons] = _DEFAULT_COMPARISON
    return config

def _build_postlude(args, sample_count, file_count):
    postlude = \
'''
watermelon_init.README
======================

Created files and dirs
----------------------
    {inputs_dir}
        source fastq dir | sample count | fastq file count
        {source_fastq_dir} | {sample_count} samples | {file_count} files
    {analysis_dir}
    {config_file}

You need to review config file: {config_file}:
-------------------------------
1) Review/adjust samples names
    Note: if you change names in config, also change the sample dir names in the input dir
2) Add a sample group for each sample
3) Add comparisons
4) Review genome and references
5) Review alignment options
6) Review trimming options

When the config file looks good:
--------------------------------
$ cd {analysis_dir}
# start a screen session:
$ screen -S watermelon{job_suffix}
# to validate the config and check the execution plan:
$ watermelon --dry-run -c {config_basename}
# to run:
$ watermelon -c {config_basename}
'''.format(inputs_dir=args.inputs_dir,
           source_fastq_dir=args.source_fastq_dir,
           sample_count=sample_count,
           file_count=file_count,
           analysis_dir=args.analysis_dir,
           config_file=args.config_file,
           config_basename=os.path.basename(args.config_file),
           job_suffix=args.job_suffix,)
    return postlude

def _make_top_level_dirs(args):
    _mkdir(args.analysis_dir)
    if args.is_source_fastq_external:
        _mkdir(args.inputs_dir)

def _write_config_file(config_filename, config_dict):
    tmp_config_dict = dict(config_dict)
    ordered_config_dict = OrderedDict()
    named_keys = [_CONFIG_KEYS.input_dir,
                  _CONFIG_KEYS.samples,
                  _CONFIG_KEYS.comparisons,
                  _CONFIG_KEYS.genome,
                  _CONFIG_KEYS.references]
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
        'source_fastq_dir',
        type=str,
        help=('Path to dir which contains samples dirs (each sample dir '
              'contains one or more fastq.gz files). The sample dir names are used to '
              'initialize the config template. If the source_fastq_dir is outside the '
              'working dir, init will make local inputs dir containing symlinks to the '
              'original files.'),
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
    args.inputs_dir = realpath('inputs', '00-multiplexed_reads')
    args.is_source_fastq_external = _is_source_fastq_external(args.source_fastq_dir,
                                                              args.x_working_dir)

    return args

def main(sys_argv):
    try:
        args = _parse_command_line_args(sys_argv)
        _CommandValidator().validate_args(args)

        (samples, file_count) = _initialize_samples(args.source_fastq_dir)
        _make_top_level_dirs(args)
        if args.is_source_fastq_external:
            _populate_inputs_dir(args.inputs_dir, samples)

        with open(args.x_template_config, 'r') as template_config_file:
            template_config = yaml.load(template_config_file)
        with open(args.x_genome_references, 'r') as genome_references_file:
            genome_references = yaml.load(genome_references_file)[args.genome_build]
        config_dict = _make_config_dict(template_config,
                                        genome_references,
                                        args.inputs_dir,
                                        samples)
        _write_config_file(args.config_file, config_dict)

        postlude = _build_postlude(args, len(samples), file_count)
        print(postlude)
        with open(README_FILENAME, 'w') as readme:
            print(postlude, file=readme)
    except _UsageError as usage_error:
        message = "watermelon-init usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'watermelon-init --help'.", file=sys.stderr)
        sys.exit(1)
    except Exception: #pylint: disable=broad-except
        print("An unexpected error occurred", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        exit(1)

_setup_yaml()
if __name__ == '__main__':
    main(sys.argv[1:])