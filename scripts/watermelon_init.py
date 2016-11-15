#!/usr/bin/env python
from __future__ import print_function, absolute_import, division
import argparse
from collections import OrderedDict
import datetime
import functools
import glob
import os
import sys

import yaml

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

def setup_yaml():
  """ http://stackoverflow.com/a/8661021 """
  represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
  yaml.add_representer(OrderedDict, represent_dict_order)    
setup_yaml()

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
Created files and dirs:
    {inputs_dir}
        source fastq dir | sample count | fastq file count
        {source_fastq_dir} | {sample_count} samples | {file_count} files
    {analysis_dir}
    {deliverables_dir}
    {config_file}

You need to review config file: {config_file}:
1) Review/adjust samples names
    Note: if you change names in config, also change the sample dir names in the input dir
2) Add a sample group for each sample
3) Add comparisons
4) Review genome and references
5) Review alignment options
6) Review trimming options

When the config file looks good:
$ cd {analysis_dir}
# start a screen session:
$ screen -S watermelon{job_suffix}
# to validate the config and check the execution plan:
$ watermelon --dry-run {config_basename}
# to run:
$ watermelon -c {config_basename}
'''.format(inputs_dir=args.inputs_dir,
           source_fastq_dir=args.source_fastq_dir,
           sample_count=sample_count,
           file_count=file_count,
           analysis_dir=args.analysis_dir,
           deliverables_dir=args.deliverables_dir,
           config_file=args.config_file,
           config_basename=os.path.basename(args.config_file),
           job_suffix=args.job_suffix,)
    return postlude

def _make_top_level_dirs(args):
    _mkdir(args.analysis_dir)
    _mkdir(args.deliverables_dir)
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
    args.deliverables_dir = realpath('deliverables{}'.format(args.job_suffix))
    args.inputs_dir = realpath('inputs', '00-multiplexed_reads')

    return args

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    _make_top_level_dirs(args)
    (samples, file_count) = _initialize_samples(args.source_fastq_dir)
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
    print(_build_postlude(args, len(samples), file_count))


if __name__ == '__main__':
    main(sys.argv[1:])