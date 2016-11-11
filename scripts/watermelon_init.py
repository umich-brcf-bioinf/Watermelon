#!/usr/bin/env python
from __future__ import print_function, absolute_import, division
import argparse
import datetime
import functools
import os
import sys

DESCRIPTION = \
'''Creates template config file and directories for a watermelon rnaseq job.'''

_TODAY = datetime.date.today()
_DEFAULT_JOB_SUFFIX = '_{:02d}_{:02d}'.format(_TODAY.month, _TODAY.day)
# TODO: cgates: pull this from config dir?
_GENOME_BUILD_OPTIONS = ('hg19', 'mm10', 'rn5')

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

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--genome_build',
        type=str,
        choices=_GENOME_BUILD_OPTIONS,
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
        '--working_dir',
        type=str,
        default=os.getcwd(),
        help=argparse.SUPPRESS)

    args = parser.parse_args(sys_argv)

    realpath = functools.partial(os.path.join, args.working_dir)
    args.analysis_dir = realpath('analysis{}'.format(args.job_suffix))
    args.configfile = os.path.join(args.analysis_dir,
                                   'config{}.yaml'.format(args.job_suffix))
    args.deliverables_dir = realpath('deliverables{}'.format(args.job_suffix))
    args.inputs_dir = realpath('inputs', '00-multiplexed_reads')

    return args

def _build_postlude(args, sample_count, file_count):
    postlude = \
'''
Created:
    {inputs_dir}
        source fastq dir | sample count | file count
        {source_fastq_dir} | {sample_count} samples | {file_count} files
    {analysis_dir}
    {deliverables_dir}
You need to review {configfile}:
1) Review/adjust samples names
    Note: if you change names in config, also change the sample dir names in the input dir
2) Add a sample group for each sample
3) Add comparisons
4) Review alignment options
5) Review trimming options

When the configfile looks good:
$ cd {analysis_dir}
$ screen -S watermelon{job_suffix} #starts a screen session
$ watermelon --dry-run {configfile} # to validate the config and check the execution plan
$ watermelon {configfile} # to run
'''.format(inputs_dir=args.inputs_dir,
           source_fastq_dir=args.source_fastq_dir,
           sample_count=sample_count,
           file_count=file_count,
           analysis_dir=args.analysis_dir,
           deliverables_dir=args.deliverables_dir,
           configfile=args.configfile,
           job_suffix=args.job_suffix,)
    return postlude

def _make_top_level_dirs(args):
    _mkdir(args.analysis_dir)
    _mkdir(args.deliverables_dir)
    _mkdir(args.inputs_dir)

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)

if __name__ == '__main__':
    main(sys.argv[1:])