#!/usr/bin/env python
"""Creates config file, input directories, and validates a sample description file for a watermelon rnaseq job.

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

2) watermelon_init creates an analysis directory which contains a template
   watermelon config file. The template config file contains directory
   information, reference lists, run parameters, etc.,
   These must be reviewed and edited to:
   a) adjust genome and references to match the experiment
   b) adjust report_info with details from the experiment
   c) adjust trimming, alignment, and fastq_screen options
   d) specify diffex parameters, and add  DESeq2 calls and contrasts
      based on the example pheno.Gend stanza (can add as many similar stanzas
      as required by the experiment)

3) watermelon_init creates a readme file that lists basic info about when/how it was run,
   what it did, and what the user has to do to prepare the template config
"""
# pylint: disable=locally-disabled,no-member
from __future__ import print_function, absolute_import, division
import argparse
from collections import defaultdict, OrderedDict, Mapping
import datetime
import errno
import functools
import getpass
import glob
import gzip
import itertools
import os
import re
import ruamel_yaml
import shutil
import socket
import subprocess
import sys
import time
import traceback
import warnings
import yaml

import pandas as pd

from . import __version__ as WAT_VER
from tests import testing_utils # Using gunzip_glob, maybe should move this to helper module?

DESCRIPTION = (
    """Creates template config file and directories for a watermelon rnaseq job."""
)

_PATH_SEP = "^"


class _InputValidationError(Exception):
    """Raised for problematic input data."""

    def __init__(self, msg, *args):
        super(_InputValidationError, self).__init__(msg, *args)


class _UsageError(Exception):
    """Raised for malformed command or invalid arguments."""

    def __init__(self, msg, *args):
        super(_UsageError, self).__init__(msg, *args)


class _CommandValidator(object):
    # pylint: disable=locally-disabled,too-few-public-methods
    def __init__(self):
        pass

    def validate_args(self, args):
        # pylint: disable=locally-disabled, missing-docstring
        validations_list = [
            self._validate_overwrite_check,
            self._validate_genomebuild,
        ]
        if not args.count_matrix: # If count matrix argument is present, fastq dir validation is skipped
            validations_list.append(self._validate_source_fastq_dirs)

        for validation in validations_list:
            validation(args)

    def validate_inputs(self, input_summary):
        for validation in [
            self._validate_runs_have_fastq_files,
            self._validate_samples_have_fastq_files,
        ]:
            validation(input_summary)

    def validate_sample_sheet(self, args, input_summary):
        '''
        Verify that the samplesheet exists, that it contains a sample column, and that
        the sample names in the samplesheet correspond to those in the input directories
        '''
        sheet_file = args.sample_sheet

        if not os.path.isfile(sheet_file):
            msg = ("The specified sample sheet file {} cannot be found.".format(sheet_file))
            raise _InputValidationError(msg)

        try:
            samplesheet = pd.read_csv(sheet_file, comment='#', dtype='object').set_index(args.x_sample_column, drop=True)
        except KeyError:
            msg = "The sample sheet must have a column labeled '{}'".format(args.x_sample_column)
            raise _InputValidationError(msg)

        sheet_samples = set(samplesheet.index)
        input_samples = set(input_summary.samples)
        if not sheet_samples.issubset(input_samples):
            msg = ("The following sample names cannot be found in the input directories: "
            + str(sheet_samples.difference(input_samples)))
            raise _InputValidationError(msg)

    @staticmethod
    def _validate_genomebuild(args):
        with open(args.x_genome_references) as ref_yaml:
            ref_dict=yaml.load(ref_yaml, Loader=yaml.SafeLoader)

        genome_options = list(ref_dict.keys()) + ['TestData'] #Add TestData as valid options (special case which is set up before _make_config_dict)
        if args.genome_build not in genome_options:
            msg='genome {} is not found in {}.\nMust be one of {}'.format(args.genome_build, args.x_genome_references, [x for x in genome_options])
            raise _UsageError(msg)
        if args.genome_build in ['Other', 'TestData']: # This is a placeholder genome build, thus files don't need to exist. TestData isn't built until later
            pass
        else:
            reference_files = ref_dict[args.genome_build]['references'].values()
            if not any([os.path.isfile(x) for x in reference_files]):
                msg_fmt = "Could not read any of the following references for genome build {}\n{}"
                msg = msg_fmt.format(args.genome_build, ref_dict[args.genome_build]['references'])
                raise _InputValidationError(msg)



    @staticmethod
    def _validate_samples_have_fastq_files(input_summary):
        if input_summary.samples_missing_fastq_files:
            print(_build_input_summary_text(input_summary), file=sys.stderr)
            msg = (
                "Some samples missing fastq files: [{}]. "
                "Review details above against inputs and try again."
            ).format(", ".join(input_summary.samples_missing_fastq_files))
            raise _InputValidationError(msg)

    @staticmethod
    def _validate_runs_have_fastq_files(input_summary):
        if input_summary.runs_missing_fastq_files:
            print(_build_input_summary_text(input_summary), file=sys.stderr)
            msg = (
                "Some runs missing fastq files: [{}]. "
                "Review details above against inputs and try again."
            ).format(", ".join(input_summary.runs_missing_fastq_files))
            raise _InputValidationError(msg)

    @staticmethod
    def _validate_source_fastq_dirs(args):
        bad_dirs = [x for x in args.source_fastq_dirs if not os.path.isdir(x)]
        if bad_dirs:
            msg = (
                "Specified source_fastq_dir(s) [{}] not a dir or cannot be "
                "read. Review inputs and try again."
            ).format(",".join(sorted(bad_dirs)))
            raise _UsageError(msg)

    @staticmethod
    def _validate_overwrite_check(args):
        existing_files = {}
        if "input_dir" in vars(args) and os.path.exists(args.input_dir):
            existing_files["dirs: input"] = args.input_dir
        if "config_file" in vars(args) and os.path.exists(args.config_file):
            existing_files["config_file"] = args.config_file
        if os.path.exists("Watermelon"):
            existing_files["dirs: Watermelon"] = os.path.realpath("Watermelon")

        if existing_files:
            file_types = [file_type for file_type in sorted(existing_files.keys())]
            file_names = [existing_files[key] for key in file_types]
            msg_format = (
                "{file_types} [{file_names}] exists (will not overwrite). "
                "Remove these dir(s)/files(s) and try again."
            )
            msg = msg_format.format(
                file_types=",".join(file_types), file_names=",".join(file_names)
            )
            raise _UsageError(msg)


class _Linker(object):
    """Attempts to hard link file falling back to symlink as necessary.
    Only copies files."""

    def __init__(self):
        self.links = defaultdict(int)

    def _link(self, source, dest):
        try:
            os.link(source, dest)
            self.links["hardlinked"] += 1
        except OSError as e:
            if e.errno != errno.EXDEV:
                raise e
            os.symlink(source, dest)
            self.links["symlinked"] += 1

    def copytree(self, source, dest):
        shutil.copytree(source, dest, copy_function=self._link)

    @property
    def results(self):
        x = ", ".join(
            ["{} {}".format(count, type) for type, count in sorted(self.links.items())]
        )
        return "linked {} source files: {}".format(sum(self.links.values()), x)


def _timestamp():
    return datetime.datetime.fromtimestamp(time.time()).strftime("%Y-%m-%d %H:%M:%S")


README_FILENAME = "watermelon.README"

_SCRIPTS_DIR = os.path.realpath(os.path.dirname(__file__))
_WATERMELON_ROOT = os.path.dirname(_SCRIPTS_DIR)
_CONFIG_DIR = os.path.join(_WATERMELON_ROOT, "config")
_DEFAULT_TEMPLATE_CONFIG = os.path.join(_CONFIG_DIR, "template_config.yaml")
_DEFAULT_SAMPLE_COLUMN = "sample"
_DEFAULT_GENOME_REFERENCES = os.path.join(_CONFIG_DIR, "genome_references.yaml")

try:
    _DEFAULT_ANALYST_INFO = pd.read_csv('/nfs/turbo/umms-brcfpipeline/pipelines/analyst_info.csv', comment="#", index_col='username') # This will be accessible on greatlakes
except:
    msg = f"\nCould not read analyst_info csv. Using default placeholder values"
    warnings.warn(msg)
    _DEFAULT_ANALYST_INFO = pd.DataFrame() # An empty dataframe used as a simple check below

_FASTQ_GLOBS = ["*.fastq", "*.fastq.gz"]

_TODAY = datetime.date.today()
_DEFAULT_JOB_SUFFIX = "_{:02d}_{:02d}".format(_TODAY.month, _TODAY.day)
_CONFIG_PRELUDE = """# config created {timestamp}
# To do:
# ------
# 1) Review genome and references
# 2) Review report information
# 3) Review alignment, trimming, and fastq_screen options
# 4) Modify the diffex options - use the 'model_one' stanza as an example, and
#    set up the DESeq2 calls and results calls for the comparisons that match the analysis.
""".format(
    timestamp=_timestamp()
)


def _setup_yaml():
    """ http://stackoverflow.com/a/8661021 """
    represent_dict_order = lambda self, data: self.represent_mapping(
        "tag:yaml.org,2002:map", data.items()
    )
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
        raise OSError(
            "a file with the same name as the desired "
            "dir, '%s', already exists." % newdir
        )
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)


def _copy_and_overwrite(source, dest):
    if os.path.exists(dest):
        shutil.rmtree(dest)
    source = os.path.join(source, "")  # Add trailing slash for rsync's sake
    # Exclude .git/ and /envs/built (speed)
    subprocess.run(
        ["rsync", "-O", "-rlt", "--exclude", ".*", "--exclude", "envs/built", source, dest]
    )

def _setup_test_data(datapath):
    #Unzip the human chr22 example references
    data_gz_glob = os.path.join(datapath , '*.gz')
    testing_utils.gunzip_glob(data_gz_glob)
    #These refs should now exist. If not, it'll be caught during snakemake run w/ init_references
    refs = {
        'genome' : 'GRCh38',
        'references' : {
            'fasta' : os.path.join(datapath, 'Homo_sapiens.GRCh38.dna_sm.chr22.fa'),
            'gtf' : os.path.join(datapath, 'Homo_sapiens.GRCh38.98.chr22.gtf'),
            'annotation_tsv' : os.path.join(datapath, 'Homo_sapiens.GRCh38.98_annotation.tsv')
        }
    }
    return refs


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
        if k in dct and isinstance(dct[k], dict) and isinstance(merge_dct[k], Mapping):
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
    # TWS - FIXME: This is inflexible, prone to breakage, ugly, etc.
    j = os.path.join
    watermelon_runs_dir = os.path.abspath(j(args.tmp_input_dir, args.input_runs_dir))
    watermelon_samples_dir = os.path.abspath(
        j(args.tmp_input_dir, args.input_samples_dir)
    )
    _mkdir(watermelon_samples_dir)
    for run_dir in [
        d
        for d in os.listdir(watermelon_runs_dir)
        if os.path.isdir(j(watermelon_runs_dir, d))
    ]:
        for sample_dir in [
            d
            for d in os.listdir(j(watermelon_runs_dir, run_dir))
            if os.path.isdir(j(watermelon_runs_dir, run_dir, d))
        ]:
            new_sample_dir = j(watermelon_samples_dir, sample_dir)
            _mkdir(new_sample_dir)
            files = os.listdir(j(watermelon_runs_dir, run_dir, sample_dir))
            if files:
                for sample_file in files:
                    source = j(watermelon_runs_dir, run_dir, sample_dir, sample_file)
                    link = j(new_sample_dir, run_dir + _PATH_SEP + sample_file)
                    os.link(source, link)
            else:
                from pathlib import Path

                Path(j(new_sample_dir, run_dir + _PATH_SEP + "empty_dir")).touch()


def _is_fastq_file(file_name):
    return file_name.endswith(".fastq") or file_name.endswith(".fastq.gz")


def _build_input_summary(args):
    input_samples_dir = os.path.join(args.tmp_input_dir, args.input_samples_dir)
    rows = []
    for root, dirs, files in os.walk(input_samples_dir):
        for name in files:
            file_name = os.path.basename(name)
            sample_name = os.path.basename(root)
            run_name = os.path.dirname(file_name.replace(_PATH_SEP, os.path.sep))
            rows.append((run_name, sample_name, file_name))

    df = pd.DataFrame(data=rows, columns=["run", "sample", "file"])
    df["run_index"] = df["run"].rank(method="dense").apply(lambda x: chr(64 + int(x)))
    df["is_fastq_file"] = df["file"].apply(_is_fastq_file)

    runs_df = (
        df.drop_duplicates(subset=["run_index", "run"])[["run_index", "run"]]
        .set_index("run_index")
        .sort_index(axis=0)
    )
    sample_run_counts_df = df.pivot_table(
        values="is_fastq_file", index=["sample"], columns="run_index", aggfunc="sum"
    )
    sample_run_counts_df = sample_run_counts_df.fillna(0).astype(int, errors="ignore")

    total_sample_count = len(sample_run_counts_df)
    total_file_count = sum(sum(sample_run_counts_df.values))
    samples = sorted(sample_run_counts_df.index)
    sample_run_counts_df.loc["run_total"] = sample_run_counts_df.sum()
    sample_run_counts_df["sample_total"] = sample_run_counts_df.sum(axis=1)
    samples_missing_fastq_files = sorted(
        sample_run_counts_df[sample_run_counts_df["sample_total"] == 0].index
    )

    run_counts_df = sample_run_counts_df[
        sample_run_counts_df.index == "run_total"
    ].transpose()
    runs_missing_fastq_files = sorted(
        run_counts_df[run_counts_df["run_total"] == 0].index
    )

    input_summary = argparse.Namespace(
        runs_df=runs_df,
        sample_run_counts_df=sample_run_counts_df,
        total_sample_count=total_sample_count,
        total_file_count=total_file_count,
        total_run_count=len(runs_df),
        samples=samples,
        samples_missing_fastq_files=samples_missing_fastq_files,
        runs_missing_fastq_files=runs_missing_fastq_files,
    )
    return input_summary

def _make_config_dict(template_config, genome_references, args):
    # Due to ruamel_yaml, template_config is an OrderedDict, and has preserved comments
    config = template_config
    # Make changes to existing data
    # Fill in job suffix for output dirs
    for out_dir in config["dirs"]:
        config["dirs"][out_dir] = config["dirs"][out_dir].format(
            job_suffix=args.job_suffix
        )
    # If "Other" or "TestData" genome_build is specified, remove fastq_screen stanza from config
    if args.genome_build == "Other" or args.genome_build == "TestData":
        config.pop('fastq_screen', None)

    # add input dir
    config["dirs"]["input"] = os.path.join(args.input_dir, args.input_samples_dir)
    # Add in more needed keys
    # Watermelon version
    config["watermelon_version"] = WAT_VER
    # Samplesheet
    samplesheet_path = os.path.abspath(args.sample_sheet)
    samplesheet_path = re.sub('/ccmb/BioinfCore/ActiveProjects/', '/nfs/med-bfx-activeprojects/', samplesheet_path) #TODO: Can remove this line after comp5/6 mounts are fixed

    config["samplesheet"] = samplesheet_path

    # Genome / references
    _dict_merge(config, genome_references)
    # Email params
    user = getpass.getuser()
    config["email"] = {
        "subject": "watermelon" + args.job_suffix,
        "to": f"{user}@umich.edu",
    }

    # Analyst name
    if config.get("report_info"):
        try:
            config["report_info"]["analyst_name"] = _DEFAULT_ANALYST_INFO.at[user, 'name']
        except KeyError:
            # _DEFAULT_ANALYST_INFO may be empty (default file not available)
            # Or the user is not in that DataFrame
            config["report_info"]["analyst_name"] = 'Analyst'


    # If count matrix argument given, Add count matrix to config, remove all unnecessary stuff
    if args.count_matrix:
        config["count_matrix"] = args.count_matrix
        config.move_to_end("count_matrix", last=False)
        config.pop("fastq_screen", None)
        config.pop("trimming_options", None)
        config.pop("alignment_options", None)
        config["dirs"].pop("alignment_output", None)
        config["dirs"].pop("input", None)
        config["references"].pop("fasta", None)
        config["references"].pop("gtf", None)


    # Reorder these added keys (move to top)
    # resulting order: email, watermelon_version, samplesheet, genome, references, template_config stuff
    config.move_to_end("references", last=False)
    config.move_to_end("genome", last=False)
    config.move_to_end("samplesheet", last=False)
    config.move_to_end("watermelon_version", last=False)
    config.move_to_end("email", last=False)

    return config


def _build_input_summary_text(input_summary):
    pd.set_option("display.max_colwidth", 300)
    text = """Input summary:
--------------
Found {total_file_count} files for {total_sample_count} samples across {total_run_count} run(s).

Run dir(s):
-----------
{run_details}

Sample x run fastq file counts:
-------------------------------
{sample_run_counts}
""".format(
        total_file_count=input_summary.total_file_count,
        total_sample_count=input_summary.total_sample_count,
        total_run_count=input_summary.total_run_count,
        run_details=input_summary.runs_df.to_string(justify="left"),
        sample_run_counts=input_summary.sample_run_counts_df.to_string(justify="left"),
    )
    return text


def _build_postlude(args, linker_results, input_summary_text):
    input_dir_relative = os.path.relpath(args.input_dir, args.x_working_dir)
    postlude = """
watermelon_init.README
======================
watermelon_init invocation:
{watermelon_init_invocation}

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
2) Review report information
3) Review alignment, trimming, and fastq_screen options
4) Modify the diffex options - use the 'model_one' stanza as an example, and
   set up the DESeq2 calls and results calls for the comparisons that match the analysis.

When the config & samplesheet look good:
--------------------------------
# Start a screen session:
$ screen -S watermelon{job_suffix}
# Activate the conda environment:
$ conda activate watermelon
# To validate the config and check the execution plan:
$ snakemake --dryrun --printshellcmds --configfile {config_basename} --snakefile {snakefile_path}
#
# To run on comp5/6:
# Singularity must be available to snakemake, for environment management under the hood
$ module load singularity/3.5.2
$ snakemake --configfile {config_basename} --snakefile {snakefile_path} --profile {profile_path_comps}
#
# To run on the greatlakes cluster:
# Singularity must be available to snakemake, for environment management under the hood
$ module load singularity/3.5.2
$ snakemake --configfile {config_basename} --snakefile {snakefile_path} --profile {profile_path_gl}
""".format(
        watermelon_init_invocation=" ".join(sys.argv[:]),
        linker_results=linker_results,
        input_summary_text=input_summary_text,
        working_dir=args.x_working_dir,
        input_dir_relative=input_dir_relative,
        input_runs_dir=args.input_runs_dir,
        input_samples_dir=args.input_samples_dir,
        config_file=args.config_file,
        config_relative=os.path.relpath(args.config_file, args.x_working_dir),
        samplesheet_relative=os.path.relpath(args.sample_sheet, args.x_working_dir),
        config_basename=os.path.basename(args.config_file),
        job_suffix=args.job_suffix,
        snakefile_path="Watermelon/rnaseq.snakefile",
        profile_path_comps="Watermelon/config/profile-comp5-6",
        profile_path_gl="Watermelon/config/profile-greatlakes"
    )
    return postlude


def _write_config_file(config_filename, config_dict):
    tmp_config_dict = dict(config_dict)
    ordered_config_dict = OrderedDict()
    named_keys = [
        "dirs",
        "main_factors",
        "phenotypes",
        "samples",
        "comparisons",
        "genome",
        "references",
    ]
    for key in named_keys:
        if key in tmp_config_dict:
            ordered_config_dict[key] = tmp_config_dict.pop(key)
    for key in sorted(tmp_config_dict):
        ordered_config_dict[key] = tmp_config_dict[key]
    with open(config_filename, "w") as config_file:
        print(_CONFIG_PRELUDE, file=config_file)
        yaml.dump(ordered_config_dict, config_file, default_flow_style=False, indent=4)


def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--genome_build",
        type=str,
        help="Config template will based on this genome",
        required=True,
    )
    parser.add_argument(
        "--job_suffix",
        type=str,
        default=_DEFAULT_JOB_SUFFIX,
        help=(
            "={} This suffix will be appended to the names of analysis/deliverables "
            "dirs and configfile. "
            "(Useful if rerunning a job with revised input fastq or revised config; "
            "can help differentiate simultaneous watermelon jobs in "
            "top/ps.)"
        ).format(_DEFAULT_JOB_SUFFIX),
        required=False,
    )
    parser.add_argument(
        "--sample_sheet",
        type=str,
        required=True,
        help=(
            "A CSV file containing sample names and phenotype/characteristic "
            "information which correspond to these samples. Watermelon_init will verify that "
            "sample names in this file match with those found in the input directories. "
        ),
    )
    parser.add_argument(
        "--count_matrix",
        type=str,
        default='',
        help=(
            "A count matrix file (Optional). If a count matrix file is supplied, watermelon_init will "
            "create a configuration file suitable for such an analysis. Watermelon_init will verify that "
            "column names in the count matrix match the samples in the sample sheet. "
        ),
    )
    parser.add_argument(
        "--AGC",
        default=False,
        action='store_true',
        help=(
            "An optional flag which enables the creation of a configuration file that can be directly used by the UMich Advanced Genomics Core without modification."
        ),
    )
    parser.add_argument(
        "source_fastq_dirs",
        type=str,
        nargs="*",
        default=None,
        help=(
            "One or more paths to run dirs. Each run dir should contain "
            "samples dirs which should in turn contain one or more fastq.gz "
            "files. The watermelon sample names will be derived from the "
            "sample directories."
        ),
    )
    parser.add_argument(
        "--x_working_dir", type=str, default=os.getcwd(), help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--x_template_config",
        type=str,
        default=_DEFAULT_TEMPLATE_CONFIG,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--x_genome_references",
        type=str,
        default=_DEFAULT_GENOME_REFERENCES,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--x_analyst_info",
        type=str,
        default=_DEFAULT_ANALYST_INFO,
        help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--x_sample_column",
        type=str,
        default=_DEFAULT_SAMPLE_COLUMN,
        help=argparse.SUPPRESS,
    )

    args = parser.parse_args(sys_argv)

    args.x_working_dir = re.sub('/ccmb/BioinfCore/ActiveProjects/', '/nfs/med-bfx-activeprojects/', args.x_working_dir) #TODO: Can remove this line after comp5/6 mounts are fixed

    realpath = functools.partial(os.path.join, args.x_working_dir)
    args.config_file = "config{}.yaml".format(args.job_suffix)
    args.tmp_input_dir = realpath(".inputs")
    args.input_runs_dir = "00-source_runs"
    args.input_samples_dir = "01-source_samples"
    args.input_dir = realpath("inputs")
    return args


def main(sys_argv):
    """See DESCRIPTION"""
    try:
        args = _parse_command_line_args(sys_argv)
        _CommandValidator().validate_args(args)

        print("Copying Watermelon source to " + os.getcwd())
        _copy_and_overwrite(
            source=_WATERMELON_ROOT, dest=os.path.join(os.getcwd(), "Watermelon")
        )

        if not args.count_matrix:
            print("Linking inputs")
            if os.path.isdir(args.tmp_input_dir):
                shutil.rmtree(args.tmp_input_dir)
            linker = _Linker()
            _link_run_dirs(args, linker)
            _merge_sample_dirs(args)
            input_summary = _build_input_summary(args)
            _CommandValidator().validate_inputs(input_summary)
            os.rename(args.tmp_input_dir, args.input_dir)
        else: # TWS FIXME: This is a temporary solution for differing procedures for watermelon_init. Refactor for ovarhaul
            samplesheet = pd.read_csv(args.sample_sheet, comment='#', dtype='object').set_index("sample", drop=True)
            samples = list(samplesheet.index)
            input_summary = argparse.Namespace(
                runs_df=None,
                sample_run_counts_df=None,
                total_sample_count=None,
                total_file_count=None,
                total_run_count=None,
                samples=samples,
                samples_missing_fastq_files=None,
                runs_missing_fastq_files=None,
            )

        _CommandValidator().validate_sample_sheet(args, input_summary)
        print("Generating example config")
        with open(args.x_template_config, 'r') as template_config_file:
            template_config = ruamel_yaml.round_trip_load(template_config_file) #Use ruamel_yaml to preserve comments
        with open(args.x_genome_references, 'r') as genome_references_file:
            all_genome_refs = yaml.load(genome_references_file, yaml.SafeLoader)

        if args.genome_build == "TestData":
            test_data_dir = os.path.join(args.x_working_dir, 'Watermelon', 'data')
            selected_genome_refs = _setup_test_data(test_data_dir)
        else:
            selected_genome_refs = all_genome_refs[args.genome_build]

        config_dict = _make_config_dict(template_config,
                                        selected_genome_refs,
                                        args)

        if args.AGC:
            # If AGC flag is passed, remove diffex section,
            # In report, change Bioinformatics Core to Advanced Genomics Core,
            # Don't include follow-up text
            config_dict.pop("diffex", None)
            curr_txt = config_dict["report_info"]["acknowledgement_text"]
            new_txt = re.sub("Bioinformatics", "Advanced Genomics", curr_txt)
            config_dict["report_info"]["acknowledgement_text"] = new_txt
            config_dict["report_info"]["include_follow_up"] = False

        # _write_config_file(args.config_file, config_dict)
        with open(args.config_file, "w") as config_file:
            print(_CONFIG_PRELUDE, file=config_file)
            ruamel_yaml.round_trip_dump(
                config_dict, config_file, default_flow_style=False, indent=4
            )

        if not args.count_matrix:
            postlude = _build_postlude(
                args, linker.results, _build_input_summary_text(input_summary)
            )
        else:
            postlude = _build_postlude(args, "", "")
        print(postlude)
        with open(README_FILENAME, "w") as readme:
            print(postlude, file=readme)
    except _UsageError as usage_error:
        message = "watermelon_init usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'watermelon_init --help'.", file=sys.stderr)
        sys.exit(1)
    except _InputValidationError as input_validation_error:
        print(input_validation_error, file=sys.stderr)
        sys.exit(1)
    except Exception:  # pylint: disable=locally-disabled,broad-except
        print("An unexpected error occurred", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        exit(1)


_setup_yaml()
if __name__ == "__main__":
    main(sys.argv[1:])
