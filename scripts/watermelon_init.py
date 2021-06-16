#!/usr/bin/env python3
'''
Produces a config.yaml to use with snakemake

Intakes the following:
* sample sheet. Contains sample IDs, paths of samples' fastq files, and columns for phenotype information
* genome build. Name of a set of reference files (e.g. GRCh38)
* project ID. For BFX this is project name, for AGC this is service request ID. No spaces
* optional output prefix. Defaults to relative path analysis_projectID/{alignment_results,deliverables,report}
* config type. align_qc or diffex. Creates a different config depending on this

Determines at runtime:
* analyst name / email. For AGC, use the organization as the name. For emailing, send to the individual
* Watermelon version
* References are readable
* README.txt is readable. If not, use templated default report description text
* Input dirs contain fastqs that are readable

Uses template values for:
* Trimming options
* Fastq_screen options
* Default report description (if README.txt not available)
'''

# Need to add/test:

# Analyst info check (csv available or not)
# Error handling if fail loading template
# Test an invalid key (genome build) in set up refs

# Validations:
# Fastqs in sample dirs, readable
# Genome build is valid, files are readable

# Don't think it makes sense to require sample dir name matches samplesheet
# anymore since the paths are given


import argparse
import copy
import getpass
import os
import ruamel_yaml
import subprocess
import sys
import warnings

import pandas as pd
from collections import OrderedDict
from collections.abc import Mapping

import rnaseq_snakefile_helper as helper

import pdb # TWS DEBUG


# Set default values
# Find path knowing current location of Watermelon/scripts/watermelon_init.py
_WATERMELON_ROOT = os.path.realpath(os.path.dirname(os.path.dirname(__file__)))
_DEFAULT_GENOME_REFERENCES = os.path.join(_WATERMELON_ROOT, "config", "genome_references.yaml")
_DEFAULT_SAMPLE_COL = "sample"

try:
    _DEFAULT_ANALYST_INFO = pd.read_csv("/nfs/turbo/umms-brcfpipeline/pipelines/analyst_info.csv", comment="#", index_col="username") # This will be accessible on greatlakes
except:
    msg = f"\nCould not read analyst_info csv. Using default placeholder values"
    warnings.warn(msg)
    _DEFAULT_ANALYST_INFO = pd.DataFrame() # An empty dataframe used as a simple check below


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


def _no_overwrite_check(check_vals):
    exists = []
    for val in check_vals:
        if os.path.exists(val):
            exists.append(val)
    if exists:
        msg_fmt = "\n\nThe following file(s)/dir(s) exist (will not overwrite). Remove these and try again:\n{}\n"
        msg = msg_fmt.format(exists)
        raise RuntimeError(msg)


def _set_up_dirs(type, project_id, fmt_str="analysis_{projid}/{resulttype}"):
    dirs = ruamel_yaml.comments.CommentedMap() # Use this instead of an OrderedDict since it prints pretty (no !!omap)
    if type == "align_qc":
        for resulttype in ["alignment_results", "deliverables", "report"]:
            dirs[resulttype] = fmt_str.format(projid=project_id, resulttype=resulttype)
    elif type == "diffex":
        for resulttype in ["diffex_results", "deliverables", "report"]:
            dirs[resulttype] = fmt_str.format(projid=project_id, resulttype=resulttype)
    return dirs


def _set_up_email(projid, user):
    email = {
        "subject": "watermelon_" + projid,
        "to": "{}@umich.edu".format(user),
    }
    return email


def _set_up_refs(grefsfn, gbuild, type):
    with open(grefsfn, "r") as reffh:
        refs = ruamel_yaml.round_trip_load(reffh)
    try:
        our_refs = refs[gbuild]
    except KeyError:
        msg_fmt = "\n\nGenome build {} not found in {}:\n {}\n"
        raise UsageError(msg_fmt.format(gbuild, grefsfn, refs.keys()))
    # If diffex type, only need annotation TSV. If align_qc, keep it all
    if type == "diffex":
        for ref in ["fasta", "gtf", "rsem_star_index"]:
            our_refs["references"].pop(ref)
    return our_refs


def get_pipe_version(vfn):
    with open(vfn, "r") as vfh:
        # version_info.smk gives access to VER_INFO with software version info
        # and also to ENV_INFO - a singularity info dict (not currently used in init)
        WORKFLOW_BASEDIR = _WATERMELON_ROOT # WORKFLOW_BASEDIR needed in snakefile
        _locals = locals() # capturing locals gives access to WORKFLOW_BASEDIR in snakefile
        exec(vfh.read(), _locals, globals())
        version = VER_INFO["watermelon"]
        return version


def get_template(args, pipe_root):
    if args.x_alt_template:
        warnings.warn("Using alternate template given. Ignoring type")
        fpath = args.x_alt_template
    else:
        fname = "template_{}.yaml".format(args.type)
        fpath = os.path.join(pipe_root, "config", fname)
    with open(fpath, "r") as tfh:
        template_dict = ruamel_yaml.round_trip_load(tfh)
        return template_dict


def make_config_dict(template_config, args, version):
    # Due to ruamel_yaml, template_config is an OrderedDict, and has preserved comments
    config = template_config

    # If "Other" or "TestData" genome_build is specified, remove fastq_screen stanza from config
    if args.genome_build == "Other" or args.genome_build == "TestData":
        config.pop("fastq_screen", None)

    # Set up some things before adding them
    refs = _set_up_refs(args.x_genome_references, args.genome_build, args.type)
    dirs = _set_up_dirs(args.type, args.project_id)
    user = getpass.getuser() # used for email and for analyst name below
    email = _set_up_email(args.project_id, user)

    # Merge the genome reference dict with the config
    # Inserts top-level keys "genome" & "references" into config
    # Integrates "fastq_screen" ref vals with template vals
    _dict_merge(config, refs)
    # desired order: email, watermelon_version, samplesheet, genome, references, dirs, template_config stuff
    # Focus on the 2nd half (working backwards)
    config.insert(0, "dirs", dirs) # Insert dirs on top
    config.move_to_end("references", last=False) # This is how to move existing keys to the top
    config.move_to_end("genome", last=False)
    # 1st half (working backwards)
    # If non-default sample col is given, add it to the config
    if args.x_sample_col != _DEFAULT_SAMPLE_COL:
        config.insert(0, "sample_col", args.x_sample_col)
    config.insert(0, "samplesheet", os.path.abspath(args.sample_sheet))
    config.insert(0, "watermelon_version", version)
    config.insert(0, "email", email)

    # Add analyst name to report info
    if config.get("report_info"):
        try:
            config["report_info"]["analyst_name"] = _DEFAULT_ANALYST_INFO.at[user, "name"]
        except KeyError:
            # _DEFAULT_ANALYST_INFO may be empty (default file not available)
            # Or the user is not in that DataFrame
            config["report_info"]["analyst_name"] = "Analyst"

    return config

def validate_fastq_dirs(config_dict, sample_col, fq_col):
    with open(config_dict["samplesheet"], "r") as ssfh:
        samplesheet = pd.read_csv(ssfh, comment="#", dtype="object")
        if not sample_col in samplesheet.columns:
            msg = "The sample sheet must have a column labeled '{}'".format(sample_col)
            raise RuntimeError(msg)
        if not fq_col in samplesheet.columns:
            msg = "\n\nThe sample sheet must have a column labeled '{}'\n".format(fq_col)
            raise RuntimeError(msg)

        samples = samplesheet[sample_col]
        fq_dirs = samplesheet[fq_col]
        fastq_dict = OrderedDict()
        for sample, dir in zip(samples, fq_dirs):
            if os.path.exists(dir):
                # The helper function has the following validations:
                # * Fastq file(s) must be present
                # * Can't mix gz and fastq files (must be all same type)
                fastq_dict[sample] = helper.get_sample_fastq_paths(dir)
                # TODO: could easily write this to a file...
            else:
                msg = "\n\nDirectory {} does not contain fastq files.\n".format(dir)
                raise RuntimeError(msg)


def validate_genomes(ref_dict, type):
    # Validate that references are readable
    # Make a copy of the dict (don't actually want to manipulate it)
    ref_dict = copy.copy(ref_dict)
    # Diffex only uses annotation_tsv, not the others
    if type == "diffex":
        for k in ["fasta", "gtf", "rsem_star_index"]:
            ref_dict.pop(k)
    for k,v in ref_dict.items():
        cant_read = []
        if k == "rsem_star_index":
            v = v + ".chrlist" # Given index prefix, test for chromosome list
        if not os.access(v, os.R_OK): # Test for read access
            cant_read.append(v)
        if cant_read:
            msg = "\n\nCan't read references:\n{}\n"
            raise RuntimeError(msg)


def write_stuff(config_dict, config_fn, wat_dir=_WATERMELON_ROOT):
    # Make sure we don't overwrite things
    _no_overwrite_check(["Watermelon", config_fn])
    print("Writing config to working dir...")
    # Write the config file
    with open(config_fn, "w") as config_fh:
        ruamel_yaml.round_trip_dump(
            config_dict, config_fh, default_flow_style=False, indent=4
        )
    # Copy Watermelon to project dir
    print("Done.\nCopying Watermelon to working dir...")
    source = os.path.join(wat_dir, "")  # Add trailing slash for rsync (syntax important)
    dest = "Watermelon"
    subprocess.run( # Exclude .git/ and /envs/built (speed)
        ["rsync", "-O", "-rlt", "--exclude", ".*", "--exclude", "envs/built", source, dest]
    )

    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="watermelon_init.py", description = "Produces a config.yaml to use with snakemake. Can produce a config for either alignment/QC or differential expression (see below).")
    parser.add_argument("-g", "--genome_build", required=True, help="Reference files for this build are added into the config.")
    parser.add_argument("-o", "--output_prefix", default="analysis_{project_id}", help="Optional output prefix. Defaults to relative path %default%.")
    parser.add_argument("-p", "--project_id", required=True, help="ID used for the project. Can represent project ID, service request ID, etc. Only alphanumeric and underscores allowed. No spaces or special characters.")
    parser.add_argument("-s", "--sample_sheet", required=True, help="CSV file that contains sample IDs, paths containing samples\" fastq files, and optional columns for phenotype information (required for diffex).")
    parser.add_argument("-t", "--type", default="align_qc", choices=["align_qc", "diffex"], help="Type of config file to produce. Uses the appropriate template, and has the appropriate values.")
    parser.add_argument("--x_analyst_info", type=str, default=_DEFAULT_ANALYST_INFO, help=argparse.SUPPRESS)
    parser.add_argument("--x_alt_template", type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--x_genome_references", type=str, default=_DEFAULT_GENOME_REFERENCES, help=argparse.SUPPRESS)
    parser.add_argument("--x_working_dir", type=str, default=os.getcwd(), help=argparse.SUPPRESS)
    parser.add_argument("--x_sample_col", type=str, default="sample", help=argparse.SUPPRESS)
    parser.add_argument("--x_input_col", type=str, default="input_dir", help=argparse.SUPPRESS)

    args = parser.parse_args()

    # Which template depends on type of config (args.type)
    template_cfg_dict = get_template(args, _WATERMELON_ROOT)

    # Ingest version_info.smk snakefile to get pipeline version
    version = get_pipe_version(os.path.join(_WATERMELON_ROOT, "version_info.smk"))

    # Make the config
    config_dict = make_config_dict(template_cfg_dict, args, version)

    # Perform config validations
    validate_genomes(config_dict["references"], args.type)
    # Some validations depend on type of config
    if args.type == "align_qc":
        validate_fastq_dirs(config_dict, args.x_sample_col, args.x_input_col)
    elif args.type == "diffex":
        foo = "bar"

    write_stuff(
        config_dict=config_dict,
        config_fn="config_{}.yaml".format(args.project_id),
        wat_dir=_WATERMELON_ROOT
    )
