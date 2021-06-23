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

# Validations:
# Fastqs in sample dirs, readable

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
from natsort import natsorted

try:
    import scripts.rnaseq_snakefile_helper as helper # works in pytest
except:
    import rnaseq_snakefile_helper as helper # works when called as script

import pdb # TWS DEBUG


# Set default values
# Find path knowing current location of Watermelon/scripts/watermelon_init.py
_WATERMELON_ROOT = os.path.realpath(os.path.dirname(os.path.dirname(__file__)))
_DEFAULT_GENOME_REFERENCES = os.path.join(_WATERMELON_ROOT, "config", "genome_references.yaml")
_DEFAULT_SAMPLE_COL = "sample"
_DEFAULT_ANALYST_INFO = "/nfs/turbo/umms-brcfpipeline/pipelines/analyst_info.csv"


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


def _prompt_yes_no(prompt):
    while True:
        value = input(prompt).lower().strip()
        if value not in ['yes', 'y', 'no', 'n']:
            print("Response must one of: yes, y, no, n")
            continue
        else:
            break
    return value in ['yes', 'y']


def _no_overwrite_check(check_vals):
    exists = []
    for val in check_vals:
        if os.path.exists(val):
            exists.append(val)
    if exists:
        prompt = "The following file(s)/dir(s) exist:\n{}\nOverwrite? (yes/no): ".format(exists)
        if not _prompt_yes_no(prompt):
            msg = "\n\nThe following file(s)/dir(s) exist and will not overwrite:\n{}\n".format(exists)
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
        raise RuntimeError(msg_fmt.format(gbuild, grefsfn, refs.keys()))
    # If diffex type, only need annotation TSV. If align_qc, keep it all
    if type == "diffex":
        for ref in ["fasta", "gtf", "rsem_star_index"]:
            our_refs["references"].pop(ref)
    return our_refs


def _get_analyst_name(analyst_info_csv, user):
    try:
        analyst_info = pd.read_csv(analyst_info_csv, comment="#", index_col="username") # This will be accessible on greatlakes
        analyst_name = analyst_info.at[user, "name"]
    except OSError:
        msg = "\n\nCould not read analyst_info csv. Using default placeholder values.\n"
        warnings.warn(msg)
        analyst_name = "Analyst"
    except KeyError:
        msg = f"\n\nCouldn't find {user} in {analyst_info_csv}. Using default placeholder values.\n"
        warnings.warn(msg)
        analyst_name = "Analyst"

    return analyst_name


def _generate_samplesheet(fq_parent_dirs, fname):
    print("Generating sample sheet {} ...".format(fname))
    samples_dict = OrderedDict()
    for parent_dir in fq_parent_dirs:
        if not os.path.exists(parent_dir):
            msg = "\n\nInput run directory {} does not exist.\n".format(parent_dir)
            raise RuntimeError(msg)
        # list of tuples
        sample_dirs = [(d.name, os.path.abspath(d)) for d in os.scandir(parent_dir) if d.is_dir()]
        sample_dirs = natsorted(sample_dirs)
        if not sample_dirs:
            msg = "\n\nInput run directoy {} contains no subdirectories.\n".format(parent_dir)
            raise RuntimeError(msg)
        # Add these entries to samples_dict
        for d in sample_dirs:
            samples_dict[d[0]] = d[1]

    # Making a dataframe to write out to csv
    samples_df = pd.DataFrame.from_dict(samples_dict, orient='index', columns=['input_dir'])
    samples_df.index.name = "sample"
    samples_df.reset_index(inplace=True)
    # Make sure we don't overwrite the samplesheet
    _no_overwrite_check([fname])
    samples_df.to_csv(fname, index=False)
    print("Done.")


def get_pipe_version(vfn):
    with open(vfn, "r") as vfh:
        # version_info.smk gives access to VER_INFO with software version info
        # and also to ENV_INFO - a singularity info dict (not currently used in init)
        WORKFLOW_BASEDIR = _WATERMELON_ROOT # WORKFLOW_BASEDIR needed in snakefile
        _locals = locals() # capturing locals gives access to WORKFLOW_BASEDIR in snakefile
        exec(vfh.read(), _locals, globals()) # Read in the snakefile
        version = VER_INFO["watermelon"]
        return version


def get_template(args, pipe_root):
    if args.x_alt_template:
        warnings.warn("Using alternate template given. Ignoring type")
        fpath = args.x_alt_template
    else:
        fname = "template_{}.yaml".format(args.type)
        fpath = os.path.join(pipe_root, "config", fname)
    try:
        with open(fpath, "r") as tfh:
            template_dict = ruamel_yaml.round_trip_load(tfh)
            return template_dict
    except:
        msg = f"\n\nCould not read template config {fpath}."
        raise RuntimeError(msg)


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
    if args.sample_sheet:
        config.insert(0, "samplesheet", os.path.abspath(args.sample_sheet))
    # Auto-generate samplesheet if input run dirs specified
    elif args.input_run_dirs:
        ss_fn = "samplesheet.csv"
        _generate_samplesheet(args.input_run_dirs, ss_fn)
        config.insert(0, "samplesheet", os.path.abspath(ss_fn))
    config.insert(0, "watermelon_version", version)
    config.insert(0, "email", email)

    # Add analyst name to report info
    analyst_name = _get_analyst_name(args.x_analyst_info, user)
    if config.get("report_info"):
        config["report_info"]["analyst_name"] = analyst_name

    return config

def validate_fastq_dirs(ssfh, sample_col, fq_col):
    samplesheet = pd.read_csv(ssfh, comment="#", dtype="object")
    if not sample_col in samplesheet.columns:
        msg = "\n\nThe sample sheet must have a column labeled '{}'\n".format(sample_col)
        raise RuntimeError(msg)
    if not fq_col in samplesheet.columns:
        msg = "\n\nThe sample sheet must have a column labeled '{}'\n".format(fq_col)
        raise RuntimeError(msg)

    samples = samplesheet[sample_col]
    fq_dirs = samplesheet[fq_col]
    # fastq_dict = OrderedDict() # If wanted to store this info somewhere
    for sample, dir in zip(samples, fq_dirs):
        if os.path.exists(dir):
            # The helper function has the following validations:
            # * Fastq file(s) must be present
            # * Can't mix gz and fastq files (must be all same type)
            fastqs = helper.get_sample_fastq_paths(dir)
            if not fastqs:
                msg = "\n\nNo fastq files found in sample directory {}.\n"
                raise RuntimeError(msg)
            # else:
            #     fastq_dict[sample] = fastqs
            #     TODO: could easily write this to a file...
        else:
            msg = "\n\nSample directory {} cannot be read.\n".format(dir)
            raise RuntimeError(msg)


def validate_genomes(ref_dict, type):
    # Validate that references are readable
    # Make a copy of the dict (don't actually want to manipulate it)
    ref_dict = copy.copy(ref_dict)
    # Diffex only uses annotation_tsv, not the others
    if type == "diffex":
        for k in ["fasta", "gtf", "rsem_star_index"]:
            ref_dict.pop(k)
    cant_read = []
    for k,v in ref_dict.items():
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
    parser.add_argument("-t", "--type", default="align_qc", choices=["align_qc", "diffex"], help="Type of config file to produce. Uses the appropriate template, and has the appropriate values.")
    samples_in = parser.add_mutually_exclusive_group(required=True)
    samples_in.add_argument("-s", "--sample_sheet", help="CSV file that contains sample IDs, paths containing samples\' fastq files, and optional columns for phenotype information (required for diffex).")
    samples_in.add_argument("-i", "--input_run_dirs", type=str, nargs="*", help="One or more paths to run dirs. Each run dir should contain samples dirs which should in turn contain one or more fastq.gz files. The sample names will be derived from the sample directories")
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
        with open(config_dict["samplesheet"]) as ssfh:
            validate_fastq_dirs(ssfh, args.x_sample_col, args.x_input_col)
    elif args.type == "diffex":
        foo = "bar"

    write_stuff(
        config_dict=config_dict,
        config_fn="config_{}.yaml".format(args.project_id),
        wat_dir=_WATERMELON_ROOT
    )
