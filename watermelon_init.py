#!/usr/bin/env python3
'''
Produces a config.yaml to use with snakemake
'''

import argparse
import copy
import getpass
import os
import re
import ruamel_yaml
import subprocess
import sys
import warnings

import pandas as pd
from collections import OrderedDict
from collections.abc import Mapping
from natsort import natsorted

import scripts.rnaseq_snakefile_helper as helper # works on CL and in pytest

import pdb # TWS DEBUG


# Set default values
# Find path knowing current location of Watermelon/watermelon_init
_WATERMELON_ROOT = os.path.realpath(os.path.dirname(__file__))
_DEFAULT_GENOME_REFERENCES = os.path.join(_WATERMELON_ROOT, "config", "genome_references.yaml")
_DEFAULT_ACKNOWLEDGEMENT = os.path.join(_WATERMELON_ROOT, "report", "default_acknowledge.txt")
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
    type_dirs = {
        'align_qc': ['alignment_results', 'deliverables', 'report'],
        'diffex': ['diffex_results', 'deliverables', 'report']
        }
    dirs = ruamel_yaml.comments.CommentedMap() # Use this instead of an OrderedDict since it prints pretty (no !!omap)
    for resulttype in type_dirs.get(type, []):
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
    # TestData doesn't exist yet. It will when Watermelon is copied to cwd
    if gbuild == "TestData":
        datapath = os.path.join(os.getcwd(), "Watermelon", "data")
        our_refs = {
            'genome' : 'GRCh38',
            'references' : {
                'fasta' : os.path.join(datapath, 'Homo_sapiens.GRCh38.dna_sm.chr22.fa'),
                'gtf' : os.path.join(datapath, 'Homo_sapiens.GRCh38.98.chr22.gtf'),
                'annotation_tsv' : os.path.join(datapath, 'Homo_sapiens.GRCh38.98_annotation.tsv'),
                'coverage_bed' : os.path.join(datapath, 'Homo_sapiens.GRCh38.98.chr22.bed')
            }
        }
    else:
        try:
            our_refs = refs[gbuild]
        except KeyError:
            msg_fmt = "\n\nGenome build {} not found in {}:\n {}\n"
            raise RuntimeError(msg_fmt.format(gbuild, grefsfn, refs.keys()))
    # If diffex type, only need annotation TSV. If align_qc, keep it all
    if type == "diffex":
        for ref in ["fasta", "gtf", "rsem_star_index"]:
            our_refs["references"].pop(ref, None)
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


def generate_samplesheet(fq_parent_dirs):
    samples_dict = OrderedDict()
    for parent_dir in fq_parent_dirs:
        if not os.path.exists(parent_dir):
            msg = "\n\nInput run directory {} does not exist.\n".format(parent_dir)
            raise RuntimeError(msg)
        # list of tuples
        sample_dirs = [(d.name, os.path.abspath(d)) for d in os.scandir(parent_dir) if d.is_dir()]
        sample_dirs = natsorted(sample_dirs)
        if not sample_dirs:
            msg = "\n\nInput run directory {} contains no subdirectories.\n".format(parent_dir)
            raise RuntimeError(msg)
        # Add these entries to samples_dict
        for d in sample_dirs:
            samples_dict[d[0]] = d[1]
    # Making a dataframe to write out to csv
    samples_df = pd.DataFrame.from_dict(samples_dict, orient='index', columns=['input_dir'])
    samples_df.index.name = "sample"
    samples_df.reset_index(inplace=True)
    return samples_df


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

    # Set up some things before adding them
    refs = _set_up_refs(args.x_genome_references, args.genome_build, args.type)
    dirs = _set_up_dirs(args.type, args.project_id)
    user = getpass.getuser() # used for email and for analyst name below
    email = _set_up_email(args.project_id, user)

    # Merge the genome reference dict with the config
    # Inserts top-level keys "genome" & "references" into config
    # Integrates "fastq_screen" ref vals with template vals
    _dict_merge(config, refs)
    # If "Other" or "TestData" genome_build is specified, or if diffex type config,
    # remove fastq_screen section from template config
    if args.genome_build == "Other" or args.genome_build == "TestData" or args.type == "diffex":
        config.pop("fastq_screen", None)
    # Now add things for the desired order:
    # email, watermelon_version, samplesheet, genome, references, dirs, template_config stuff
    # Focus on the 2nd half (working backwards)
    config.insert(0, "dirs", dirs) # Insert dirs on top
    config.move_to_end("references", last=False) # This is how to move existing keys to the top
    config.move_to_end("genome", last=False)
    # 1st half (working backwards)
    # If non-default sample col is given, add it to the config
    if args.x_sample_col != _DEFAULT_SAMPLE_COL:
        config.insert(0, "sample_col", args.x_sample_col)
    # Set samplesheet absolute path based on how inputs were given
    # Also get sequencing info from a file accordingly
    if args.input_run_dirs:
        # default 'samplesheet.csv' if it will be auto-generated later
        ssfp = os.path.join(os.getcwd(), "samplesheet.csv")
    elif args.sample_sheet:
        ssfp = os.path.abspath(args.sample_sheet)
    if args.type == "diffex":
        if not args.count_matrix:
            msg = "\n\nCannot create diffex config without count_matrix argument.\n"
            raise RuntimeError(msg)
        config.insert(0, "count_matrix", args.count_matrix)
    config.insert(0, "samplesheet", ssfp)
    config.insert(0, "watermelon_version", version)
    config.insert(0, "email", email)

    # Add analyst name to report info
    # Update report project_name with args.project_id
    analyst_name = _get_analyst_name(args.x_analyst_info, user)
    if config.get("report_info"):
        config["report_info"]["analyst_name"] = analyst_name
        config["report_info"]["project_name"] = args.project_id
        # Insert sequencing info from file
        with open(args.background_info, "r") as fh:
            lines = fh.readlines()
            linestring = "".join([l.strip() for l in lines])
            config["report_info"]["prep_description"] = linestring
        # Insert acknowledgement text from file
        with open(args.x_acknowledgement_text, "r") as fh:
            lines = fh.readlines()
            linestring = "".join(lines) # Preserve newlines to separate text
            if args.AGC:
                linesting = re.sub("Bioinformatics", "Advanced Genomics", linestring)
            config["report_info"]["acknowledgement_text"] = linestring
        # Set (other) AGC-specific report details
        if args.AGC:
            config["report_info"]["include_follow_up"] = False
            config["report_info"]["include_pct_dups"] = False

    return config


def validate_count_matrix(counts_fp):
    counts = pd.read_csv(counts_fp, sep="\t").set_index("gene_id", drop=True)
    if not all(counts.dtypes == "int64"):
        msg = "\n\nThe count matrix must only contain columns of numeric count data\n"
        raise RuntimeError(msg)


def validate_fastq_dirs(ss_df, sample_col, fq_col):
    if not sample_col in ss_df.columns:
        msg = "\n\nThe sample sheet must have a column labeled '{}'\n".format(sample_col)
        raise RuntimeError(msg)
    if not fq_col in ss_df.columns:
        msg = "\n\nThe sample sheet must have a column labeled '{}'\n".format(fq_col)
        raise RuntimeError(msg)

    samples = ss_df[sample_col]
    fq_dirs = ss_df[fq_col]
    # fastq_dict = OrderedDict() # If wanted to store this info somewhere
    for sample, dir in zip(samples, fq_dirs):
        if os.path.exists(dir):
            # The helper function assists/provides the following validations:
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


def validate_genomes(ref_dict):
    # Validate that references are readable
    # Make a copy of the dict (don't actually want to manipulate it)
    ref_dict = copy.copy(ref_dict)
    ref_dict.pop("Other", None) # This is a placeholder and doesn't need validation
    cant_read = []
    for k,v in ref_dict.items():
        if k == "rsem_star_index":
            v = v + ".chrlist" # Given index prefix, test for chromosome list
        if not os.access(v, os.R_OK): # Test for read access
            cant_read.append(v)
    if cant_read:
        msg = "\n\nCan't read references:\n{}\n"
        raise RuntimeError(msg)


def write_stuff(config_dict, config_fn, wat_dir, ss_df=None):
    # Make sure we don't overwrite things
    check_list = ["Watermelon", config_fn]
    if isinstance(ss_df, pd.DataFrame):
        # Add samplesheet to overwrite check list if auto-generated
        check_list.append(os.path.relpath(config_dict["samplesheet"]))
    _no_overwrite_check(check_list)
    print("Writing config to working dir...")
    # Write the config file
    with open(config_fn, "w") as config_fh:
        ruamel_yaml.round_trip_dump(
            config_dict, config_fh, default_flow_style=False, indent=4
        )
    # Write auto-generated samplesheet if present
    if isinstance(ss_df, pd.DataFrame):
        print("Done.\nWriting auto-generated samplesheet...")
        ss_df.to_csv(config_dict["samplesheet"], index=False)
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
    parser.add_argument("-b", "--background_info", type=str, default=os.path.join(_WATERMELON_ROOT, "report", "default_seq_info.txt"), help="A text file containing background information (e.g. sequencing and/or preparation procedure). This text will be used in the report. Defaults to 'Watermelon/report/default_seq_info.txt'")
    parser.add_argument("-c", "--count_matrix", type=str, help="Count data used as input for diffex pipeline. Contains a column of gene IDs and additional columns of count data from aligned reads.")
    parser.add_argument("-g", "--genome_build", type=str, required=True, help="Reference files for this build are added into the config.")
    parser.add_argument("-o", "--output_prefix", type=str, default="analysis_{project_id}", help="Optional output prefix. Defaults to relative path 'analysis_{project_id}'")
    parser.add_argument("-p", "--project_id", type=str, required=True, help="ID used for the project. Can represent project ID, service request ID, etc. Only alphanumeric and underscores allowed. No spaces or special characters.")
    parser.add_argument("-t", "--type", type=str, default="align_qc", choices=["align_qc", "diffex"], help="Type of config file to produce. Uses the appropriate template, and has the appropriate values.")
    samples_in = parser.add_mutually_exclusive_group(required=True)
    samples_in.add_argument("-s", "--sample_sheet", type=str, help="CSV file that contains sample IDs, paths containing samples\' fastq files, and optional columns for phenotype information (required for diffex).")
    samples_in.add_argument("-i", "--input_run_dirs", type=str, nargs="*", help="One or more paths to run dirs. Each run dir should contain samples dirs which should in turn contain one or more fastq.gz files. The sample names will be derived from the sample directories.")
    parser.add_argument("--AGC", default=False, action="store_true", help="An optional flag which enables configuration options particularly needed by the UMich Advanced Genomics Core.")
    parser.add_argument("--x_analyst_info", type=str, default=_DEFAULT_ANALYST_INFO, help=argparse.SUPPRESS)
    parser.add_argument("--x_alt_template", type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--x_genome_references", type=str, default=_DEFAULT_GENOME_REFERENCES, help=argparse.SUPPRESS)
    parser.add_argument("--x_acknowledgement_text", type=str, default=_DEFAULT_ACKNOWLEDGEMENT, help=argparse.SUPPRESS)
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

    # Auto-generate samplesheet if input run dirs specified
    if args.input_run_dirs:
        samplesheet_df = generate_samplesheet(args.input_run_dirs)
    else:
        # Read in the given samplesheet for validation
        samplesheet_df = pd.read_csv(config_dict["samplesheet"], comment="#", dtype="object")

    # Perform config validations
    if args.genome_build != "Other" and args.genome_build != "TestData":
        validate_genomes(config_dict["references"])
    # Some validations depend on type of config
    if args.type == "align_qc":
        validate_fastq_dirs(samplesheet_df, args.x_sample_col, args.x_input_col)
    elif args.type == "diffex":
        validate_count_matrix(config_dict["count_matrix"])

    # Write the outputs
    write_kwargs = { # Using kwargs for optional samplesheet_df if it's to be written
        "config_dict": config_dict,
        "config_fn": "config_{}.yaml".format(args.project_id),
        "wat_dir": _WATERMELON_ROOT
    }
    if args.input_run_dirs:
        write_kwargs["ss_df"] = samplesheet_df
    write_stuff(**write_kwargs)

    if args.genome_build == "TestData":
        # These need to be unzipped after watermelon is copied
        data_gz_glob = os.path.join(os.getcwd(), "Watermelon", "data", "*.gz")
        helper.gunzip_glob(data_gz_glob)
