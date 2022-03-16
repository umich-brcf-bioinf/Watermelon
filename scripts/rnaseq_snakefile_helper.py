'''Functions that support the RNASeq Snakemake snakefile.
'''

from collections import defaultdict
from functools import partial
import glob
import gzip
import os
import pandas as pd
import re
import subprocess
import warnings

from snakemake import workflow


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

def _filter_dict_by_keys(dict_in, keep_keys):
    all_keys = set(dict_in.keys())
    keep_keys = set(keep_keys)
    if keep_keys - all_keys:
        #watermelon_init validation should ensure that this doesn't happen for samplesheet/input_dir discrepancies, but you never know...
        msg_fmt = 'Some key(s) given but not found in dict:\n{}'
        raise RuntimeError(msg_fmt.format(keep_keys - all_keys))
    if all_keys - keep_keys:
        msg_fmt = 'Filtering the following keys from dict:\n{}'
        warnings.warn(msg_fmt.format(all_keys - keep_keys))
    #Use dict comprehension to create filtered dict: https://stackoverflow.com/a/3420156
    dict_out = {k: dict_in[k] for k in all_keys.intersection(keep_keys)}
    return dict_out


def gunzip_file(source_file, dest_file):
    with gzip.GzipFile(source_file, 'rb') as inF, \
         open(dest_file, 'wb') as outF:
        data = inF.read()
        outF.write(data)
        os.remove(source_file)


def gunzip_glob(source_file_pattern):
    for source_filename in glob.glob(source_file_pattern):
        dest_filename = source_filename.rstrip('.gz')
        gunzip_file(source_filename, dest_filename)


def email(email_config, subject_prefix, msg="", attachment=""):
    if not email_config:
        print(subject_prefix, msg)
    else:
        command = "echo '{msg}' | mutt -s '{subject_prefix}{subject}' {attachment} {to}".format(
                to=email_config['to'],
                subject_prefix=subject_prefix,
                subject=email_config['subject'],
                attachment=attachment,
                msg=msg,
                )
        subprocess.run(command, shell=True)


def phenotypes_with_replicates(phenotype_sample_list):
    with_replicates = []
    for (label, values) in phenotype_sample_list.items():
        max_sample_count = max([len(samples) for (value, samples) in values.items()])
        if max_sample_count > 1:
            with_replicates.append(label)
    return sorted(with_replicates)


def phenotype_sample_list(samplesheet):
    '''Translates samplesheet df into nested dict of phenotypes and their samples.
    Specifically {phenotype_label : {phenotype_value : [list of samples] } }
    '''

    phenotype_dict = defaultdict(partial(defaultdict, list))

    # {phenotype_label : {sample: phenotype_value } }
    samplesheet_dict = samplesheet.to_dict(orient='dict')
    for label in samplesheet_dict:
        for samp in samplesheet_dict[label]:
            value = samplesheet_dict[label][samp]
            if value:
                phenotype_dict[label][value].append(samp)
    # {phenotype_label : {phenotype_value : [list of samples] } }
    return phenotype_dict


def sample_bnames_from_filenames(samplesheet, capture_regex, bname_fmt):
    # Uses the samplesheet to get sample names and fastq directories,
    # then creates a dict with file basenames like {sample}_R{readnum} for each sample in the dict.
    # Returns:
    #     bnames_dict: list - dict with samples as keys and list of basenames as values
    bnames_dict = defaultdict(list)
    for row in samplesheet.itertuples(index=True, name=None):
        sample = row[0]
        fq_list = get_sample_fastq_paths(row[1])
        captured_groups = set()
        for file in fq_list:
            basename = os.path.basename(file)
            rmatch = re.match(capture_regex, basename)
            if not rmatch:
                msg_fmt = 'File {} did not match regular expression {}. Cannot capture grouping information from filename.'
                raise RuntimeError(msg_fmt.format(basename, capture_regex))
            else:
                captured = rmatch.group(1)
                captured_groups.add(captured)
        sample_bnames = [bname_fmt.format(sample, group) for group in captured_groups]
        bnames_dict[sample] = sorted(sample_bnames)
    return bnames_dict


def gather_basenames(sample_bnames_dict, sample_list):
    '''Gather a list of basenames from sample_bnames_dict, for all samples in sample_list'''
    basenames = []
    filtered_dict = _filter_dict_by_keys(sample_bnames_dict, sample_list)
    #Gather basenames from filtered basenames_dict
    for v in filtered_dict.values():
        basenames.extend(v)
    return sorted(basenames)


def get_sample_fastq_paths(sample_fastq_glob):
    #Uses file glob to gather list of all .fastq or .fastq.gz files in a sample's fastq dir.
    #Raises an error if sample has mixed .fastq.gz and .fastq files
    fastqs = glob.glob(sample_fastq_glob)
    ext_gz = set([x.endswith(".gz") for x in fastqs]) # Set from the True/False values in the list
    #Each sample should have only one or the other, plaintext or gz, not mixed
    if len(ext_gz) == 2:
        msg_fmt = "{} contains a mixture of fastq and fastq.gz files. Each sample must have either gzipped or plaintext fastqs, not both."
        raise RuntimeError(msg_fmt.format(sample_fastq_glob))
    elif len(ext_gz) == 1:
        return(sorted(fastqs))
    elif len(ext_gz) == 0:
        return None


def fastqs_to_concat(samplesheet, capture_regex):
    # Uses the samplesheet to get sample names and fastq directories
    # Then gets grouping information for each directory of fastq files,
    # and creates a dict of input files to concatenate based on that grouping
    # Returns:
    #     cat_fq_dict: Dict with sample and group representing top-level and 1st-level keys, repsectively
    cat_fq_dict = dict()
    for row in samplesheet.itertuples(name="samplesheet"):
        grouping_dict = defaultdict(list)
        sample = row.Index
        fq_list = get_sample_fastq_paths(getattr(row, "input_glob"))
        for file in fq_list:
            basename = os.path.basename(file)
            rmatch = re.match(capture_regex, basename)
            if not rmatch:
                msg_fmt = "File {} did not match regular expression {}. Cannot capture grouping information from filename."
                raise RuntimeError(msg_fmt.format(basename, capture_regex))
            else:
                captured = rmatch.group(1)
                grouping_dict[captured].append(file)
        cat_fq_dict[sample] = grouping_dict
    return cat_fq_dict


def detect_paired_end_bool(fastqs):
    if len(fastqs) == 2:
        return(True)
    elif len(fastqs) == 1:
        return(False)
    else:
        msg_format = 'Found {} fastqs ({}); expected either 1 or 2 fastq files/sample'
        raise ValueError(msg_format.format(len(fastqs), ','.join(fastqs)))


def diffex_models(diffex_config):
    not_factors = ['count_min_cutoff']
    model_names = [k for k in diffex_config.keys() if k not in not_factors]
    return(model_names)


def diffex_model_info(diffex_config):
    '''Returns two dictionaries: one with contrasts, and one (nested) with linear fold-change and p-val cutoff, for each model.
    Keyed by model name. Former is used to setup targets in snakefile, both used in report.
    Examples:
    cont_dict = { 'model_foo' : ['val1_v_val2', 'val3_v_val4'] }
    info_dict = { 'model_foo' : {
        'linear_fold_change' : 1.5,
        'adjustedPValue' : 0.05
    }'''
    info_dict = {}
    cont_dict = {}
    for model in diffex_models(diffex_config):
        info_dict[model] = {
            'model' : diffex_config[model]['DESeq2']['design'],
            'factor_name' : diffex_config[model]['DESeq2']['factor_name'],
            'linear_fold_change' : diffex_config[model]['linear_fold_change'],
            'adjustedPValue' : diffex_config[model]['adjustedPValue']
        }
        cont_dict[model] = diffex_config[model]['contrasts']
    return(info_dict, cont_dict)


def expand_model_contrast_filenames(model_contrasts_format, contrast_dict):
    paths = []
    for model in contrast_dict:
        conts = contrast_dict[model]
        paths.extend(workflow.expand(model_contrasts_format, model_name=model, contrast=conts))
    return(paths)
