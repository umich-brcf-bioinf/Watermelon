'''Functions that support the RNASeq Snakemake snakefile.
'''

from collections import defaultdict
from functools import partial
import glob
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

class PhenotypeManager(object):
    '''Interprets a subset of the config to help answer questions around how
    samples map to phenotype labels and values and vice versa.'''
    def __init__(self,
                 config={}):
        self.samplesheet = pd.read_csv(config["samplesheet"], comment='#', dtype='string', keep_default_na=False).set_index("sample", drop=True)
        #sample_phenotype_value_dict : {sample : { pheno_label: pheno_value } }
        self.sample_phenotype_value_dict = self.samplesheet.to_dict(orient='index')

    @property
    def phenotype_sample_list(self):
        '''Translates config phenotypes/samples into nested dict of phenotypes.
        Specifically {phenotype_label : {phenotype_value : [list of samples] } }
        '''

        phenotype_dict = defaultdict(partial(defaultdict, list))

        # {phenotype_label : {sample: phenotype_value } }
        samplesheet_dict = self.samplesheet.to_dict(orient='dict')
        for label in samplesheet_dict:
            for samp in samplesheet_dict[label]:
                value = samplesheet_dict[label][samp]
                if value:
                    phenotype_dict[label][value].append(samp)
        # {phenotype_label : {phenotype_value : [list of samples] } }
        return phenotype_dict

    @property
    def phenotypes_with_replicates(self):
        with_replicates = []
        for (label, values) in self.phenotype_sample_list.items():
            max_sample_count = max([len(samples) for (value, samples) in values.items()])
            if max_sample_count > 1:
                with_replicates.append(label)
        return sorted(with_replicates)

class InputFastqManager(object):
    def __init__(self, input_dir, capture_regex):
        self.input_dir = input_dir
        self.input_paths_dict = self._get_fastq_paths_dict_from_input_dir()
        self.fastqs_to_concat_dict = self._fastqs_to_concat_from_filenames(capture_regex=capture_regex)
        self.sample_bnames_dict = self._sample_bnames_from_filenames(capture_regex=capture_regex, bname_fmt='{}_R{}')


    def _get_sample_fastq_paths(self, sample_fastq_dir):
        #Uses file glob to gather list of all .fastq or .fastq.gz files in a sample's fastq dir.
        #Raises an error if sample has mixed gzipped & plaintext fastqs
        fastq_glob = os.path.join(sample_fastq_dir, '*.fastq')
        fastq_gz_glob = os.path.join(sample_fastq_dir, '*.fastq.gz')
        fastqs = glob.glob(fastq_glob)
        fastq_gzs = glob.glob(fastq_gz_glob)
        #Each sample should have only one or the other, plaintext or gz, not mixed
        if fastqs and fastq_gzs:
            msg_fmt = "{} contains a mixture of fastq and fastq.gz files. Each sample must have either gzipped or plaintext fastqs, not both."
            raise RuntimeError(msg_fmt.format(sample_fastq_dir))
        elif fastqs:
            return(fastqs)
        elif fastq_gzs:
            return(fastq_gzs)
        else:
            return None

    def _get_fastq_paths_dict_from_input_dir(self):
        #Creates a dict of fastq paths {'sample1': ['path/to/sample1_R1.fastq', 'path/to/sample1_R2.fastq']} from samples in input_dir
        #By calling _get_sample_fastq_paths for each subdir of self.input_dir
        input_dir = self.input_dir
        fastq_paths_dict = {}
        #Get list of subdirs
        potential_sample_dirs = next(os.walk(input_dir))[1] #https://stackoverflow.com/a/25705093
        for d in potential_sample_dirs:
            fqs = self._get_sample_fastq_paths(os.path.join(input_dir, d))
            if fqs:
                fastq_paths_dict[d] = sorted(fqs)
            else:
                msg = 'Sample directory {} did not contain .fastq or .fastq.gz files. Skipping'.format(
                    os.path.join(input_dir, d)
                )
                warnings.warn(msg)
        if fastq_paths_dict == {}:
            msg = 'No .fastq or .fastq.gz files in any subdirectories {} of input directory {}. Exiting'.format(
                potential_sample_dirs, input_dir
            )
            raise RuntimeError(msg)
        else:
            return fastq_paths_dict

    def _fastqs_to_concat_from_filenames(self, capture_regex):
        # Uses filenames in input_paths_dict, to get grouping information for each sample's input files,
        # and gathers a list of input files to concatenate based on that grouping for each sample in the dict
        # Arguments:
        #     capture_regex: string - regular expression used to capture group information
        # Returns:
        #     cat_fq_dict: Nested dict with sample and group representing top-level and 1st-level keys, repsectively
        input_files = self.input_paths_dict

        def inner_work(sample_input_files):
            sample_cat_dict = defaultdict(list)
            for file in sample_input_files:
                basename = os.path.basename(file)
                rmatch = re.match(capture_regex, basename)
                if not rmatch:
                    msg_fmt = 'File {} did not match regular expression {}. Cannot capture grouping information from filename.'
                    raise RuntimeError(msg_fmt.format(basename, capture_regex))
                else:
                    captured = rmatch.group(1)
                    sample_cat_dict[captured].append(file)
            return sample_cat_dict

        cat_fq_dict = {key: inner_work(self.input_paths_dict[key]) for key in self.input_paths_dict}
        return cat_fq_dict

    def _sample_bnames_from_filenames(self, capture_regex, bname_fmt):
        # Uses filenames in input_paths_dict, to get grouping information for a sample's input files,
        # Then creates a list of file basenames like {sample}_R{readnum} for each sample in the dict.
        # Arguments:
        #     capture_regex: string - regular expression used to capture group information
        #     bname_fmt: string - format with two groups to be filed in e.g. {}_R{}
        # Returns:
        #     bnames_dict: list - dict with samples as keys and list of basenames as values

        def inner_work(sample, sample_input_files):
            captured_groups = set()
            for file in sample_input_files:
                basename = os.path.basename(file)
                rmatch = re.match(capture_regex, basename)
                if not rmatch:
                    msg_fmt = 'File {} did not match regular expression {}. Cannot capture grouping information from filename.'
                    raise RuntimeError(msg_fmt.format(basename, capture_regex))
                else:
                    captured = rmatch.group(1)
                    captured_groups.add(captured)
            sample_bnames = [bname_fmt.format(sample, group) for group in captured_groups]
            return sorted(sample_bnames)

        bnames_dict = {key: inner_work(key, self.input_paths_dict[key]) for key in self.input_paths_dict}
        return bnames_dict

    def gather_basenames(self, sample_list):
        '''Gather a list of basenames from sample_bnames_dict, for all samples in sample_list'''
        basenames = []
        filtered_dict = _filter_dict_by_keys(self.sample_bnames_dict, sample_list)
        #Gather basenames from filtered basenames_dict
        for v in filtered_dict.values():
            basenames.extend(v)
        return sorted(basenames)


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
