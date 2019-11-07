'''Functions that support the RNASeq Snakemake snakefile.
'''
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import collections
from collections import defaultdict
import csv
from functools import partial
import glob
import hashlib
import os
from os.path import join
from os.path import isfile
import pandas as pd
import sys
import yaml

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

def init_references(config_references):
    def existing_link_target_is_different(link_name, link_path):
        original_abs_path = os.path.realpath(link_name)
        new_abs_path = os.path.realpath(link_path)
        return original_abs_path != new_abs_path
    if not os.path.exists("references"):
        os.mkdir("references")
    os.chdir("references")
    references = config_references if config_references else {}
    for link_name, link_path in references.items():
        if not os.path.exists(link_path):
            msg_fmt = 'ERROR: specified config reference files/dirs [{}:{}] cannot be read'
            msg = msg_fmt.format(link_name, link_path)
            raise ValueError(msg)
        elif not os.path.exists(link_name):
            os.symlink(link_path, link_name)
        elif existing_link_target_is_different(link_name, link_path):
            os.remove(link_name)
            os.symlink(link_path, link_name)
        else:
            pass #link matches existing link
    os.chdir("..")

class PhenotypeManager(object):
    '''Interprets a subset of the config to help answer questions around how
    samples map to phenotype labels and values and vice versa.'''
    def __init__(self,
                 config={}):
        self.samplesheet = pd.read_csv(config["sample_description_file"], comment='#', keep_default_na=False).set_index("sample", drop=True)
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

    @property
    def comparison_values(self):
        '''Returns dict of {pheno_label : [val1,val2,val3,val4] }'''
        phenotype_label_values = {}
        for phenotype, comparison_list in self.comparisons.items():
            comparison_list = list(map(str.strip, comparison_list))
            unique_conditions = set()
            for group in comparison_list:
                unique_conditions.update(group.split(self.comparison_infix))
            values = sorted(unique_conditions)
            phenotype_label_values[phenotype] = values
        return phenotype_label_values

    def _phenotypes_comparisons(self, include_phenotypes_without_replicates=True):
        if include_phenotypes_without_replicates:
            include = lambda p: True
        else:
            with_replicates = set(self.phenotypes_with_replicates)
            include = lambda p: p in with_replicates
        pheno_comps = [(p, c) for p,c in self.comparisons.items() if include(p)]
        phenotypes=[]
        comparisons=[]
        for pheno, comps in sorted(pheno_comps):
            for comp in sorted(comps):
                phenotypes.append(pheno)
                comparisons.append(comp.strip())
        PhenotypesComparisons = collections.namedtuple('PhenotypeComparisons',
                                                       'phenotypes comparisons')
        return PhenotypesComparisons(phenotypes=phenotypes, comparisons=comparisons)

def detect_paired_end_bool(fastqs):
    if len(fastqs) == 2:
        return(True)
    elif len(fastqs) == 1:
        return(False)
    else:
        msg_format = 'Found {} fastqs ({}); expected either 1 or 2 fastq files/sample'
        raise ValueError(msg_format.format(len(fastqs), ','.join(fastqs)))

def _get_sample_reads(fastq_base_dir, samples):
    sample_reads = {}
    read_suffix = {0: "", 1:"_SE", 2:"_PE"}
    is_fastq = lambda fn: fn.endswith('.fastq') or fn.endswith('.fastq.gz')
    for sample in samples:
        found_reads = []
        for read in ['R1', 'R2']:
            sample_dir = join(fastq_base_dir, sample, '')
            read_present = list(filter(is_fastq, glob.glob(sample_dir + '*_{}*'.format(read))))
            if read_present:
                found_reads.append(read)
        found_reads = list(map(lambda x: x+read_suffix[len(found_reads)], found_reads))
        sample_reads[sample] = found_reads
    return sample_reads

def flattened_sample_reads(fastq_base_dir, samples):
    sample_reads = _get_sample_reads(fastq_base_dir, samples)
    return sorted([(sample,read) for (sample,reads) in sample_reads.items() for read in reads])

def expand_sample_read_endedness(sample_read_endedness_format,
                                 all_flattened_sample_reads,
                                 sample=None):
    if not all_flattened_sample_reads:
        return []
    samples, reads = zip(*[(s, r) for s, r in all_flattened_sample_reads if not sample or s == sample])
    return workflow.expand(sample_read_endedness_format,
                           zip,
                           sample=samples,
                           read_endedness=reads)

def expand_read_stats_if_paired(read_stats_filename_format,
                                flattened_sample_reads,
                                sample):
    result = []
    if len([s for s,r in flattened_sample_reads if s == sample]) > 1:
        result.append(read_stats_filename_format.format(sample=sample))
    return result

def diffex_models(diffex_config):
    not_factors = ['adjustedPValue', 'linear_fold_change']
    model_names = [k for k in diffex_config.keys() if k not in not_factors]
    return(model_names)

def diffex_contrasts(diffex_config):
    '''Returns contrasts dict in the form of
    {model_name: ['val1_v_val2', 'val3_v_val4']}
    '''
    cont_dict = {}
    for model in diffex_models(diffex_config):
        cont_dict[model] = diffex_config[model]['contrasts']
    return(cont_dict)

def expand_model_contrast_filenames(model_contrasts_format, contrast_dict):
    paths = []
    for model in contrast_dict:
        conts = contrast_dict[model]
        paths.extend(workflow.expand(model_contrasts_format, model_name=model, contrast=conts))
    return(paths)
