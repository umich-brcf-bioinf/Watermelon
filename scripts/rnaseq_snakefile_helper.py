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

from scripts import watermelon_config
from scripts.watermelon_config import CONFIG_KEYS
from scripts.watermelon_config import DEFAULT_COMPARISON_INFIX
from scripts.watermelon_config import DEFAULT_PHENOTYPE_DELIM

HISAT2_NAME = 'HISAT2'
STRINGTIE_NAME = 'stringtie'
RSEM_NAME = 'rsem'

STRAND_CONFIG_PARAM = { 'fr-unstranded' : {HISAT2_NAME: '', STRINGTIE_NAME: '', RSEM_NAME: ''},
                        'fr-firststrand' : {HISAT2_NAME : '--rna-strandness RF', STRINGTIE_NAME: '--rf', RSEM_NAME: '--strandedness reverse'},
                        'fr-secondstrand' : {HISAT2_NAME : '--rna-strandness FR', STRINGTIE_NAME: '--fr', RSEM_NAME: '--strandedness forward'},
                        'unstranded' : {HISAT2_NAME: '', STRINGTIE_NAME: '', RSEM_NAME: ''},
                        'reverse_forward' : {HISAT2_NAME : '--rna-strandness RF', STRINGTIE_NAME: '--rf', RSEM_NAME: '--strandedness reverse'},
                        'forward_reverse' : {HISAT2_NAME : '--rna-strandness FR', STRINGTIE_NAME: '--fr', RSEM_NAME: '--strandedness forward'},
                        }

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
                 config={},
                 delimiter=DEFAULT_PHENOTYPE_DELIM,
                 comparison_infix=DEFAULT_COMPARISON_INFIX):
        #self.phenotype_labels_string = config.get(CONFIG_KEYS.phenotypes, None)
        self.samplesheet = pd.read_csv(config["sample_description_file"]).set_index("sample", drop=True)
        #self.phenotype_labels_string = self.samplesheet.columns
        self.sample_phenotype_value_dict = self.samplesheet.to_dict(orient='index')
        self.comparisons = config.get(CONFIG_KEYS.comparisons, None)
        self.delimiter = delimiter
        self.comparison_infix = comparison_infix

    @property
    def phenotype_sample_list(self):
        '''Translates config phenotypes/samples into nested dict of phenotypes.
        Specifically {phenotype_label : {phenotype_value : [list of samples] } }

        Strips all surrounding white space.

        phenotypes_string : delimited phenotype labels (columns)
        sample_phenotype_value_dict : {sample_id : delimited phenotype_value_string} (rows)
        '''
        def check_labels_match_values(phenotype_labels,
                                      sample,
                                      phenotype_values):
            if len(phenotype_labels) != len(phenotype_values):
                msg_fmt = 'expected {} phenotype values but sample {} had {}'
                msg = msg_fmt.format(len(phenotype_labels),
                                     sample,
                                     len(phenotype_values))
                raise ValueError(msg)

        def check_phenotype_labels(labels):
            for i,label in enumerate(labels):
                if label == '':
                    msg_fmt = 'label of phenotype {} is empty'
                    msg = msg_fmt.format(i + 1)
                    raise ValueError(msg)

        phenotype_dict = defaultdict(partial(defaultdict, list))
        #phenotype_labels = list(map(str.strip, self.phenotype_labels_string.split(self.delimiter)))
        phenotype_labels = list(self.samplesheet.columns)
        check_phenotype_labels(phenotype_labels)
        sorted_sample_phenotypes_items = sorted([(k, v) for k,v in self.sample_phenotype_value_dict.items()])
        for sample, phenotype_values in sorted_sample_phenotypes_items:
            #sample_phenotype_values = list(map(str.strip, phenotype_values.split(self.delimiter)))
            sample_phenotypes = dict(zip(phenotype_labels, phenotype_values))
            check_labels_match_values(phenotype_labels,
                                      sample,
                                      phenotype_values)
            for label, value in sample_phenotypes.items():
                if value != '':
                    phenotype_dict[label][value].append(sample)
        return phenotype_dict

    @property
    def phenotypes_with_replicates(self):
        with_replicates = []
        for (label, values) in self.phenotype_sample_list.items():
            max_sample_count = max([len(samples) for (value, samples) in values.items()])
            if max_sample_count > 1:
                with_replicates.append(label)
        return sorted(with_replicates)

    def separated_comparisons(self, output_delimiter):
        '''Returns a dict of {phenotype_label: 'val1_v_val2,val3_v_val4'}'''
        comparisons = {}
        for phenotype, comparison in self.comparisons.items():
            comp_string = output_delimiter.join(list(map(str.strip, comparison)))
            comparisons[phenotype]= comp_string
        return comparisons

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

    def concatenated_comparison_values(self, output_delimiter):
        '''Returns dict of {pheno_label : 'val1,val2,val3,val4' }'''
        phenotype_label_values = {}
        for phenotype_label, values in self.comparison_values.items():
            phenotype_label_values[phenotype_label] = output_delimiter.join(values)
        return phenotype_label_values

    @property
    def phenotypes_comparisons_all_tuple(self):
        '''Returns named tuple of phenotype, comparison for all phenotypes'''
        return self._phenotypes_comparisons(True)

    @property
    def phenotypes_comparisons_replicates_tuple(self):
        '''Returns named tuple of phenotype, comparison for phenotypes with replicates'''
        return self._phenotypes_comparisons(False)

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

def _strand_option(program_name, strand_option):
    if not strand_option in STRAND_CONFIG_PARAM:
        msg_format = ('ERROR: config:alignment_options:library_type={} is '
                      'not valid. Valid library_type options are: {}')
        msg = msg_format.format(strand_option, ','.join(STRAND_CONFIG_PARAM.keys()))
        raise ValueError(msg)
    return STRAND_CONFIG_PARAM[strand_option][program_name]

def strand_option_hisat2(config):
    config_strand_option = config["alignment_options"]["library_type"]
    return _strand_option(HISAT2_NAME, config_strand_option)

def strand_option_stringtie(config):
    config_strand_option = config["alignment_options"]["library_type"]
    return _strand_option(STRINGTIE_NAME, config_strand_option)

def strand_option_rsem(config):
    config_strand_option = config["alignment_options"]["library_type"]
    return _strand_option(RSEM_NAME, config_strand_option)

def hisat_detect_paired_end(fastqs):
    if len(fastqs) == 2:
        return('-1 ' + fastqs[0] + ' -2 ' + fastqs[1])
    elif len(fastqs) == 1:
        return('-U ' + fastqs[0])
    else:
        msg_format = 'Found {} fastqs ({}); expected either 1 or 2 fastq files/sample'
        raise ValueError(msg_format.format(len(fastqs), ','.join(fastqs)))

def detect_paired_end_bool(fastqs):
    if len(fastqs) == 2:
        return(True)
    elif len(fastqs) == 1:
        return(False)
    else:
        msg_format = 'Found {} fastqs ({}); expected either 1 or 2 fastq files/sample'
        raise ValueError(msg_format.format(len(fastqs), ','.join(fastqs)))

def cutadapt_options(trim_params):
    run_trimming_options = 0
    for option, value in trim_params.items():
        if not isinstance(value, int):
            msg_format = "ERROR: config:trimming_options '{}={}' must be integer"
            msg = msg_format.format(option, value)
            raise ValueError(msg)
        if value != 0:
            run_trimming_options = 1
    return run_trimming_options

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

def transform_config(config):
    watermelon_config.transform_config(config)

def diffex_models(diffex_config):
    not_models = ['adjustedPValue', 'fold_change']
    model_names = [k for k in diffex_config.keys() if k not in not_models]
    return(model_names)

def expand_DESeq2_model_contrasts(model_contrasts_format, diffex_config):
    paths = []
    for model in diffex_models(diffex_config):
        if 'DESeq2' in diffex_config[model]:
            contrasts = diffex_config[model]['DESeq2']['results']['contrasts']
            #Create contrast info for filenames
            #Split on ^, throw away first item, then join remaining two with _v_
            cont_fn = map(lambda x: "_v_".join(x.split("^")[1:]), contrasts)
            paths.extend(workflow.expand(model_contrasts_format, model_name=model, contrast=cont_fn))
    return(paths)
