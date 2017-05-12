'''Functions that support DESeq2.
'''
from __future__ import print_function, absolute_import, division
from collections import defaultdict
import os
import sys

import yaml

from scripts.watermelon_config import CONFIG_KEYS, DEFAULT_PHENOTYPE_DELIM, MAIN_FACTOR_TRUE, DEFAULT_COMPARISON_INFIX
from scripts.watermelon_config import split_config_list
from scripts.rnaseq_snakefile_helper import PhenotypeManager

_COMBINATORIC_GROUP = 'combinatoric_group'
_CONTRASTS_HEADER = ['factor','test_level','reference_level', 'base_file_name']
_SAMPLE_METADATA_HEADER = 'sample_name'
_REPLICATE_SUFFIX = '^replicate'

def _build_combinatoric_group_value(phenotype_labels,
                                    main_factors,
                                    replicated_phenotype_labels,
                                    phenotype_values):
    phenotype_label_value = dict([(label, value) for label, value in zip(phenotype_labels, phenotype_values)])
    main_factors = [label for label, main_factor in zip(phenotype_labels, main_factors) if main_factor == MAIN_FACTOR_TRUE]
    replicated_main_factors = [label for label in main_factors if label in replicated_phenotype_labels]
    combinatoric_phenotype_values = []
    for phenotype_label in replicated_main_factors:
        combinatoric_phenotype_values.append(phenotype_label_value[phenotype_label])
    combinatoric_string = DEFAULT_PHENOTYPE_DELIM.join(combinatoric_phenotype_values)
    return combinatoric_string



def _build_sample_metadata_list(config, pheno_with_replicates):
    replicate_count = defaultdict(int)
    def _replicate_number(label, value):
        if value == '':
            return ''
        label_value = (label, value)
        replicate_count[label_value] += 1
        return str(replicate_count[label_value])

    lines = []
    main_factors = split_config_list(config[CONFIG_KEYS.main_factors])
    phenotype_labels = split_config_list(config[CONFIG_KEYS.phenotypes])
    replicated_phenotype_labels = [p for p in phenotype_labels if p in pheno_with_replicates]

    header = [_SAMPLE_METADATA_HEADER]
    for label in replicated_phenotype_labels:
        header.append(label)
        header.append(label + _REPLICATE_SUFFIX)
    header.append(_COMBINATORIC_GROUP)
    header.append(_COMBINATORIC_GROUP + _REPLICATE_SUFFIX)
    
    lines.append(header)
    for sample_name, phenotype_values_string in sorted(config[CONFIG_KEYS.samples].items()):
        sample_line = [sample_name]
        phenotype_values = split_config_list(phenotype_values_string)
        for label, value in zip(phenotype_labels, phenotype_values):
            if label in replicated_phenotype_labels:
                sample_line.append(value)
                sample_line.append(_replicate_number(label, value))
        combinatoric_group_value = _build_combinatoric_group_value(phenotype_labels,
                                                                   main_factors,
                                                                   replicated_phenotype_labels,
                                                                   phenotype_values)
        sample_line.append(combinatoric_group_value)
        sample_line.append(_replicate_number(_COMBINATORIC_GROUP, combinatoric_group_value))
        lines.append(sample_line)
    return lines

def _build_contrasts_list(config, phenos_with_replicates):
    lines = []
    lines.append(_CONTRASTS_HEADER)
    for pheno_label, comparison_strings in sorted(config[CONFIG_KEYS.comparisons].items()):
        comparisons = sorted([x.split(DEFAULT_COMPARISON_INFIX) for x in comparison_strings])
        for test_level, reference_level in comparisons:
            comparison_name = test_level + DEFAULT_COMPARISON_INFIX + reference_level
            base_file_name = comparison_name
            if pheno_label in phenos_with_replicates:
                lines.append([pheno_label, test_level, reference_level, base_file_name])
    return lines

def _write_tab_delim_file(lines, output_filename):
    with open(output_filename, 'w') as output_file:
        for line in lines:
            print('\t'.join(line), file=output_file)

def build_sample_metadata(config, pheno_with_replicates, sample_metadata_filename):
    _write_tab_delim_file(_build_sample_metadata_list(config, pheno_with_replicates),
                          sample_metadata_filename)

def build_contrasts(config, phenos_with_replicates, contrasts_filename):
    _write_tab_delim_file(_build_contrasts_list(config, phenos_with_replicates), contrasts_filename)

def main(config_filename, sample_metadata_filename, contrasts_filename):
    with open(config_filename, 'r') as config_file:
        config = yaml.load(config_file)
    phenos_with_replicates = PhenotypeManager(config).phenotypes_with_replicates

    build_sample_metadata(config, phenos_with_replicates, sample_metadata_filename)
    build_contrasts(config, phenos_with_replicates, contrasts_filename)

if __name__ == '__main__':
    config_file = sys.argv[1]
    sample_metadata_filename = sys.argv[2]
    contrasts_filename = sys.argv[3]
    main(config_file, sample_metadata_filename, contrasts_filename)
    print('done', file=sys.stderr)