'''Functions that support DESeq2.
'''
from __future__ import print_function, absolute_import, division
import os
import sys

import yaml

from scripts.watermelon_config import CONFIG_KEYS, DEFAULT_PHENOTYPE_DELIM, MAIN_FACTOR_TRUE, DEFAULT_COMPARISON_INFIX

_COMBINATORIC_GROUP = 'combinatoric_group'
_CONTRASTS_HEADER = ['factor','test_level','reference_level','directory_name', 'base_file_name']
_SAMPLE_METADATA_HEADER = 'phenotype'

def _split_config_list(config_string):
    config_list = []
    for i in config_string.split(DEFAULT_PHENOTYPE_DELIM):
        config_list.append(i.strip())
    return config_list

def _build_combinatoric_group_value(phenotype_labels, main_factors, phenotype_values):
    main_factors = [label for label, main_factor in zip(phenotype_labels, main_factors) if main_factor == MAIN_FACTOR_TRUE]
    phenotype_label_value = dict([(label, value) for label, value in zip(phenotype_labels, phenotype_values)])
    main_factor_phenotype_values = []
    for phenotype_label in main_factors:
        main_factor_phenotype_values.append(phenotype_label_value[phenotype_label])
    combinatoric_string = DEFAULT_PHENOTYPE_DELIM.join(main_factor_phenotype_values)
    return combinatoric_string

def _build_sample_metadata_list(config):
    lines = []
    main_factors = _split_config_list(config[CONFIG_KEYS.main_factors])
    phenotype_labels = _split_config_list(config[CONFIG_KEYS.phenotypes])
    header = [_SAMPLE_METADATA_HEADER]
    header.extend(phenotype_labels)
    header.append(_COMBINATORIC_GROUP)
    lines.append(header)
    for sample_name, phenotype_values_string in sorted(config[CONFIG_KEYS.samples].items()):
        sample_line = [sample_name]
        phenotype_values = _split_config_list(phenotype_values_string)
        sample_line.extend(phenotype_values)
        combinatoric_group_value = _build_combinatoric_group_value(phenotype_labels,
                                                                   main_factors,
                                                                   phenotype_values)
        sample_line.append(combinatoric_group_value)
        lines.append(sample_line)
    return lines

def _build_contrasts_list(config, output_base_path):
    lines = []
    lines.append(_CONTRASTS_HEADER)
    for pheno_label, comparison_strings in sorted(config[CONFIG_KEYS.comparisons].items()):
        comparisons = sorted([x.split(DEFAULT_COMPARISON_INFIX) for x in comparison_strings])
        for test_level, reference_level in comparisons:
            comparison_name = test_level + DEFAULT_COMPARISON_INFIX + reference_level
            directory_name = os.path.join(output_base_path, pheno_label)
            base_file_name = comparison_name
            lines.append([pheno_label, test_level, reference_level, directory_name, base_file_name])
    return lines

def _write_tab_delim_file(lines, output_filename):
    with open(output_filename, 'w') as output_file:
        for line in lines:
            print('\t'.join(line), file=output_file)

def build_sample_metadata(config, sample_metadata_filename):
    _write_tab_delim_file(_build_sample_metadata_list(config), sample_metadata_filename)

def build_contrasts(config, comparison_file_prefix, contrasts_filename):
    _write_tab_delim_file(_build_contrasts_list(config, comparison_file_prefix), contrasts_filename)

def main(config_filename, sample_metadata_filename, comparison_file_prefix, contrasts_filename):
    with open(config_filename, 'r') as config_file:
        config = yaml.load(config_file)
    build_sample_metadata(config, sample_metadata_filename)
    build_contrasts(config, comparison_file_prefix, contrasts_filename)

if __name__ == '__main__':
    config_file = sys.argv[1]
    sample_metadata_filename = sys.argv[2]
    contrast_comparison_file_prefix = sys.argv[3]
    contrasts_filename = sys.argv[4]
    main(config_file, sample_metadata_filename, contrast_comparison_file_prefix, contrasts_filename)
    print('done', file=sys.stderr)