'''Functions that support DESeq2.
'''
from __future__ import print_function, absolute_import, division

from scripts.watermelon_config import CONFIG_KEYS, DEFAULT_PHENOTYPE_DELIM, MAIN_FACTOR_TRUE

SAMPLE_METADATA_HEADER = 'phenotype'
COMBINATORIC_GROUP = 'combinatoric_group'

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
    header = [SAMPLE_METADATA_HEADER]
    header.extend(phenotype_labels)
    header.append(COMBINATORIC_GROUP)
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
    return []

def build_sample_metadata(config, sample_metadata_file_name):
    with open(sample_metadata_file_name, 'w') as metadata_file:
        lines = _build_sample_metadata_list(config)
        for line in lines:
            print('\t'.join(line), file=metadata_file)

def build_contrasts_file(config, output_base_path, sample_metadata_file):
    pass

