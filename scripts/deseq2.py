'''Functions that support DESeq2.
'''
from __future__ import print_function, absolute_import, division

from scripts.watermelon_config import CONFIG_KEYS, DEFAULT_PHENOTYPE_DELIM

SAMPLE_METADATA_HEADER = 'phenotype'

def _split_config_list(config_string):
    config_list = []
    for i in config_string.split(DEFAULT_PHENOTYPE_DELIM):
        config_list.append(i.strip())
    return config_list

def _build_sample_metadata_list(config):
    lines = []
    phenotype_labels = _split_config_list(config[CONFIG_KEYS.phenotypes])
    header = [SAMPLE_METADATA_HEADER]
    header.extend(phenotype_labels)
    lines.append(header)
    for sample_name, phenotype_values in sorted(config[CONFIG_KEYS.samples].items()):
        sample_line = [sample_name]
        sample_line.extend(_split_config_list(phenotype_values))
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

