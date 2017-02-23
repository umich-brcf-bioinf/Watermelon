'''Functions that support DESeq2.
'''
from __future__ import print_function, absolute_import, division

def _split_config_list(config_string):
    config_list = []
    for i in config_string.split('^'):
        config_list.append(i.strip())
    return config_list

def _build_sample_metadata_list(config):
    lines = []

    phenotype_labels = _split_config_list(config['phenotypes'])
    header = ['phenotype']
    header.extend(phenotype_labels)
    
    lines.append(header)
    for sample_name, phenotype_values in sorted(config['samples'].items()):
        sample_line = [sample_name]
        sample_line.extend(_split_config_list(phenotype_values))
        lines.append(sample_line)
    return lines

def _build_contrasts_list(config, output_base_path):
    return []

def build_sample_metadata(config, sample_metadata_file):
    pass

def build_contrasts_file(config, output_base_path, sample_metadata_file):
    pass

