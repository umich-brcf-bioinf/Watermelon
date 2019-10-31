#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

from os.path import dirname, exists, join, realpath
import scripts.watermelon_init as watermelon_init
import scripts.annotate as annotate
import unittest
import yaml

_TEST_DIR = realpath(dirname(__file__))
_CONFIG_DIR = realpath(join(_TEST_DIR, '..', 'config'))
_GENOME_CONFIG_FILE_NAME = 'genome_references.yaml'
_GENOME_CONFIG_PATH = join(_CONFIG_DIR, _GENOME_CONFIG_FILE_NAME)

class GenomeTest(unittest.TestCase):

    def test_genomes_references_exist(self):
        with open(_GENOME_CONFIG_PATH, 'r') as config_file:
            genome_config = yaml.load(config_file, Loader=yaml.SafeLoader)
        missing = []
        if 'Other' in genome_config.keys():
            del genome_config['Other'] # TWS - Genome 'Other' is a placeholder - doesn't need to exist
        for genome, values in genome_config.items():
            references = values['references']
            for file_type, file_path in references.items():
                if not exists(file_path):
                    missing.append('{}:{}:{}'.format(genome, file_type, file_path))
        self.assertEqual([],
                         missing,
                         '{}: some references do not exist'.format(_GENOME_CONFIG_FILE_NAME))

    def test_genomes_references_well_formed(self):
        expected_references = set(['fasta', 'gtf', 'annotation_tsv'])
        with open(_GENOME_CONFIG_PATH, 'r') as config_file:
            genome_config = yaml.load(config_file, Loader=yaml.SafeLoader)
        invalid_references = {}
        for genome, values in genome_config.items():
            references = values['references']
            if set(references.keys()) != expected_references:
                invalid_references[genome] = references.keys()
        self.assertEqual({},
                         invalid_references,
                         '{}: reference has missing/extra keys'.format(_GENOME_CONFIG_FILE_NAME))
