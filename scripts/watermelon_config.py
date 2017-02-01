from argparse import Namespace
CONFIG_KEYS = Namespace(phenotypes='phenotypes',
                        samples='samples',
                        comparisons='comparisons')
DEFAULT_COMPARISON_INFIX = '_v_'
DEFAULT_PHENOTYPE_DELIM = '^'

REQUIRED_FIELDS = set([CONFIG_KEYS.phenotypes,
                       CONFIG_KEYS.samples,
                       CONFIG_KEYS.comparisons])
