from argparse import Namespace
CONFIG_KEYS = Namespace(comparisons='comparisons',
                        main_factors='main_factors',
                        phenotypes='phenotypes',
                        samples='samples',
                        )
DEFAULT_COMPARISON_INFIX = '_v_'
DEFAULT_PHENOTYPE_DELIM = '^'

MAIN_FACTOR_TRUE = 'yes'

REQUIRED_FIELDS = set([CONFIG_KEYS.phenotypes,
                       CONFIG_KEYS.samples,
                       CONFIG_KEYS.comparisons])
