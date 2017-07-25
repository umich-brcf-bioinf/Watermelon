from argparse import Namespace
CONFIG_KEYS = Namespace(input_dir='input_dir',
                        comparisons='comparisons',
                        main_factors='main_factors',
                        phenotypes='phenotypes',
                        samples='samples',
                        genome='genome',
                        references='references',
                        )
DEFAULT_COMPARISON_INFIX = '_v_'
DEFAULT_PHENOTYPE_DELIM = '^'

MAIN_FACTOR_TRUE = 'yes'
MAIN_FACTOR_FALSE = 'no'

MAIN_FACTOR_VALID_VALUES = set([MAIN_FACTOR_TRUE, MAIN_FACTOR_FALSE])

REQUIRED_FIELDS = set([CONFIG_KEYS.comparisons,
                       CONFIG_KEYS.main_factors,
                       CONFIG_KEYS.phenotypes,
                       CONFIG_KEYS.samples,])

def split_config_list(config_string, delim=DEFAULT_PHENOTYPE_DELIM):
    if type(config_string) != str:
        if config_string:
            config_string = MAIN_FACTOR_TRUE
        else:
            config_string = MAIN_FACTOR_FALSE
    config_list = []
    for i in config_string.split(delim):
        config_list.append(i.strip())
    return config_list