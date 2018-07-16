from argparse import Namespace
CONFIG_KEYS = Namespace(dirs='dirs',
                        dirs_input='input',
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

def transform_config(config):
    '''Accepts old or new config dict updating obsolete keys in place.'''
    def _check_not_both(config_keys, dirs):
        old_new_pairs = []
        for old_key, new_key in config_keys.items():
            if old_key in config and new_key in dirs:
                old_new_pairs.append((old_key, 'dirs:' + new_key))
        if old_new_pairs:
            msg = ('found both old and new style dirs '
                   '({}); remove old config line(s) ({}) and try again')
            pairs = '; '.join(sorted(['{},{}'.format(x,y) for x,y in old_new_pairs]))
            old = '; '.join(sorted([x for x,_ in old_new_pairs]))
            raise ValueError(msg.format(pairs, old))

    config_keys = {'input_dir': 'input',
                   'alignment_output_dir': 'alignment_output',
                   'diffex_output_dir': 'diffex_output',
                   'deliverables_output_dir': 'deliverables_output',
                }
    dirs = config.get('dirs', {})
    _check_not_both(config_keys, dirs)
    for old_key, new_key in config_keys.items():
        if old_key in config:
            print('remapping ' + old_key)
            dirs[new_key] = config.pop(old_key)
    if dirs:
        config['dirs'] = dirs
