import gzip
import os
import yaml
from collections.abc import Mapping
from glob import glob

# https://stackoverflow.com/a/32357112
def recursive_update(d, u):
    '''recursively update the keys of a dictionary
    d is original dict, u is dict of updated items'''
    for k, v in u.items():
        if isinstance(d, Mapping):
            if isinstance(v, Mapping):
                r = recursive_update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = u[k]
        else:
            d = {k: u[k]}
    return d

# https://stackoverflow.com/a/30921635
def ddict2dict(d):
    '''recursively convert default dict type to dict'''
    for k, v in d.items():
        if isinstance(v, dict):
            d[k] = ddict2dict(v)
    return dict(d)

def create_modified_config(example_file, modified_file, replacements, rm_keys=None):
    with open(example_file, 'r') as example_config_file:
        example_config = yaml.load(example_config_file, Loader=yaml.SafeLoader)
    recursive_update(example_config, replacements)
    if rm_keys:
        [example_config.pop(x, None) for x in rm_keys] #https://stackoverflow.com/a/32616904
    with open(modified_file, 'w') as modified_config_file:
        yaml.dump(example_config, modified_config_file, default_flow_style=False, indent=4)
