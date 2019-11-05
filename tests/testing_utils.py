import collections
import gzip
import os
import yaml
from glob import glob

# https://stackoverflow.com/a/32357112
def recursive_update(d, u):
    for k, v in u.items():
        if isinstance(d, collections.Mapping):
            if isinstance(v, collections.Mapping):
                r = recursive_update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = u[k]
        else:
            d = {k: u[k]}
    return d

def create_modified_config(example_file, modified_file, replacements):
    with open(example_file, 'r') as example_config_file:
        example_config = yaml.load(example_config_file, Loader=yaml.SafeLoader)
    recursive_update(example_config, replacements)
    with open(modified_file, 'w') as modified_config_file:
        yaml.dump(example_config, modified_config_file, default_flow_style=False, indent=4)

def gunzip(source_file_pattern):
    def _gunzip_file(source_file, dest_file):
        with gzip.GzipFile(source_file, 'rb') as inF, \
             open(dest_file, 'wb') as outF:
            data = inF.read()
            outF.write(data)
            os.remove(source_file)
    for source_filename in glob(source_file_pattern):
        dest_filename = source_filename.rstrip('.gz')
        _gunzip_file(source_filename, dest_filename)
