'''Creates checksum files that represent the top-level keys of a Snakemake snakefile.
'''
from __future__ import print_function, absolute_import, division

import collections
import glob
import hashlib
import os
from os.path import join
from os.path import isfile


FILE_EXTENSION = '.watermelon.md5'
'''Appended to all checksum filenames; must be very distinctive to ensure the module can
unambiguously identify and safely remove extraneous checksums.'''

def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)

def _checksum_filepath(checksum_dir, key):
    checksum_filename = '{}{}'.format(key, FILE_EXTENSION)
    return join(checksum_dir, checksum_filename)

def _checksum_matches(checksum_filepath, new_checksum):
    checksums_match = False
    if os.path.exists(checksum_filepath):
        with open(checksum_filepath, 'r') as file:
            old_checksum = file.readlines()
        checksums_match = len(old_checksum) == 1 and old_checksum[0] == new_checksum
    return checksums_match

def _reset_checksum_file(checksum_filepath, new_checksum):
    with open(checksum_filepath, 'w') as file:
        file.write(new_checksum)

def _remove_extra_checksum_files(checksum_dir, valid_checksum_files):
    '''Removes extraneous checksum files (checking for proper name and checksum prefix).'''
    wildcard = join(checksum_dir, '*{}'.format(FILE_EXTENSION))
    all_checksum_files = set(glob.glob(wildcard))
    extra_checksum_files = all_checksum_files - valid_checksum_files
    for filename in extra_checksum_files:
        with open(filename, 'r') as file:
            lines = file.readlines()
        if len(lines) == 1 and lines[0].startswith(FILE_EXTENSION):
            os.remove(filename)

def _build_checksum(value):
    '''Creates a checksum based on the supplied value. Typically just the md5 of the
    str representation, but if the value is itself a dict will create a str representation
    of an ordered dict. (This extra step avoids non-deterministic behavior around how
    dicts are ordered between/within python2/3, but the naive implementation assumes that
    the config will only be two dicts deep at most.)
    '''
    if isinstance(value, dict):
        value = collections.OrderedDict(sorted(value.items()))
    return FILE_EXTENSION + ':' + hashlib.md5(str(value).encode('utf-8')).hexdigest()

def reset_checksums(checksum_dir, config):
    '''Create, examine, or update a file for each top-level key in the config dict.
    Each value yields a checksum which is compared with the checksum stored in the file.
    If the file is absent, or the checksum doesn't match a new checksum file is created.
    Checksum files which are not found in the config keys are removed.'''
    _mkdir(checksum_dir)
    checksum_files = set()
    for key, value in config.items():
        checksum_filepath = _checksum_filepath(checksum_dir, key)
        new_checksum = _build_checksum(value)
        if not _checksum_matches(checksum_filepath, new_checksum):
            _reset_checksum_file(checksum_filepath, new_checksum)
        checksum_files.add(checksum_filepath)
    _remove_extra_checksum_files(checksum_dir, checksum_files)
