#!/bin/env python
from __future__ import print_function, absolute_import, division

import glob
import hashlib
import os
from os.path import join
from os.path import isfile


FILE_EXTENSION = '.watermelon.md5'
'''Appended to all checksum filenames; should be *very* distinctive as this util will
*remove* checksums not found in the config.'''

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
    wildcard = join(checksum_dir, '*{}'.format(FILE_EXTENSION))
    all_checksum_files = set(glob.glob(wildcard))
    extra_checksum_files = all_checksum_files - valid_checksum_files
    for filename in extra_checksum_files:
        with open(filename, 'r') as file:
            first_line = file.readline()
        if first_line.startswith(FILE_EXTENSION):
            os.remove(filename)

def reset_checksums(checksum_dir, config):
    _mkdir(checksum_dir)
    checksum_files = set()
    for key, value in config.iteritems():
        checksum_filepath = _checksum_filepath(checksum_dir, key)
        new_checksum = FILE_EXTENSION + ':' + hashlib.md5(str(value)).hexdigest()
        if not _checksum_matches(checksum_filepath, new_checksum):
            _reset_checksum_file(checksum_filepath, new_checksum)
        checksum_files.add(checksum_filepath)
    _remove_extra_checksum_files(checksum_dir, checksum_files)
