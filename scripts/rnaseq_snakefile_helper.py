'''Creates checksum files that represent the top-level keys of a Snakemake snakefile.
'''
from __future__ import print_function, absolute_import, division

import collections
from collections import defaultdict
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

def _checksum_reset_file(checksum_filepath, new_checksum):
    with open(checksum_filepath, 'w') as file:
        file.write(new_checksum)

def _checksum_remove_extra_files(checksum_dir, valid_checksum_files):
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

def checksum_reset_all(checksum_dir, config):
    '''Create, examine, or update a file for each top-level key in the config dict.
    Each value yields a checksum which is compared with the checksum stored in the file.
    If the file is absent, or the checksum doesn't match, a new checksum file is created.
    Checksum files which are not found in the config keys are removed.'''
    _mkdir(checksum_dir)
    checksum_files = set()
    for key, value in config.items():
        checksum_filepath = _checksum_filepath(checksum_dir, key)
        new_checksum = _build_checksum(value)
        if not _checksum_matches(checksum_filepath, new_checksum):
            _checksum_reset_file(checksum_filepath, new_checksum)
        checksum_files.add(checksum_filepath)
    _checksum_remove_extra_files(checksum_dir, checksum_files)

def init_references(config_references):
    def existing_link_target_is_different(link_name, link_path):
        original_abs_path = os.path.realpath(link_name)
        new_abs_path = os.path.realpath(link_path)
        return original_abs_path != new_abs_path
    if not os.path.exists("references"):
        os.mkdir("references")
    os.chdir("references")
    references = config_references if config_references else {}
    for link_name, link_path in references.items():
        if not os.path.exists(link_path):
            msg_fmt = 'ERROR: specified config reference files/dirs [{}:{}] cannot be read'
            msg = msg_fmt.format(link_name, link_path)
            raise ValueError(msg)
        elif not os.path.exists(link_name):
            os.symlink(link_path, link_name)
        elif existing_link_target_is_different(link_name, link_path):
            os.remove(link_name)
            os.symlink(link_path, link_name)
        else:
            pass #link matches existing link
    os.chdir("..")

def cuffdiff_conditions(comparison_infix, explicit_comparisons):
    unique_conditions = set()
    for comparison in explicit_comparisons:
        unique_conditions.update(comparison.split(comparison_infix))
    multi_condition_comparison = comparison_infix.join(sorted(unique_conditions))
    return(multi_condition_comparison)


def check_strand_option(library_type,strand_option):

    tuxedo_library_type = { 0 : "fr-unstranded", 1 : "fr-firststrand", 2 : "fr-secondstrand"}
    htseq_library_type = { 0 : "no", 1 : "yes", 2 : "reverse"}
    options = range(0,3)
# 
#     TUXEDO = 'tuxedo'
#     strand_config_param = { 0 : {'tuxedo': 'fr-unstranded', 'htseq': 'no'}, }
# 
#     try:
#         param_value = strand_config_param[config_param][library_type]
#     except KeyError:
#         raise KeyError('whatup with yout config')

    if not strand_option in options:
        msg_format = "ERROR: 'alignment_options:library_type' in config file is '{}'. Valid library_type options are: 0 (un-stranded), 1 (first-strand), or 2 (second-strand)."
        msg = msg_format.format(strand_option)
        raise ValueError(msg)

    if library_type == "tuxedo":
        return tuxedo_library_type[strand_option]
    elif library_type == "htseq":
        return htseq_library_type[strand_option]
    else:
        msg_format = "ERROR: Unknown library_type {}."
        msg = msg_format.format(library_type)
        raise ValueError(msg)


def cuffdiff_labels(comparison_infix, underbar_separated_comparisons):
    return underbar_separated_comparisons.replace(comparison_infix, ",")

def cuffdiff_samples(comparison_infix,
                     underbar_separated_comparisons,
                     sample_name_group,
                     sample_file_format):
    group_sample_names = defaultdict(list)
    for actual_sample_name, group in sample_name_group.items():
        group_sample_names[group].append(sample_file_format.format(sample_placeholder=actual_sample_name))
    group_sample_names = dict(group_sample_names)
    params = []
    for group in underbar_separated_comparisons.split(comparison_infix):
        params.append(','.join(sorted(group_sample_names[group])))
    return ' '.join(params)

def cutadapt_options(trim_params):
    run_trimming_options = 0
    for option, value in trim_params.items():
        if not isinstance(value, int):
            msg_format = "ERROR: Config trimming_options '{}' must be integer"
            msg = msg_format.format(value)
            raise ValueError(msg)
        run_trimming_options += value
    return run_trimming_options


def tophat_options(alignment_options):
    options = ""
    if not isinstance(alignment_options["transcriptome_only"], bool):
        raise ValueError("config alignment_options:transcriptome_only must be boolean")
    if alignment_options["transcriptome_only"]:
        options += " --transcriptome-only "
    else:
        options += " --no-novel-juncs "  # used in Legacy for transcriptome + genome alignment
    return options
