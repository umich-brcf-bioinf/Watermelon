'''Functions that support the RNASeq Snakemake snakefile.
'''
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import collections
from collections import defaultdict
from functools import partial
import glob
import hashlib
import os
from os.path import join
from os.path import isfile


CONFIG_KEYS = Namespace(phenotypes='phenotypes',
                        samples='samples',
                        comparisons='comparisons')
DEFAULT_COMPARISON_INFIX = '_v_'
DEFAULT_PHENOTYPE_DELIM = '^'

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

class ChecksumManager(object):
    def __init__(self, checksum_dir, file_extension):
        self.checksum_dir = checksum_dir
        self.file_extension = file_extension

    def _checksum_filepath(self, prefix, key):
        checksum_filename = '{}-{}{}'.format(prefix, key, self.file_extension)
        return join(self.checksum_dir, checksum_filename)

    @staticmethod
    def _checksum_matches(checksum_filepath, new_checksum):
        checksums_match = False
        if os.path.exists(checksum_filepath):
            with open(checksum_filepath, 'r') as file:
                old_checksum = file.readlines()
            checksums_match = len(old_checksum) == 1 and old_checksum[0] == new_checksum
        return checksums_match

    @staticmethod
    def _checksum_reset_file(checksum_filepath, new_checksum):
        with open(checksum_filepath, 'w') as file:
            file.write(new_checksum)

    def _checksum_remove_extra_files(self, valid_checksum_files):
        '''Removes extraneous checksum files (checking for proper name and checksum prefix).'''
        wildcard = join(self.checksum_dir, '*{}'.format(self.file_extension))
        all_checksum_files = set(glob.glob(wildcard))
        extra_checksum_files = all_checksum_files - valid_checksum_files
        for filename in extra_checksum_files:
            with open(filename, 'r') as file:
                lines = file.readlines()
            if len(lines) == 1 and lines[0].startswith(self.file_extension):
                os.remove(filename)

    @staticmethod
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

    def reset_checksum_for_dict(self, the_name, the_dict):
        '''Checksum of each value is compared with the checksum stored in the file; if
        the checksum doesn't match, create a new checksum file.'''
        checksum_files = set()
        for key, value in the_dict.items():
            checksum_filepath = self._checksum_filepath(the_name, key)
            new_checksum = self._build_checksum(value)
            if not self._checksum_matches(checksum_filepath, new_checksum):
                self._checksum_reset_file(checksum_filepath, new_checksum)
            checksum_files.add(checksum_filepath)
        return checksum_files

    def reset_checksums(self, **kwargs):
        '''Create, examine, or update a file for each top-level key in the dictionaries.
        Stray checksum files are removed.'''
        _mkdir(self.checksum_dir)
        checksum_files = set()
        for the_name, the_dict in kwargs.items():
            checksum_files.update(tuple(self.reset_checksum_for_dict(the_name, the_dict)))
        self._checksum_remove_extra_files(checksum_files)

def checksum_reset_all(checksum_dir, **kwargs):
    ChecksumManager(checksum_dir, FILE_EXTENSION).reset_checksums(**kwargs)

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

class PhenotypeManager(object):
    def __init__(self,
                 config={},
                 delimiter=DEFAULT_PHENOTYPE_DELIM,
                 comparison_infix=DEFAULT_COMPARISON_INFIX):
        self.phenotype_labels_string = config.get(CONFIG_KEYS.phenotypes, None)
        self.sample_phenotype_value_dict = config.get(CONFIG_KEYS.samples, None)
        self.comparisons = config.get(CONFIG_KEYS.comparisons, None)
        self.delimiter = delimiter
        self.comparison_infix = comparison_infix

    @property
    def phenotype_sample_list(self):
        '''Translates config phenotypes/samples into nested dict of phenotypes.
        Specifically {phenotype_label : {phenotype_value : [list of samples] } }

        Strips all surrounding white space.
    
        phenotypes_string : delimited phenotype labels (columns)
        sample_phenotype_value_dict : {sample_id : delimited phenotype_value_string} (rows)
        '''
        def check_labels_match_values(phenotype_labels,
                                      sample,
                                      phenotype_values):
            if len(phenotype_labels) != len(phenotype_values):
                msg_fmt = 'expected {} phenotype values but sample {} had {}'
                msg = msg_fmt.format(len(phenotype_labels),
                                     sample,
                                     len(phenotype_values))
                raise ValueError(msg)

        def check_phenotype_labels(labels):
            for i,label in enumerate(labels):
                if label == '':
                    msg_fmt = 'label of phenotype {} is empty'
                    msg = msg_fmt.format(i + 1)
                    raise ValueError(msg)

        phenotype_dict = defaultdict(partial(defaultdict, list))
        phenotype_labels = list(map(str.strip, self.phenotype_labels_string.split(self.delimiter)))
        check_phenotype_labels(phenotype_labels)
        sorted_sample_phenotypes_items = sorted([(k, v) for k,v in self.sample_phenotype_value_dict.items()])
        for sample, phenotype_values in sorted_sample_phenotypes_items:
            sample_phenotype_values = map(str.strip, phenotype_values.split(self.delimiter))
            sample_phenotypes = dict(zip(phenotype_labels, sample_phenotype_values)) 
            check_labels_match_values(phenotype_labels,
                                      sample,
                                      sample_phenotype_values)
            for label, value in sample_phenotypes.items():
                if value != '':
                    phenotype_dict[label][value].append(sample)
        return phenotype_dict

    def separated_comparisons(self, output_delimiter):
        '''Returns a dict of {phenotype_label: 'val1_v_val2,val3_v_val4'}'''
        comparisons = {}
        for phenotype, comparison in self.comparisons.items():
            comp_string = output_delimiter.join(list(map(str.strip, comparison.split())))
            comparisons[phenotype]= comp_string
        return comparisons

    @property
    def comparison_values(self):
        '''Returns dict of {pheno_label : [val1,val2,val3,val4] }'''
        phenotype_label_values = {}
        for phenotype, comparison in self.comparisons.items():
            comparison_list = list(map(str.strip, comparison.split()))
            unique_conditions = set()
            for group in comparison_list:
                unique_conditions.update(group.split(self.comparison_infix))
            values = sorted(unique_conditions) 
            phenotype_label_values[phenotype] = values
        return phenotype_label_values

    def concatenated_comparison_values(self, output_delimiter):
        '''Returns dict of {pheno_label : 'val1,val2,val3,val4' }'''
        phenotype_label_values = {}
        for phenotype_label, values in self.comparison_values.items():
            phenotype_label_values[phenotype_label] = output_delimiter.join(values)
        return phenotype_label_values

    @property
    def phenotypes_comparisons_tuple(self):
        '''Returns named tuple of phenotype, comparison for all comparisons'''
        phenotype_comparison = dict([(k, v.split()) for k,v in self.comparisons.items()])
        phenotypes=[]
        comparisons=[]
        for k, v in phenotype_comparison.items():
            for c in v:
                phenotypes.append(k)
                comparisons.append(c)
        PhenotypesComparisons = collections.namedtuple('PhenotypeComparisons',
                                                       'phenotypes comparisons')

        return PhenotypesComparisons(phenotypes=phenotypes, comparisons=comparisons)

    def cuffdiff_samples(self,
                         phenotype,
                         sample_file_format):
        '''Returns a list of sample files grouped by phenotype value.
           Each sample group represents the comma separated samples for a phenotype value.
           Sample groups are separated by a space.
           
           sample_file_format is format string with placeholder {sample_placeholder}'''

        group_separator = ' '
        file_separator = ','
        sample_name_group = self.phenotype_sample_list[phenotype]

        group_sample_names = defaultdict(list)
        for phenotype, sample_list in sample_name_group.items():
            for sample_name in sample_list:
                sample_file = sample_file_format.format(sample_placeholder=sample_name)
                group_sample_names[phenotype].append(sample_file)
        group_sample_names = dict(group_sample_names)

        params = []
        for group in self.comparison_values:
            params.append(file_separator.join(sorted(group_sample_names[group])))
        return group_separator.join(params)


#TODO: (cgates): is this used?
# def cuffdiff_conditions(comparison_infix, explicit_comparisons):
#     unique_conditions = set()
#     for comparison in explicit_comparisons:
#         unique_conditions.update(comparison.split(comparison_infix))
#     multi_condition_comparison = comparison_infix.join(sorted(unique_conditions))
#     return(multi_condition_comparison)

def check_strand_option(library_type,strand_option):
    strand_config_param = { 'fr-unstranded' : {'tuxedo': 'fr-unstranded', 'htseq': 'no'}, 
                            'fr-firststrand' : {'tuxedo' : 'fr-firststrand', 'htseq': 'yes'}, 
                            'fr-secondstrand' : {'tuxedo' : 'fr-secondstrand', 'htseq': 'reverse'} }

    if not strand_option in strand_config_param.keys():
        msg_format = "ERROR: 'alignment_options:library_type' in config file is '{}'. Valid library_type options are: 'fr-unstranded', 'fr-firststrand', 'fr-secondstrand'."
        msg = msg_format.format(strand_option)
        raise ValueError(msg)

    try:
        param_value = strand_config_param[strand_option][library_type]
        return param_value
    except KeyError:
        raise KeyError('invalid library type option: ', library_type)

def cuffdiff_samples(comparison_infix,
                     underbar_separated_comparisons,
                     sample_name_group,
                     sample_file_format):
    group_sample_names = defaultdict(list)
    for phenotype, sample_list in sample_name_group.items():
        for sample_name in sample_list:
            group_sample_names[phenotype].append(sample_file_format.format(sample_placeholder=sample_name))
    group_sample_names = dict(group_sample_names)

    params = []
    for group in underbar_separated_comparisons.split(comparison_infix):
        params.append(','.join(sorted(group_sample_names[group])))
    return ' '.join(params)


def cutadapt_options(trim_params):
    run_trimming_options = 0
    for option, value in trim_params.items():
        if not isinstance(value, int):
            raise ValueError(msg)
            msg_format = "ERROR: Config trimming_options '{}' must be integer"
            msg = msg_format.format(value)
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
