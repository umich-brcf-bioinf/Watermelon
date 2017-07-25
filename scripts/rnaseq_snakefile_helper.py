'''Functions that support the RNASeq Snakemake snakefile.
'''
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import collections
from collections import defaultdict
import csv
from functools import partial
import glob
import hashlib
import os
from os.path import join
from os.path import isfile
import sys

from snakemake import workflow

from scripts.watermelon_config import CONFIG_KEYS
from scripts.watermelon_config import DEFAULT_COMPARISON_INFIX
from scripts.watermelon_config import DEFAULT_PHENOTYPE_DELIM

FILE_EXTENSION = '.watermelon.md5'
'''Appended to all checksum filenames; must be very distinctive to ensure the module can
unambiguously identify and safely remove extraneous checksums.'''

TOPHAT_NAME = 'tophat'
HTSEQ_NAME = 'htseq'

STRAND_CONFIG_PARAM = { 'fr-unstranded' : {TOPHAT_NAME: 'fr-unstranded', HTSEQ_NAME: 'no'},
                        'fr-firststrand' : {TOPHAT_NAME : 'fr-firststrand', HTSEQ_NAME: 'reverse'},
                        'fr-secondstrand' : {TOPHAT_NAME : 'fr-secondstrand', HTSEQ_NAME: 'yes'} }

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
    '''Interprets a subset of the config to help answer questions around how
    samples map to phenotype labels and values and vice versa.'''
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
            sample_phenotype_values = list(map(str.strip, phenotype_values.split(self.delimiter)))
            sample_phenotypes = dict(zip(phenotype_labels, sample_phenotype_values))
            check_labels_match_values(phenotype_labels,
                                      sample,
                                      sample_phenotype_values)
            for label, value in sample_phenotypes.items():
                if value != '':
                    phenotype_dict[label][value].append(sample)
        return phenotype_dict

    @property
    def phenotypes_with_replicates(self):
        with_replicates = []
        for (label, values) in self.phenotype_sample_list.items():
            max_sample_count = max([len(samples) for (value, samples) in values.items()])
            if max_sample_count > 1:
                with_replicates.append(label)
        return sorted(with_replicates)

    def separated_comparisons(self, output_delimiter):
        '''Returns a dict of {phenotype_label: 'val1_v_val2,val3_v_val4'}'''
        comparisons = {}
        for phenotype, comparison in self.comparisons.items():
            comp_string = output_delimiter.join(list(map(str.strip, comparison)))
            comparisons[phenotype]= comp_string
        return comparisons

    @property
    def comparison_values(self):
        '''Returns dict of {pheno_label : [val1,val2,val3,val4] }'''
        phenotype_label_values = {}
        for phenotype, comparison_list in self.comparisons.items():
            comparison_list = list(map(str.strip, comparison_list))
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
    def phenotypes_comparisons_all_tuple(self):
        '''Returns named tuple of phenotype, comparison for all phenotypes'''
        return self._phenotypes_comparisons(True)

    @property
    def phenotypes_comparisons_replicates_tuple(self):
        '''Returns named tuple of phenotype, comparison for phenotypes with replicates'''
        return self._phenotypes_comparisons(False)

    def _phenotypes_comparisons(self, include_phenotypes_without_replicates=True):
        if include_phenotypes_without_replicates:
            include = lambda p: True
        else:
            with_replicates = set(self.phenotypes_with_replicates)
            include = lambda p: p in with_replicates
        pheno_comps = [(p, c) for p,c in self.comparisons.items() if include(p)]
        phenotypes=[]
        comparisons=[]
        for pheno, comps in sorted(pheno_comps):
            for comp in sorted(comps):
                phenotypes.append(pheno)
                comparisons.append(comp.strip())
        PhenotypesComparisons = collections.namedtuple('PhenotypeComparisons',
                                                       'phenotypes comparisons')
        return PhenotypesComparisons(phenotypes=phenotypes, comparisons=comparisons)


    def cuffdiff_samples(self,
                         phenotype_label,
                         sample_file_format):
        '''Returns a list of sample files grouped by phenotype value.
           Each sample group represents the comma separated samples for a phenotype value.
           Sample groups are separated by a space.

           sample_file_format is format string with placeholder {sample_placeholder}'''

        group_separator = ' '
        file_separator = ','

        sample_name_group = self.phenotype_sample_list[phenotype_label]

        group_sample_names = defaultdict(list)
        for phenotype_value, sample_list in sample_name_group.items():
            for sample_name in sample_list:
                sample_file = sample_file_format.format(sample_placeholder=sample_name)
                group_sample_names[phenotype_value].append(sample_file)
        group_sample_names = dict(group_sample_names)

        params = []
        for group in self.comparison_values[phenotype_label]:
            params.append(file_separator.join(sorted(group_sample_names[group])))

        return group_separator.join(params)

def _strand_option(tuxedo_or_htseq, strand_option):
    if not strand_option in STRAND_CONFIG_PARAM:
        msg_format = ('ERROR: config:alignment_options:library_type={} is '
                      'not valid. Valid library_type options are: {}')
        msg = msg_format.format(strand_option, ','.join(STRAND_CONFIG_PARAM.keys()))
        raise ValueError(msg)
    return STRAND_CONFIG_PARAM[strand_option][tuxedo_or_htseq]

def strand_option_tophat(config_strand_option):
    return _strand_option(TOPHAT_NAME, config_strand_option)

def strand_option_htseq(config_strand_option):
    return _strand_option(HTSEQ_NAME, config_strand_option)

def cutadapt_options(trim_params):
    run_trimming_options = 0
    for option, value in trim_params.items():
        if not isinstance(value, int):
            msg_format = "ERROR: config:trimming_options '{}={}' must be integer"
            msg = msg_format.format(option, value)
            raise ValueError(msg)
        if value != 0:
            run_trimming_options = 1
    return run_trimming_options

def tophat_options(alignment_options):
    options = ""
    if not isinstance(alignment_options["transcriptome_only"], bool):
        raise ValueError("config:alignment_options:transcriptome_only must be True or False")
    if alignment_options["transcriptome_only"]:
        options += " --transcriptome-only "
    else:
        options += " --no-novel-juncs "  # used in Legacy for transcriptome + genome alignment
    return options

def _get_sample_reads(fastq_base_dir, samples):
    sample_reads = {}
    read_suffix = {0: "", 1:"_SE", 2:"_PE"}
    is_fastq = lambda fn: fn.endswith('.fastq') or fn.endswith('.fastq.gz')
    for sample in samples:
        found_reads = []
        for read in ['R1', 'R2']:
            sample_dir = join(fastq_base_dir, sample, '')
            read_present = list(filter(is_fastq, glob.glob(sample_dir + '*_{}*'.format(read))))
            if read_present:
                found_reads.append(read)
        found_reads = list(map(lambda x: x+read_suffix[len(found_reads)], found_reads))
        sample_reads[sample] = found_reads
    return sample_reads

def flattened_sample_reads(fastq_base_dir, samples):
    sample_reads = _get_sample_reads(fastq_base_dir, samples)
    return sorted([(sample,read) for (sample,reads) in sample_reads.items() for read in reads])

def expand_sample_read_endedness(sample_read_endedness_format,
                                 all_flattened_sample_reads,
                                 sample=None):
    if not all_flattened_sample_reads:
        return []
    samples, reads = zip(*[(s, r) for s, r in all_flattened_sample_reads if not sample or s == sample])
    return workflow.expand(sample_read_endedness_format,
                           zip,
                           sample=samples,
                           read_endedness=reads)

def tophat_paired_end_flags(read_stats_filename=None):
    def read_values(filename, expected_headers):
        with open(filename, mode='r') as infile:
            try:
                row = next(csv.DictReader(infile, delimiter='\t'))
            except StopIteration:
                raise ValueError('missing data row in [{}]'.format(filename))
        missing_headers = sorted(expected_headers - row.keys())
        if missing_headers:
            msg = 'missing [{}] in header of [{}]'.format(','.join(missing_headers),
                                                          read_stats_filename)
            raise ValueError(msg)
        return row

    def parse_values(filename, row, headers):
        parsed_row = {}
        invalid_values = {}
        for header in headers:
            value = row[header]
            try:
                parsed_row[header] = str(int(round(float(value))))
            except ValueError:
                invalid_values[header]=value
        if invalid_values:
            header_values = ','.join(['{}:{}'.format(h,v) for h,v in sorted(invalid_values.items())])
            msg = 'invalid number(s) [{}] in {}'.format(header_values, read_stats_filename)
            raise ValueError(msg)
        return parsed_row

    header_flag = {'insert_std_dev': '--mate-std-dev',
                   'inner_mate_dist': '--mate-inner-dist'}
    result = []
    if read_stats_filename and os.path.exists(read_stats_filename):
        raw_values = read_values(read_stats_filename, header_flag.keys())
        parsed_values = parse_values(read_stats_filename,
                                     raw_values,
                                     header_flag.keys())
        for header, flag in sorted(header_flag.items()):
            result.append(flag)
            result.append(parsed_values[header])
    return ' '.join(result)

def expand_read_stats_if_paired(read_stats_filename_format,
                                flattened_sample_reads,
                                sample):
    result = []
    if len([s for s,r in flattened_sample_reads if s == sample]) > 1:
        result.append(read_stats_filename_format.format(sample=sample))
    return result
