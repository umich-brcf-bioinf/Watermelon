#!/usr/bin/env python

from collections import Counter,defaultdict
from functools import partial
from itertools import chain
import jsonschema
import os
import re
import sys
#from snakemake import utils # Using utils.validate, keeping the utils name to prevent confusion

import yaml
import pandas as pd

if __package__ == None or __package__ == "": # If not loaded as a module
    # helper is colocated with this validator
    import rnaseq_snakefile_helper
    import __init__ as WAT_VER
else:
    from . import rnaseq_snakefile_helper
    from .__init__ import __version__ as WAT_VER

_HEADER_RULE = '=' * 70 + '\n'
_SECTION_RULE = '-' * 70 + '\n'

_LEGAL_NAME_REGEX = r'^[A-Za-z][A-Za-z0-9_\.]*$'
_RESERVED_NAMES = set(['T', 'F', 'NAN'])


def _is_name_well_formed(name):
    if not re.match(_LEGAL_NAME_REGEX, name):
        return False
    return True

def _is_name_reserved(name):
    return name in _RESERVED_NAMES


'''Returns contrast string dict in the form of
{factor_name: ['val1_v_val2', 'val3_v_val4']}'''
def _DESeq2_factor_contrasts(diffex_config):
    cont_dict = {}
    for model in rnaseq_snakefile_helper.diffex_models(diffex_config):
        factor = diffex_config[model]['DESeq2']['factor_name']
        cont_dict[factor] = diffex_config[model]['contrasts']
    return(cont_dict)

'''Returns contrast values dict in the form of
{factor_name: ['val1', 'val2', 'val3', 'val4']}
'''
def _DESeq2_contrast_vals(contrast_dict):
    contrast_vals_dict = defaultdict(set)
    for pheno in contrast_dict:
        cont_strs = contrast_dict[pheno]
        cont_vals = chain.from_iterable(map(lambda x: x.split("_v_"), cont_strs))
        contrast_vals_dict[pheno].update(cont_vals)
    return(contrast_vals_dict)

class _WatermelonConfigFailure(Exception):
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(_WatermelonConfigFailure, self).__init__(error_msg)


class _WatermelonConfigWarning(Exception):
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(_WatermelonConfigWarning, self).__init__(error_msg)


class _ValidationCollector(object):
    def __init__(self, log):
        self.failures = []
        self.warnings = []
        self._log = log

    def check(self, validation):
        try:
            validation()
        except _WatermelonConfigWarning as warning:
            self.warnings.append(warning)
        except _WatermelonConfigFailure as failure:
            self.failures.append(failure)

    @property
    def passed(self):
        if self.warnings or self.failures:
            return False
        else:
            return True

    @property
    def summary_result(self):
        results = []
        summary = 'OK'
        if self.warnings:
            summary = 'WARNING'
            results.append('{} warnings'.format(len(self.warnings)))
        if self.failures:
            summary = 'FAILED'
            results.insert(0, '{} failures'.format(len(self.failures)))
        if results:
            return '{} ({}):'.format(summary, ', '.join(results))
        else:
            return summary

    def log_results(self):
        self._log.write(_HEADER_RULE)
        self._log.write('config validation: {}\n'.format(self.summary_result))
        if not self.passed:
            self._log.write(_SECTION_RULE)
        for warning in self.warnings:
            self._log.write('warning: {}\n'.format(warning))
        for failure in self.failures:
            self._log.write('failure: {}\n'.format(failure))

    def ok_to_proceed(self, prompt_to_override_warnings):
        result = True
        if self.failures:
            self._log.write(_HEADER_RULE)
            self._log.write('There were config file validation failures; please review/revise the config and try again.\n')
            result = False
        elif self.warnings:
            self._log.write(_HEADER_RULE)
            self._log.write('There were config file validation warnings.\n')
            result = prompt_to_override_warnings()
            if result:
                self._log.write('Based on user input Watermelon will continue despite warnings above.\n')
        if not result:
            self._log.write('Watermelon stopped to review config file warnings/errors.\n')
        self._log.write(_HEADER_RULE)
        return result


class _ConfigValidator(object):
    def __init__(self, config, schema, log=sys.stderr):
        self.config = config
        self.schema = schema
        try:
            self.samplesheet = pd.read_csv(config['samplesheet'], comment='#', dtype='string')
        except:
            print('problem reading samplesheet')
        self.contrasts = _DESeq2_factor_contrasts(config['diffex'])
        self.contrast_values = _DESeq2_contrast_vals(self.contrasts)
        self._log = log
        self._PARSING_VALIDATIONS = [self._check_config_against_schema,
                                     self._check_config_matches_watermelon_version,
                                     self._check_phenotype_labels_not_unique,
                                     self._check_contrasts_not_a_pair]
        self._CONTENT_VALIDATIONS = [self._check_phenotype_labels_illegal_values,
                                     self._check_phenotype_labels_reserved_name,
                                     self._check_phenotype_values_illegal_values,
                                     self._check_phenotype_values_reserved_name,
                                     self._check_contrast_values_distinct,
                                     self._check_contrast_missing_phenotype_value,
                                     self._check_contrast_missing_phenotype_label,
                                     self._check_contrast_references_unknown_phenotype_label,
                                     self._check_contrast_references_unknown_phenotype_value,
                                     self._check_samples_excluded_from_contrast,
                                     self._check_phenotype_has_replicates]

    def validate(self):
        collector = _ValidationCollector(self._log)

        for validation in self._PARSING_VALIDATIONS:
            collector.check(validation)
            #TWS - Why is this here (below)? The result is that only
            #parsing errors/warnings are reported if any of these fail
            #if not collector.passed:
            #    return collector

        for validation in self._CONTENT_VALIDATIONS:
            collector.check(validation)
        return collector

    def _check_config_against_schema(self):
        try:
            jsonschema.validate(self.config, self.schema) #snakemake utils schema validator
        except Exception as e:
            msg = e.message #This is given by jsonschema.exceptions.ValidationError
            raise(_WatermelonConfigFailure(str(msg)))

    def _check_config_matches_watermelon_version(self):
        #TWS - Request code review here
        try:
            config_ver = self.config['watermelon_version']
        except KeyError as e:
            msg = 'Key {} not found in config. Cannot check version match'.format(e)
            raise _WatermelonConfigFailure(msg)
        if config_ver != WAT_VER:
            msg = 'Current watermelon version: {} does not match config: {}'.format(WAT_VER, config_ver)
            raise _WatermelonConfigWarning(msg)


    def _check_phenotype_has_replicates(self):
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        all_phenotypes = set(phenotype_manager.phenotype_sample_list.keys())
        with_replicates = set(phenotype_manager.phenotypes_with_replicates)
        if not with_replicates:
            msg = ('No phenotype labels have any replicates; DESeq2 '
                   'will not be run.')
            raise _WatermelonConfigWarning(msg)
        without_replicates = all_phenotypes - with_replicates
        if without_replicates:
            msg_fmt = ('Some phenotype labels ({}) have no replicates and will be '
                       'excluded from DESeq2 results.')
            raise _WatermelonConfigWarning(msg_fmt,
                                           ', '.join(sorted(without_replicates)))

    def _check_contrast_values_distinct(self):
        contrasts = self.contrasts
        nondistinct_contrasts = []
        for phenotype_label, contrast_list in contrasts.items():
            for contrast in contrast_list:
                if contrast:
                    values = contrast.strip().split("_v_")
                    values = [i for i in values if i]
                    if len(values) == 2 and values[0]==values[1]:
                        nondistinct_contrasts.append(contrast)
        if nondistinct_contrasts:
            msg = ('Some [contrast] test-control values are not distinct: ({}); '
                   'review config and try again')
            problem_str = ', '.join(sorted(nondistinct_contrasts))
            raise _WatermelonConfigFailure(msg, problem_str)

    def _check_contrasts_not_a_pair(self):
        contrasts = self.contrasts
        malformed_contrasts = []
        for phenotype_label, contrast_list in contrasts.items():
            for contrast in contrast_list:
                if not contrast:
                    malformed_contrasts.append("[empty contrast]")
                else:
                    values = contrast.strip().split("_v_")
                    values = [i for i in values if i]
                    if len(values) != 2:
                        malformed_contrasts.append(contrast)
        if malformed_contrasts:
            msg = ('Some [contrasts] are not paired: ({}); '
                   'review config and try again')
            malformed_str = ', '.join(sorted(malformed_contrasts))
            raise _WatermelonConfigFailure(msg, malformed_str)

    def _check_phenotype_labels_not_unique(self):
        pheno_labels = list(self.samplesheet.columns)
        duplicate_labels = sorted([label for label, freq in Counter(pheno_labels).items() if freq > 1])
        if duplicate_labels:
            msg = ('[phenotypes] labels must be unique; review/revise [{}]')
            raise _WatermelonConfigFailure(msg, ', '.join(duplicate_labels))

    def _check_contrast_references_unknown_phenotype_label(self):
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        contrast_pheno_labels = set(self.contrast_values.keys())
        phenotype_labels = set(phenotype_manager.phenotype_sample_list.keys())
        unknown_pheno_labels = sorted(contrast_pheno_labels - phenotype_labels)
        if unknown_pheno_labels:
            msg = ('Some [contrasts] referenced phenotype labels not found in '
                   '[phenotypes]: ({}); review config and try again.'
                  )
            raise _WatermelonConfigFailure(msg, ', '.join(unknown_pheno_labels))

    def _check_contrast_references_unknown_phenotype_value(self):
        pheno_label_values = []
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        for label, values in phenotype_manager.phenotype_sample_list.items():
            for value in values.keys():
                pheno_label_values.append((label, value))
        contrast_label_values = []

        for label, values in self.contrast_values.items():
            for value in values:
                contrast_label_values.append((label, value))
        unknown_pheno_values = sorted(set(contrast_label_values) - set(pheno_label_values))
        if unknown_pheno_values:
            msg = ('Some [contrasts] referenced phenotype values not found in '
                   '[samples]: ({}); review config and try again.'
                  )
            pheno_value_strings = ['{}:{}'.format(l,v)for (l,v) in unknown_pheno_values]
            raise _WatermelonConfigFailure(msg, ', '.join(pheno_value_strings))


    def _check_contrast_missing_phenotype_value(self):
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        contrast_pheno_values = set()
        for pheno_label, pheno_values in self.contrast_values.items():
            for pheno_value in pheno_values:
                contrast_pheno_values.add((pheno_label, pheno_value))

        pheno_values = defaultdict(list)
        pheno_samples = defaultdict(list)
        for phenoLabel, phenoValueSamples in phenotype_manager.phenotype_sample_list.items():
            for phenoValue, samples in phenoValueSamples.items():
                if (phenoLabel, phenoValue) not in contrast_pheno_values:
                    pheno_values[phenoLabel].append(phenoValue)
                    pheno_samples[phenoLabel].extend(samples)
        if pheno_values:
            label_values = [label +':' + ','.join(sorted(values)) for label, values in sorted(pheno_values.items())]
            samples = [label + ':' + ','.join(sorted(samples)) for label, samples in sorted(pheno_samples.items())]
            msg_fmt = ('Some phenotype values are not present in [contrasts].\n'
                       'Missing phenotype values: \n({})\n')
            raise _WatermelonConfigWarning(msg_fmt,
                                           ';'.join(label_values))

    def _check_contrast_missing_phenotype_label(self):
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        contrast_pheno_labels = self.contrasts.keys()
        sample_pheno_labels = phenotype_manager.phenotype_sample_list.keys()
        missing_labels = sorted(set(sample_pheno_labels) - set(contrast_pheno_labels))

        if missing_labels:
            msg_fmt = ('Some phenotype labels ({}) are not present in [contrasts].')
            raise _WatermelonConfigWarning(msg_fmt, ','.join(missing_labels))

    def _check_samples_excluded_from_contrast(self):
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        all_samples = set(phenotype_manager.sample_phenotype_value_dict.keys())
        contrast_pheno_labels = self.contrast_values.keys()
        phenotype_sample_list = phenotype_manager.phenotype_sample_list
        all_compared_samples = set()
        for (pheno_label, pheno_values) in self.contrast_values.items():
            for pheno_value in pheno_values:
                samples = phenotype_sample_list[pheno_label][pheno_value]
                all_compared_samples.update(samples)
        missing_samples = sorted(all_samples - all_compared_samples)
        if missing_samples:
            msg_fmt = 'Some samples ({}) will not be compared.'
            raise _WatermelonConfigWarning(msg_fmt,
                                           ','.join(missing_samples))

    def _check_phenotype_labels_illegal_values(self):
        bad_labels = []
        pheno_labels = list(self.samplesheet.columns)
        for label in pheno_labels:
            if not _is_name_well_formed(label):
                bad_labels.append(label)
        bad_labels = sorted(set(['<blank>' if v=='' else v for v in bad_labels]))
        if bad_labels:
            msg_fmt = ('[phenotypes] labels must begin with a letter and can contain '
                       'only (A-Za-z0-9_.); review/revise these labels [{}]')
            raise _WatermelonConfigFailure(msg_fmt, ', '.join(bad_labels))

    def _check_phenotype_labels_reserved_name(self):
        bad_labels = []
        pheno_labels = list(self.samplesheet.columns)
        for label in pheno_labels:
            if _is_name_reserved(label):
                bad_labels.append(label)
        bad_labels = sorted(set(bad_labels))
        if bad_labels:
            msg_fmt = ('[phenotypes] labels cannot be a reserved word ({}); '
                       'review/revise these labels [{}]')
            raise _WatermelonConfigFailure(msg_fmt,
                                           ', '.join(_RESERVED_NAMES),
                                           ', '.join(bad_labels))

    def _check_phenotype_values_illegal_values(self):
        label_values = []
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        for label, values in phenotype_manager.phenotype_sample_list.items():
            for value in values:
                label_values.append((label, value))
        bad_label_values = []
        for label,value in label_values:
            if value and label and not _is_name_well_formed(value):
                bad_label_values.append((label, value))
        bad_label_value_strings = ['{}:{}'.format(l,v) for l,v in sorted(set(bad_label_values))]
        if bad_label_value_strings:
            msg_fmt = ('[sample] phenotype values must begin with a letter and can contain '
                       'only (A-Za-z0-9-_.); review/revise these values [{}]')
            raise _WatermelonConfigFailure(msg_fmt,
                                           ', '.join(bad_label_value_strings))

    def _check_phenotype_values_reserved_name(self):
        label_values = []
        phenotype_manager = rnaseq_snakefile_helper.PhenotypeManager(self.config)
        for label, values in phenotype_manager.phenotype_sample_list.items():
            for value in values:
                label_values.append((label, value))
        bad_label_values = []
        for label,value in label_values:
            if value and label and _is_name_reserved(value):
                bad_label_values.append((label, value))
        bad_label_value_strings = ['{}:{}'.format(l,v) for l,v in sorted(set(bad_label_values))]
        if bad_label_value_strings:
            msg_fmt = ('[sample] phenotype values cannot be a reserved word ({}); '
                       'review/revise these values [{}]')
            raise _WatermelonConfigFailure(msg_fmt,
                                           ', '.join(_RESERVED_NAMES),
                                           ', '.join(bad_label_value_strings))

def prompt_to_override():
    while True:
        value = input('Should Watermelon ignore these problems and proceed? (yes/no): ')

        if value not in ['yes', 'no']:
            print("Response must be either 'yes' or 'no'")
            continue
        else:
            break

    return value.lower().strip() == 'yes'

def main(config_filename,
         schema_filename,
         log=sys.stderr,
         prompt_to_override=prompt_to_override):

    exit_code = 1 #Assume that it doesn't pass

    try:
        with open(config_filename, 'r') as config_file:
            config = yaml.load(config_file, Loader=yaml.SafeLoader)
    except (yaml.parser.ParserError, yaml.scanner.ScannerError, yaml.YAMLError) as e:
        log.write(_HEADER_RULE)
        log.write(('Error: Could not parse the following config file [{}]\n'
                   'Verify that YAML is valid.\n').format(config_filename))
        return exit_code
    try:
        with open(schema_filename, 'r') as schema_file:
            schema = yaml.load(schema_file, Loader=yaml.SafeLoader)
    except (yaml.parser.ParserError, yaml.scanner.ScannerError, yaml.YAMLError) as e:
        log.write(_HEADER_RULE)
        log.write(('Error: Could not parse the following schema file [{}]\n'
                   'Verify that YAML is valid.\n').format(schema_filename))
        return exit_code

    validator = _ConfigValidator(config, schema, log=log)
    validation_collector = validator.validate()
    validation_collector.log_results()
    if validation_collector.ok_to_proceed(prompt_to_override):
        exit_code = 0

    return exit_code

if __name__ == '__main__':
    config_filepath = os.path.realpath(sys.argv[1])
    schema_filepath = os.path.realpath(sys.argv[2])
    exit_code = main(config_filepath, schema_filepath)
    exit(exit_code)
