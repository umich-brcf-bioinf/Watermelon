#!/usr/bin/env python

from collections import Counter,defaultdict
from functools import partial
import os
import re
import sys

import yaml
import pandas as pd
import pdb # TWS DEBUG

from scripts.rnaseq_snakefile_helper import PhenotypeManager

_HEADER_RULE = '=' * 70 + '\n'
_SECTION_RULE = '-' * 70 + '\n'

_LEGAL_NAME_REGEX = r'^[A-Za-z][A-Za-z0-9\-_\.]*$'
_RESERVED_NAMES = set(['T', 'F', 'NAN'])


def _is_name_well_formed(name):
    if not re.match(_LEGAL_NAME_REGEX, name):
        return False
    return True

def _is_name_reserved(name):
    return name in _RESERVED_NAMES


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
        for failure in self.failures:
            self._log.write('failure: {}\n'.format(failure))
        for warning in self.warnings:
            self._log.write('warning: {}\n'.format(warning))

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
    def __init__(self, config, log=sys.stderr):
        self.config = config
        self.samplesheet = pd.read_csv(config['sample_description_file'])
        self._log = log
        self._PARSING_VALIDATIONS = [self._check_phenotype_labels_blank,
                                     self._check_phenotype_labels_not_unique,
                                     self._check_comparisons_malformed,
                                     self._check_comparisons_not_a_pair,
                                     ]
        self._CONTENT_VALIDATIONS = [self._check_phenotype_labels_illegal_values,
                                     self._check_phenotype_labels_reserved_name,
                                     self._check_phenotype_values_illegal_values,
                                     self._check_phenotype_values_reserved_name,
                                     self._check_comparison_missing_phenotype_value,
                                     self._check_comparison_missing_phenotype_label,
                                     self._check_comparison_references_unknown_phenotype_label,
                                     self._check_comparison_references_unknown_phenotype_value,
                                     self._check_samples_excluded_from_comparison,
                                     self._check_phenotype_has_replicates]

    def validate(self):
        collector = _ValidationCollector(self._log)

        for validation in self._PARSING_VALIDATIONS:
            collector.check(validation)
            if not collector.passed:
                return collector

        for validation in self._CONTENT_VALIDATIONS:
            collector.check(validation)
        return collector

    def _check_phenotype_has_replicates(self):
        phenotype_manager = PhenotypeManager(self.config)
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

    def _check_comparisons_malformed(self):
        msg = ('Config entry [comparisons] must be formatted as a dict of lists; '
               'review config an try again.')
        failure = _WatermelonConfigFailure(msg)
        comparisons = self.config['comparisons']
        if not isinstance(comparisons, dict):
            raise failure
        for value in comparisons.values():
            if not isinstance(value, list):
                raise failure

    def _check_comparison_values_distinct(self):
        comparisons = self.config['comparisons']
        nondistinct_comparisons = []
        for phenotype_label, comparison_list in comparisons.items():
            for comparison in comparison_list:
                if comparison:
                    values = comparison.strip().split("foo") #TWS FIXME
                    values = [i for i in values if i]
                    if len(values) == 2 and values[0]==values[1]:
                        nondistinct_comparisons.append(comparison)
        if nondistinct_comparisons:
            msg = ('Some [comparison] test-control values are not distinct: ({}); '
                   'review config and try again')
            problem_str = ', '.join(sorted(nondistinct_comparisons))
            raise _WatermelonConfigFailure(msg, problem_str)

    def _check_comparisons_not_a_pair(self):
        comparisons = self.config['comparisons']
        malformed_comparisons = []
        for phenotype_label, comparison_list in comparisons.items():
            for comparison in comparison_list:
                if not comparison:
                    malformed_comparisons.append("[empty comparison]")
                else:
                    values = comparison.strip().split("foo") #TWS FIXME
                    values = [i for i in values if i]
                    if len(values) != 2:
                        malformed_comparisons.append(comparison)
        if malformed_comparisons:
            msg = ('Some [comparisons] are not paired: ({}); '
                   'review config and try again')
            malformed_str = ', '.join(sorted(malformed_comparisons))
            raise _WatermelonConfigFailure(msg, malformed_str)

    def _check_phenotype_labels_blank(self):
        #pheno_labels = watermelon_config.split_config_list(self.config['phenotypes']) #TWS FIXME
        for label in pheno_labels:
            if not label:
                raise _WatermelonConfigFailure('[phenotypes] labels cannot be blank')

    def _check_phenotype_labels_not_unique(self):
        #pheno_labels = watermelon_config.split_config_list(self.config[watermelon_config.CONFIG_KEYS.phenotypes]) #TWS FIXME
        duplicate_labels = sorted([label for label, freq in Counter(pheno_labels).items() if freq > 1])
        if duplicate_labels:
            msg = ('[phenotypes] labels must be unique; review/revise [{}]')
            raise _WatermelonConfigFailure(msg, ', '.join(duplicate_labels))

    def _check_comparison_references_unknown_phenotype_label(self):
        phenotype_manager = PhenotypeManager(self.config)
        comparison_pheno_labels = set(phenotype_manager.comparison_values.keys())
        phenotype_labels = set(phenotype_manager.phenotype_sample_list.keys())
        unknown_pheno_labels = sorted(comparison_pheno_labels - phenotype_labels)
        if unknown_pheno_labels:
            msg = ('Some [comparisons] referenced phenotype labels not found in '
                   '[phenotypes]: ({}); review config and try again.'
                  )
            raise _WatermelonConfigFailure(msg, ', '.join(unknown_pheno_labels))

    def _check_comparison_references_unknown_phenotype_value(self):
        pheno_label_values = []
        phenotype_manager = PhenotypeManager(self.config)
        for label, values in phenotype_manager.phenotype_sample_list.items():
            for value in values:
                pheno_label_values.append((label, value))
        comparison_label_values = []
        for label, values in phenotype_manager.comparison_values.items():
            for value in values:
                comparison_label_values.append((label, value))
        unknown_pheno_values = sorted(set(comparison_label_values) - set(pheno_label_values))
        if unknown_pheno_values:
            msg = ('Some [comparisons] referenced phenotype values not found in '
                   '[samples]: ({}); review config and try again.'
                  )
            pheno_value_strings = ['{}:{}'.format(l,v)for (l,v) in unknown_pheno_values]
            raise _WatermelonConfigFailure(msg, ', '.join(pheno_value_strings))


    def _check_comparison_missing_phenotype_value(self):
        phenotype_manager = PhenotypeManager(self.config)
        comparison_pheno_values = set()
        for pheno_label, pheno_values in phenotype_manager.comparison_values.items():
            for pheno_value in pheno_values:
                comparison_pheno_values.add((pheno_label, pheno_value))

        pheno_values = defaultdict(list)
        pheno_samples = defaultdict(list)
        for phenoLabel, phenoValueSamples in phenotype_manager.phenotype_sample_list.items():
            for phenoValue, samples in phenoValueSamples.items():
                if (phenoLabel, phenoValue) not in comparison_pheno_values:
                    pheno_values[phenoLabel].append(phenoValue)
                    pheno_samples[phenoLabel].extend(samples)
        if pheno_values:
            label_values = [label +':' + ','.join(sorted(values)) for label, values in sorted(pheno_values.items())]
            samples = [label + ':' + ','.join(sorted(samples)) for label, samples in sorted(pheno_samples.items())]
            msg_fmt = ('Some phenotype values ({}) are not present in [comparisons]; '
                       'some samples ({}) will be excluded from comparisons for those '
                       'phenotypes.')
            raise _WatermelonConfigWarning(msg_fmt,
                                           ';'.join(label_values),
                                           ';'.join(samples))

    def _check_comparison_missing_phenotype_label(self):
        phenotype_manager = PhenotypeManager(self.config)
        comparison_pheno_labels = phenotype_manager.comparison_values.keys()
        sample_pheno_labels = phenotype_manager.phenotype_sample_list.keys()
        missing_labels = sorted(set(sample_pheno_labels) - set(comparison_pheno_labels))

        if missing_labels:
            msg_fmt = ('Some phenotype labels ({}) are not present in [comparisons].')
            raise _WatermelonConfigWarning(msg_fmt, ','.join(missing_labels))

    def _check_samples_excluded_from_comparison(self):
        phenotype_manager = PhenotypeManager(self.config)
        all_samples = set(phenotype_manager.sample_phenotype_value_dict.keys())
        comparison_pheno_labels = phenotype_manager.comparison_values.keys()
        phenotype_sample_list = phenotype_manager.phenotype_sample_list
        all_compared_samples = set()
        for (pheno_label, pheno_values) in phenotype_manager.comparison_values.items():
            for pheno_value in pheno_values:
                samples = phenotype_sample_list[pheno_label][pheno_value]
                all_compared_samples.update(samples)
        missing_samples = sorted(all_samples - all_compared_samples)
        if missing_samples:
            msg_fmt = 'Some samples ({}) will not be compared.'
            raise _WatermelonConfigWarning(msg_fmt,
                                           ','.join(missing_samples))

    def _check_main_factors_illegal_values(self):
        #actual_values = set(watermelon_config.split_config_list(self.config[watermelon_config.CONFIG_KEYS.main_factors])) #TWS FIXME
        #legal_values = watermelon_config.MAIN_FACTOR_VALID_VALUES #TWS FIXME
        illegal_values = sorted(actual_values - legal_values)
        illegal_values = ['<blank>' if v=='' else v for v in illegal_values]
        if illegal_values:
            msg_fmt = 'Expected [main_factor] values to be ({}) but found ({})';
            raise _WatermelonConfigFailure(msg_fmt,
                                           ', '.join(sorted(legal_values)),
                                           ', '.join(illegal_values))

    def _check_phenotype_labels_illegal_values(self):
        bad_labels = []
        #actual_labels = watermelon_config.split_config_list(self.config[watermelon_config.CONFIG_KEYS.phenotypes]) #TWS FIXME
        for label in actual_labels:
            if not _is_name_well_formed(label):
                bad_labels.append(label)
        bad_labels = sorted(set(['<blank>' if v=='' else v for v in bad_labels]))
        if bad_labels:
            msg_fmt = ('[phenotypes] labels must begin with a letter and can contain '
                       'only (A-Za-z0-9-_.); review/revise these labels [{}]')
            raise _WatermelonConfigFailure(msg_fmt, ', '.join(bad_labels))

    def _check_phenotype_labels_reserved_name(self):
        bad_labels = []
        #actual_labels = watermelon_config.split_config_list(self.config[watermelon_config.CONFIG_KEYS.phenotypes]) #TWS FIXME
        for label in actual_labels:
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
        phenotype_manager = PhenotypeManager(self.config)
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
        phenotype_manager = PhenotypeManager(self.config)
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
    value = input('Should Watermelon ignore these problems and proceed? (yes/no): ')
    return value.lower().strip() == 'yes'

def main(config_filename,
         log=sys.stderr,
         prompt_to_override=prompt_to_override):
    exit_code = 1
    try:
        with open(config_filename, 'r') as config_file:
            config = yaml.load(config_file, Loader=yaml.SafeLoader)
        validator = _ConfigValidator(config, log=log)
        pdb.set_breakpoint()
        validation_collector = validator.validate()
        validation_collector.log_results()
        if validation_collector.ok_to_proceed(prompt_to_override):
            exit_code = 0
    except yaml.YAMLError:
        log.write(_HEADER_RULE)
        log.write('config validation: ERROR\n')
        log.write(('Could not parse this config file [{}]; '
                   'is it valid YAML?\n').format(config_filename))
    return exit_code

if __name__ == '__main__':
    config_filepath = os.path.realpath(sys.argv[1])
    exit_code = main(config_filepath)
    exit(exit_code)
