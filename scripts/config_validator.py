#!/usr/bin/env python
from __future__ import print_function, absolute_import, division

from collections import defaultdict
from functools import partial
import sys

import yaml

from scripts.rnaseq_snakefile_helper import PhenotypeManager
from scripts import watermelon_config

_HEADER_RULE = '=' * 70 + '\n'
_SECTION_RULE = '-' * 70 + '\n'

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
        self._log = log
        self._primary_validations = [self._check_missing_required_field,
                                     self._check_samples_malformed,
                                     self._check_samples_not_stringlike,
                                     self._check_comparisons_malformed,
                                     self._check_comparisons_not_a_pair,
                                     self._check_phenotypes_samples_not_rectangular,
                                     ]
        self._secondary_validations = [self._check_comparison_missing_phenotype_value,
                                       self._check_comparison_missing_phenotype_label,
                                       self._check_samples_excluded_from_comparison,
                                       ]

    def validate(self):
        collector = _ValidationCollector(self._log)
        for validation in self._secondary_validations:
            collector.check(validation)
        return collector

    def _check_missing_required_field(self):
        missing_fields = watermelon_config.REQUIRED_FIELDS - self.config.keys()
        if missing_fields:
            msg_fmt = ('Some required fields were missing from config: ({}); '
                       'review config and try again.')
            raise _WatermelonConfigFailure(msg_fmt, ', '.join(sorted(missing_fields)))

    def _check_samples_malformed(self):
        msg = ('Config entry [samples] must be formatted as a dict of "{}" separated '
               'strings; review config and try again.'
              ).format(watermelon_config.DEFAULT_PHENOTYPE_DELIM)
        failure = _WatermelonConfigFailure(msg)
        samples = self.config[watermelon_config.CONFIG_KEYS.samples]
        if not isinstance(samples, dict):
            raise failure

    def _check_samples_not_stringlike(self):
        samples = self.config[watermelon_config.CONFIG_KEYS.samples]
        odd_samples = []
        for sample, pheno_value in samples.items():
            if not isinstance(pheno_value, (str, int, float)):
                odd_samples.append(sample)
        if odd_samples:
            msg = ('Some [samples] phenotype values could not be parsed: ({}); '
                   'review config and try again.')
            odd_samples_str = ', '.join(sorted(odd_samples))
            raise _WatermelonConfigFailure(msg, odd_samples_str)

    def _check_comparisons_malformed(self):
        msg = ('Config entry [comparisons] must be formatted as a dict of lists; '
               'review config an try again.')
        failure = _WatermelonConfigFailure(msg)
        comparisons = self.config[watermelon_config.CONFIG_KEYS.comparisons]
        if not isinstance(comparisons, dict):
            raise failure
        for value in comparisons.values():
            if not isinstance(value, list):
                raise failure

    def _check_comparisons_not_a_pair(self):
        comparisons = self.config[watermelon_config.CONFIG_KEYS.comparisons]
        malformed_comparisons = []
        for phenotype_label, comparison_list in comparisons.items():
            for comparison in comparison_list:
                if not comparison:
                    malformed_comparisons.append("[empty comparison]")
                else:
                    values = comparison.strip().split(watermelon_config.DEFAULT_COMPARISON_INFIX)
                    values = [i for i in values if i]
                    if len(values) != 2:
                        malformed_comparisons.append(comparison)
        if malformed_comparisons:
            msg = ('Some [comparisons] are not paired: ({}); '
                   'review config and try again')
            malformed_str = ', '.join(sorted(malformed_comparisons))
            raise _WatermelonConfigFailure(msg, malformed_str)

    def _check_phenotypes_samples_not_rectangular(self):
        pheno_labels = self.config[watermelon_config.CONFIG_KEYS.phenotypes]
        pheno_labels_count = len(pheno_labels.split(watermelon_config.DEFAULT_PHENOTYPE_DELIM))
        sample_pheno_values = self.config[watermelon_config.CONFIG_KEYS.samples]
        problem_samples = defaultdict(int)
        for sample, pheno_values in sample_pheno_values.items():
            pheno_values_count = len(str(pheno_values).split(watermelon_config.DEFAULT_PHENOTYPE_DELIM))
            if pheno_labels_count != pheno_values_count:
                problem_samples[sample] = pheno_values_count
        if problem_samples:
            problems = ['{} [{} values]'.format(s, v) for s,v in sorted(problem_samples.items())]
            problems_str = ', '.join(problems)
            msg = ('Some [samples] had unexpected number of phenotype values '
                   '[expected {} values]: ({}); '
                    'review config and try again.')
            raise _WatermelonConfigFailure(msg, pheno_labels_count, problems_str)

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

def prompt_to_override():
    value = input('Should Watermelon ignore these problems and proceed? (yes/no): ')
    return value.lower().strip() == 'yes'

def main(config_filename,
         log=sys.stderr,
         prompt_to_override=prompt_to_override):
    exit_code = 1
    try:
        with open(config_filename, 'r') as config_file:
            config = yaml.load(config_file)
        validator = _ConfigValidator(config, log=log)
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
    config_filepath = sys.argv[1]
    exit_code = main(config_filepath)
    exit(exit_code)