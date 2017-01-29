#!/usr/bin/env python
from __future__ import print_function, absolute_import, division

from collections import defaultdict
import sys

import yaml

from scripts.rnaseq_snakefile_helper import PhenotypeManager

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
    def __init__(self, phenotype_manager, log=sys.stderr):
        self.phenotype_manager = phenotype_manager
        self._validations = [self._comparison_missing_phenotype_value,
                             self._comparison_missing_phenotype_label,
                             self._samples_excluded_from_comparison]
        self._log = log

    def validate(self):
        collector = _ValidationCollector(self._log)
        for validation in self._validations:
            collector.check(validation)
        return collector

    def _comparison_missing_phenotype_value(self):
        comparison_pheno_values = set()
        for pheno_label, pheno_values in self.phenotype_manager.comparison_values.items():
            for pheno_value in pheno_values:
                comparison_pheno_values.add((pheno_label, pheno_value))

        pheno_values = defaultdict(list)
        pheno_samples = defaultdict(list)
        for phenoLabel, phenoValueSamples in self.phenotype_manager.phenotype_sample_list.items():
            for phenoValue, samples in phenoValueSamples.items():
                if (phenoLabel, phenoValue) not in comparison_pheno_values:
                    pheno_values[phenoLabel].append(phenoValue)
                    pheno_samples[phenoLabel].extend(samples)
        if pheno_values:
            label_values = [label +':' + ','.join(sorted(values)) for label, values in sorted(pheno_values.items())]
            samples = [label + ':' + ','.join(sorted(samples)) for label, samples in sorted(pheno_samples.items())]
            msg_fmt = ('Some phenotypes ({}) are not present in comparisons; '
                       'some samples ({}) will be excluded from comparisons for those '
                       'phenotypes.')
            raise _WatermelonConfigWarning(msg_fmt,
                                           ';'.join(label_values),
                                           ';'.join(samples))

    def _comparison_missing_phenotype_label(self):
        comparison_pheno_labels = self.phenotype_manager.comparison_values.keys()
        sample_pheno_labels = self.phenotype_manager.phenotype_sample_list.keys()
        missing_labels = sorted(set(sample_pheno_labels) - set(comparison_pheno_labels))

        if missing_labels:
            msg_fmt = ('Some phenotype labels ({}) are not present in comparisons.')
            raise _WatermelonConfigWarning(msg_fmt, ','.join(missing_labels))

    def _samples_excluded_from_comparison(self):
        all_samples = set(self.phenotype_manager.sample_phenotype_value_dict.keys())
        comparison_pheno_labels = self.phenotype_manager.comparison_values.keys()
        phenotype_sample_list = self.phenotype_manager.phenotype_sample_list
        all_compared_samples = set()
        for (pheno_label, pheno_values) in self.phenotype_manager.comparison_values.items():
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
    exit_code = 0
    with open(config_filename, 'r') as config_file:
        config = yaml.load(config_file)
    phenotype_manager = PhenotypeManager(config)
    validator = _ConfigValidator(phenotype_manager, log=log)
    validation_collector = validator.validate()
    validation_collector.log_results()
    if not validation_collector.ok_to_proceed(prompt_to_override):
        exit_code = 1
    return exit_code

if __name__ == '__main__':
    config_filepath = sys.argv[1]
    exit_code = main(config_filepath)
    exit(exit_code)