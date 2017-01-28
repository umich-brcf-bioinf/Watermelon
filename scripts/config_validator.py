from __future__ import print_function, absolute_import, division

from collections import defaultdict
import sys


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
    def __init__(self):
        self.failures = []
        self.warnings = []

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
            return '{} ({}):'.format(summary, ','.join(results))
        else:
            return summary

    def log_results(self, log):
        log.write('config validation: {}\n'.format(self.summary_result))

class _ConfigValidator(object):
    def __init__(self, phenotype_manager, stderr=sys.stderr):
        self.phenotype_manager = phenotype_manager
        self._validations = [self._comparison_missing_phenotype_value,
                             self._comparison_missing_phenotype_label,
                             self._samples_excluded_from_comparison]

    def validate(self):
        collector = _ValidationCollector()
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
            msg = msg_fmt.format(';'.join(label_values), ';'.join(samples))
            raise _WatermelonConfigWarning(msg)

    def _comparison_missing_phenotype_label(self):
        comparison_pheno_labels = self.phenotype_manager.comparison_values.keys()
        sample_pheno_labels = self.phenotype_manager.phenotype_sample_list.keys()
        missing_labels = sorted(set(sample_pheno_labels) - set(comparison_pheno_labels))

        if missing_labels:
            msg_fmt = ('Some phenotype labels ({}) are not present in comparisons.')
            msg = msg_fmt.format(','.join(missing_labels))
            raise _WatermelonConfigWarning(msg)

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
            msg = msg_fmt.format(','.join(missing_samples))
            raise _WatermelonConfigWarning(msg)



def main(config_filename):
    pass