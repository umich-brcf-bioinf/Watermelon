from __future__ import print_function, absolute_import, division

from argparse import Namespace
from collections import defaultdict
import sys


class WatermelonConfigWarning(Exception):
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(WatermelonConfigWarning, self).__init__(error_msg)


class ConfigValidator(object):
    def __init__(self, phenotype_manager, stderr=sys.stderr):
        self.phenotype_manager = phenotype_manager
        self._validations = [self._comparison_missing_phenotype_value,
                             self._comparison_missing_phenotype_label,
                             self._samples_excluded_from_comparison]

    def validate(self):
        passed = True
        warnings = []
        errors = []
        try:
            for validation in self._validations:
                validation()
        except WatermelonConfigWarning as warning:
            warnings.append(warning)
        if warnings or errors:
            passed = False
        return Namespace(passed=passed, warnings=warnings, errors=errors)

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
            raise WatermelonConfigWarning(msg)

    def _comparison_missing_phenotype_label(self):
        comparison_pheno_labels = self.phenotype_manager.comparison_values.keys()
        sample_pheno_labels = self.phenotype_manager.phenotype_sample_list.keys()
        missing_labels = sorted(set(sample_pheno_labels) - set(comparison_pheno_labels))

        if missing_labels:
            msg_fmt = ('Some phenotype labels ({}) are not present in comparisons.')
            msg = msg_fmt.format(','.join(missing_labels))
            raise WatermelonConfigWarning(msg)

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
            raise WatermelonConfigWarning(msg)