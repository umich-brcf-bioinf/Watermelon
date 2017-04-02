#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import glob
import os
import time
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures.tempdirectory import TempDirectory
import yaml

from scripts.config_validator import _ConfigValidator
from scripts.config_validator import _WatermelonConfigFailure
from scripts.config_validator import _WatermelonConfigWarning
from scripts.config_validator import _ValidationCollector
from scripts.config_validator import main

class MockValidation(object):
    def __init__(self, warning_or_error=None):
        self.validate_was_called = False
        self.warning_or_error = warning_or_error

    def validate(self, *argv):
        self.validate_was_called = True
        if self.warning_or_error:
            raise self.warning_or_error


class MockPromptOverride(object):
    def __init__(self, prompt_to_override_return=False):
        self.prompt_to_override_was_called = False
        self.prompt_to_override_return = prompt_to_override_return

    def prompt_to_override(self):
        self.prompt_to_override_was_called = True
        return self.prompt_to_override_return


class MockValidationCollector(object):
    def __init__(self):
        self.check_calls = []

    def check(self, validation):
        self.check_calls.append(validation)

def _write_config_file(path, filename, contents):
    config_filename = os.path.join(path, filename)
    with open(config_filename, 'w') as config_file:
        print(contents, file=config_file)
    return config_filename

class ConfigValidatorTest(unittest.TestCase):

    def ok(self):
        self.assertEqual(1, 1)

    def test_validate_primary_ok(self):
        validator = _ConfigValidator({})
        validation1 = MockValidation()
        validation2 = MockValidation()
        validator._PRIMARY_VALIDATIONS = [validation1.validate,
                                          validation2.validate]
        validator._SECONDARY_VALIDATIONS = []

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(True, validation_result.passed)
        self.assertEqual([], validation_result.failures)
        self.assertEqual([], validation_result.warnings)

    def test_validate_primaryFailureFailsFast(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        failure =_WatermelonConfigFailure("I'm angry")
        primaryValidation1 = MockValidation(failure)
        primaryValidation2 = MockValidation()
        secondaryValidation1 = MockValidation()

        validator._PRIMARY_VALIDATIONS = [primaryValidation1.validate,
                                          primaryValidation2.validate]
        validator._SECONDARY_VALIDATIONS = [secondaryValidation1.validate]

        validation_result = validator.validate()

        self.assertEqual(True, primaryValidation1.validate_was_called)
        self.assertEqual(False, primaryValidation2.validate_was_called)
        self.assertEqual(False, secondaryValidation1.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([], validation_result.warnings)
        self.assertEqual([failure], validation_result.failures)


    def test_validate_secondary_ok(self):
        validator = _ConfigValidator({})
        validation1 = MockValidation()
        validation2 = MockValidation()
        validator._PRIMARY_VALIDATIONS = []
        validator._SECONDARY_VALIDATIONS = [validation1.validate,
                                            validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(True, validation_result.passed)
        self.assertEqual([], validation_result.failures)
        self.assertEqual([], validation_result.warnings)

    def test_validate_secondary_warning(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        validation1 = MockValidation()
        warning =_WatermelonConfigWarning("I'm angry")
        validation2 = MockValidation(warning)
        validator._PRIMARY_VALIDATIONS = []
        validator._SECONDARY_VALIDATIONS = [validation1.validate,
                                            validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([warning], validation_result.warnings)
        self.assertEqual([], validation_result.failures)


    def test_validate_secondaryWarningContinuesProcessingValidations(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        warning =_WatermelonConfigWarning("I'm angry")
        validation1 = MockValidation(warning)
        validation2 = MockValidation()
        validator._PRIMARY_VALIDATIONS = []
        validator._SECONDARY_VALIDATIONS = [validation1.validate,
                                            validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([warning], validation_result.warnings)
        self.assertEqual([], validation_result.failures)

    def test_validate_secondaryFailure(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        validation1 = MockValidation()
        error =_WatermelonConfigFailure("I'm angry")
        validation2 = MockValidation(error)
        validator._PRIMARY_VALIDATIONS = []
        validator._SECONDARY_VALIDATIONS = [validation1.validate,
                                            validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([], validation_result.warnings)
        self.assertEqual([error], validation_result.failures)


    def test_validate_secondaryMixOfWarningsAndErrors(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        problem1 =_WatermelonConfigFailure("I'm angry")
        validation1 = MockValidation(problem1)
        problem2 =_WatermelonConfigWarning("I'm angry")
        validation2 = MockValidation(problem2)
        problem3 =_WatermelonConfigWarning("I'm angry")
        validation3 = MockValidation(problem3)
        problem4 =_WatermelonConfigFailure("I'm angry")
        validation4 = MockValidation(problem4)
        validator._PRIMARY_VALIDATIONS = []
        validator._SECONDARY_VALIDATIONS = [validation1.validate,
                                            validation2.validate,
                                            validation3.validate,
                                            validation4.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([problem2, problem3], validation_result.warnings)
        self.assertEqual([problem1, problem4], validation_result.failures)

    def test_check_missing_required_field_ok(self):
        config_string = \
'''main_factors: X
phenotypes: X
samples: X
comparisons: X
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_missing_required_field()
        self.ok()

    def test_check_missing_required_field_raises(self):
        config_string = 'foo: bar'
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               (r'Some required fields were missing from config: '
                                r'\(comparisons, main_factors, phenotypes, samples\); '
                                r'review config and try again.'),
                               validator._check_missing_required_field)

    def test_check_samples_malformed_ok(self):
        config_string = \
'''samples:
    s1: 1
    s2: b
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        validator._check_samples_malformed()
        self.ok()

    def test_check_samples_malformed_raisesIfList(self):
        config_string = \
'''samples:
   - s1
   - s2 
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               r'\[samples\].*as a dict of.*strings',
                               validator._check_samples_malformed)

    def test_check_samples_malformed_raisesIfStr(self):
        config_string = \
'''samples:
    s1
    s2 
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               r'\[samples\].*as a dict of.*strings',
                               validator._check_samples_malformed)

    def test_check_samples_not_stringlike_ok(self):
        config_string = \
'''samples:
    s1: A
    s2: 42
    s3: 42.42
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        validator._check_samples_not_stringlike()
        self.ok()

    def test_check_samples_not_stringlike_raisesIfNotStringlike(self):
        config_string = \
'''samples:
    s1:
        a: b
    s2: A
    s3: 42
    s4: 42.42
    s5:
    - b
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               r'Some \[samples\] phenotype values could not be parsed: '
                               r'\(s1, s5\); review config and try again.',
                               validator._check_samples_not_stringlike)

    def test_check_comparisons_malformed_ok(self):
        config_string = \
'''comparisons:
    foo:
    - x1_v_x2
    bar:
    - y1_v_y2
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_comparisons_malformed()
        self.ok()

    def test_check_comparisons_malformed_raisesIfNotDict(self):
        config_string = \
'''comparisons:
    - foo_v_bar
    - bar_v_foo
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               r'comparisons.*as a dict of lists',
                               validator._check_comparisons_malformed)

    def test_check_comparisons_malformed_raisesIfNotList(self):
        config_string = \
'''comparisons:
    foo:
        foo_v_bar
        bar_v_foo
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               r'comparisons.*as a dict of lists',
                               validator._check_comparisons_malformed)

    def test_check_comparisons_not_a_pair_ok(self):
        config_string = \
'''comparisons:
    foo:
    - x1_v_x2
    - x3_v_x4
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_comparisons_not_a_pair()
        self.ok()

    def test_check_comparisons_not_a_pair_raises(self):
        config_string = \
'''comparisons:
    foo:
    - ok1_v_ok2
    - bar_foo
    - x_v_
    - _v_
    - a_v_b_v_c
    - 
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               (r'\[comparisons\] are not paired: \(\[empty comparison\], '
                                r'_v_, a_v_b_v_c, bar_foo, x_v_\);'),
                               validator._check_comparisons_not_a_pair)

    def test_check_phenotypes_samples_not_rectangular_ok(self):
        config_string = \
'''phenotypes: a ^ b
samples:
    s1: 1 ^ 2
    s2: 3 ^ 4
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_phenotypes_samples_not_rectangular()
        self.ok()

    def test_check_phenotypes_samples_not_rectangular_raises(self):
        config_string = \
'''phenotypes: a ^ b
samples:
    s1: 1
    s2: 1 ^ 2
    s3: 1 ^ 2 ^ 3
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               (r'Some \[samples\] had unexpected number of phenotype values '
                                r'\[expected 2 values\]: '
                                r'\(s1 \[1 values\], s3 \[3 values\]\); '
                                r'review config and try again.'),
                               validator._check_phenotypes_samples_not_rectangular)

    def test_check_comparison_missing_phenotype_value_singleValue(self):
        config_string = \
'''phenotypes: diet
samples:
    s1: HD
    s2: LD
    s3: ND
    s4: NA
comparisons:
    diet:
    - HD_v_LD
    - HD_v_ND'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)

        expected_message = (r'\(diet:NA\) are not present in \[comparisons\];'
                            r' some samples \(diet:s4\) will be excluded from comparisons .*')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._check_comparison_missing_phenotype_value)

    def test_check_comparison_missing_phenotype_value_muiltipleValues(self):
        config_string = \
'''phenotypes: pLabelA
samples:
    s1: pVal1
    s2: pVal2
    s3: pVal3
    s4: pVal4
    s5: pVal5
comparisons:
    pLabelA:
    - pVal1_v_pVal2
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)

        expected_message = (r'\(pLabelA:pVal3,pVal4,pVal5\)'
                            r' are not present in \[comparisons\]; '
                            r'some samples \(pLabelA:s3,s4,s5\) will be excluded from comparisons .*')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._check_comparison_missing_phenotype_value)

    def test_check_comparison_missing_phenotype_value_muiltiplePhenotypes(self):
        config_string = \
'''phenotypes: pLabelA ^ pLabelB
samples:
    s1: A1 ^ B4
    s2: A2 ^ B3
    s3: A3 ^ B2
    s4: A4 ^ B1
comparisons:
    pLabelA:
    - A1_v_A2
    pLabelB:
    - B1_v_B2
'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        
        expected_message = (r'\(pLabelA:A3,A4;pLabelB:B3,B4\)'
                            r' are not present in \[comparisons\]; '
                            r'some samples \(pLabelA:s3,s4;pLabelB:s1,s2\) will be excluded from comparisons .*')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._check_comparison_missing_phenotype_value)

    def test_check_comparison_missing_phenotype_label_ok(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
    s3: LD ^ M
    s4: HD ^ F
comparisons:
    diet:
        HD_v_LD
    gender:
    - M_v_F'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_comparison_missing_phenotype_label()
        self.ok()

    def test_check_comparison_missing_phenotype_label_singleValueRaises(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
    s3: LD ^ M
    s4: HD ^ F
comparisons:
    diet:
    - HD_v_LD'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        expected_message = (r'\(gender\) are not present in \[comparisons\].')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._check_comparison_missing_phenotype_label)

    def test_check_samples_excluded_from_comparison_ok(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
comparisons:
    diet:
    - HD_v_LD'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_samples_excluded_from_comparison()
        self.ok()

    def test_check_samples_excluded_from_comparison_missingBecausePhenoLabelRaises(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
    s3:    ^ M
    s4:    ^ M
comparisons:
    diet:
    - HD_v_LD'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        expected_message = r'Some samples \(s3,s4\) will not be compared.'
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._check_samples_excluded_from_comparison)

    def test_check_samples_excluded_from_comparison_missingBecausePhenoValueRaises(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ X
    s3:    ^ F
    s4:    ^ X
comparisons:
    gender:
    - M_v_F'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        expected_message = r'Some samples \(s2,s4\) will not be compared.'
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._check_samples_excluded_from_comparison)

    def test_check_comparison_references_unknown_phenotype_label_ok(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ X
comparisons:
    diet:
    - HD_v_LD'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_comparison_references_unknown_phenotype_label()
        self.ok()

    def test_check_comparison_references_unknown_phenotype_label_raises(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ X
    s3:    ^ F
    s4:    ^ X
comparisons:
    genotype:
    - M_v_F
    diet:
    - HD_v_LD
    status:
    - Hi_v_Lo'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        expected_message = r'referenced phenotype labels not found in.*genotype, status'
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               expected_message,
                               validator._check_comparison_references_unknown_phenotype_label)

    def test_check_comparison_references_unknown_phenotype_value_ok(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ X
comparisons:
    diet:
    - HD_v_LD'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        validator._check_comparison_references_unknown_phenotype_value()
        self.ok()

    def test_check_comparison_references_unknown_phenotype_value_raises(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ X
    s3:    ^ F
    s4:    ^ X
comparisons:
    diet:
    - HD_v_LD
    - HD_v_FOO
    gender:
    - M_v_F
    - X_v_BAZ
    - HD_v_BAR'''
        config = yaml.load(config_string)
        log = StringIO()
        validator = _ConfigValidator(config, log)
        expected_message = (r'referenced phenotype values not found in.*'
                            r'diet:FOO, gender:BAR, gender:BAZ')
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               expected_message,
                               validator._check_comparison_references_unknown_phenotype_value)

    def test_main_factors_malformed_ok(self):
        config_string = \
'''
main_factors: yes ^ no
phenotypes:   foo ^ bar
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        validator._check_main_factors_malformed()
        self.ok()

    def test_main_factors_malformed_wrongNumberValues(self):
        config_string = \
'''
main_factors: yes ^ no ^ yes
phenotypes:   foo ^ bar
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        expected_message = (r'Number of values .* must match \(3!=2\)')
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               expected_message,
                               validator._check_main_factors_malformed)

    def test_main_factors_illegal_values_ok(self):
        config_string = \
'''
main_factors: yes ^ no
phenotypes:   foo ^ bar
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        validator._check_main_factors_illegal_values()
        self.ok()

    def test_main_factors_blank_values(self):
        config_string = \
'''
main_factors:     ^ nope
phenotypes:   foo ^ bar
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        expected_message = (r'Expected \[main_factor\] values .* found \(<blank>, nope\)')
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               expected_message,
                               validator._check_main_factors_illegal_values)

    def test_main_factors_illegal_values(self):
        config_string = \
'''
main_factors: yep ^ nope
phenotypes:   foo ^ bar
'''
        validator = _ConfigValidator(yaml.load(config_string), StringIO())
        expected_message = (r'Expected \[main_factor\] values .* found \(nope, yep\)')
        self.assertRaisesRegex(_WatermelonConfigFailure,
                               expected_message,
                               validator._check_main_factors_illegal_values)



    def test_init_validationsAreIncluded(self):
        validator = _ConfigValidator(config={}, log=StringIO())
        primary_validation_names = [f.__name__ for f in validator._PRIMARY_VALIDATIONS]
        self.assertEqual(['_check_missing_required_field',
                          '_check_samples_malformed',
                          '_check_samples_not_stringlike',
                          '_check_comparisons_malformed',
                          '_check_comparisons_not_a_pair',
                          '_check_phenotypes_samples_not_rectangular'],
                         primary_validation_names)
        secondary_validation_names = [f.__name__ for f in validator._SECONDARY_VALIDATIONS]
        self.assertEqual(['_check_comparison_missing_phenotype_value',
                          '_check_comparison_missing_phenotype_label',
                          '_check_comparison_references_unknown_phenotype_label',
                          '_check_comparison_references_unknown_phenotype_value',
                          '_check_samples_excluded_from_comparison'],
                         secondary_validation_names)

class ValidationCollectorTest(unittest.TestCase):
    def test_passed_true(self):
        collector = _ValidationCollector(None)
        self.assertEqual(True, collector.passed)
        passing_validation = MockValidation()
        collector.check(passing_validation.validate)
        self.assertEqual(True, collector.passed)
        self.assertEqual([], collector.warnings)
        self.assertEqual([], collector.failures)

    def test_passed_falseIfConfigWarning(self):
        collector = _ValidationCollector(None)
        self.assertEqual(True, collector.passed)
        problem = _WatermelonConfigWarning('foo')
        problem_validation = MockValidation(problem)
        collector.check(problem_validation.validate)
        self.assertEqual(False, collector.passed)
        self.assertEqual([problem], collector.warnings)
        self.assertEqual([], collector.failures)

    def test_passed_falseIfConfigError(self):
        collector = _ValidationCollector(None)
        self.assertEqual(True, collector.passed)
        problem = _WatermelonConfigFailure('error')
        problem_validation = MockValidation(problem)
        collector.check(problem_validation.validate)
        self.assertEqual(False, collector.passed)
        self.assertEqual([], collector.warnings)
        self.assertEqual([problem], collector.failures)

    def test_passed_raiseIfUnknownError(self):
        collector = _ValidationCollector(None)
        self.assertEqual(True, collector.passed)
        problem_validation = MockValidation(ValueError('foo'))
        self.assertRaisesRegex(ValueError,
                               'foo',
                               collector.check,
                               problem_validation.validate)

    def test_log_validation_results_ok(self):
        log = StringIO()
        collector = _ValidationCollector(log)

        collector.log_results()

        lines = iter(log.getvalue().split('\n'))
        self.assertRegex(next(lines), r'===')
        self.assertEqual('config validation: OK', next(lines))

    def test_log_validation_results_warning(self):
        log = StringIO()
        collector = _ValidationCollector(log)
        validation = MockValidation(_WatermelonConfigWarning('warn1'))
        collector.check(validation.validate)
        collector.log_results()

        lines = iter(log.getvalue().split('\n'))
        self.assertRegex(next(lines), r'===')
        self.assertEqual('config validation: WARNING (1 warnings):', next(lines))
        self.assertRegex(next(lines), r'---')
        self.assertEqual('warning: warn1', next(lines))

    def test_log_validation_results_failure(self):
        log = StringIO()
        collector = _ValidationCollector(log)
        validation = MockValidation(_WatermelonConfigFailure('fail1'))
        collector.check(validation.validate)
        collector.log_results()

        lines = iter(log.getvalue().split('\n'))
        self.assertRegex(next(lines), r'===')
        self.assertEqual('config validation: FAILED (1 failures):', next(lines))
        self.assertRegex(next(lines), r'---')
        self.assertEqual('failure: fail1', next(lines))

    def test_log_validation_failuresComeFirst(self):
        log = StringIO()
        collector = _ValidationCollector(log)
        validation1 = MockValidation(_WatermelonConfigFailure('fail1'))
        validation2 = MockValidation(_WatermelonConfigWarning('warn1'))
        validation3 = MockValidation(_WatermelonConfigFailure('fail2'))
        validation4 = MockValidation(_WatermelonConfigWarning('warn2'))
        collector.check(validation1.validate)
        collector.check(validation2.validate)
        collector.check(validation3.validate)
        collector.check(validation4.validate)

        collector.log_results()

        lines = iter(log.getvalue().split('\n'))
        self.assertRegex(next(lines), r'===')
        self.assertEqual('config validation: FAILED (2 failures, 2 warnings):', next(lines))
        self.assertRegex(next(lines), r'---')
        self.assertEqual('failure: fail1', next(lines))
        self.assertEqual('failure: fail2', next(lines))
        self.assertEqual('warning: warn1', next(lines))
        self.assertEqual('warning: warn2', next(lines))

    def test_ok_to_proceed_TrueIfPassed(self):
        log = StringIO()
        prompt = MockPromptOverride(False)
        collector = _ValidationCollector(log)
        self.assertEqual(True, collector.ok_to_proceed(prompt.prompt_to_override))
        self.assertEqual(False, prompt.prompt_to_override_was_called)
        self.assertRegex(log.getvalue(), r'===')

    def test_ok_to_proceed_FalseIfFailed(self):
        log = StringIO()
        prompt = MockPromptOverride(False)
        collector = _ValidationCollector(log)
        collector.failures = [1]
        self.assertEqual(False, collector.ok_to_proceed(prompt.prompt_to_override))
        self.assertEqual(False, prompt.prompt_to_override_was_called)
        lines = iter(log.getvalue().split('\n'))
        self.assertRegex(next(lines), r'===')
        self.assertRegex(next(lines), 'There were.*failures')
        self.assertRegex(next(lines), 'Watermelon stopped')

    def test_ok_to_proceed_TrueIfWarningsWhenOverride(self):
        log = StringIO()
        prompt = MockPromptOverride(True)
        collector = _ValidationCollector(log)
        collector.warnings = [1]
        self.assertEqual(True, collector.ok_to_proceed(prompt.prompt_to_override))
        self.assertEqual(True, prompt.prompt_to_override_was_called)
        lines = iter(log.getvalue().split('\n'))
        self.assertRegex(next(lines), r'===')
        self.assertRegex(next(lines), 'There were.*warnings')
        self.assertRegex(next(lines), 'Watermelon will continue despite warnings above.')

    def test_ok_to_proceed_FalseIfWarningsWhenNoOverride(self):
        log = StringIO()
        prompt = MockPromptOverride(False)
        collector = _ValidationCollector(log)
        collector.warnings = [1]
        self.assertEqual(False, collector.ok_to_proceed(prompt.prompt_to_override))
        self.assertEqual(True, prompt.prompt_to_override_was_called)
        lines = iter(log.getvalue().split('\n'))
        self.assertRegex(next(lines), r'===')
        self.assertRegex(next(lines), 'There were.*warnings')
        self.assertRegex(next(lines), 'Watermelon stopped')

    def test_main_configOk(self):
        mock_override = MockPromptOverride(False)
        log = StringIO()
        config_contents = \
'''main_factors: yes ^ no
phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
comparisons:
    diet:
    - HD_v_LD
    gender:
    - M_v_F
'''
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_filename = _write_config_file(temp_dir.path,
                                                 'config.yaml',
                                                 config_contents)

            exit_code = main(config_filename,
                             log=log,
                             prompt_to_override=mock_override.prompt_to_override)

        lines = iter(log.getvalue().split('\n'))
        self.assertEquals(0, exit_code)
        self.assertEquals(False, mock_override.prompt_to_override_was_called)
        self.assertRegex(next(lines), r'===')
        self.assertEquals('config validation: OK', next(lines))

    def test_main_invalidYaml(self):
        mock_override = MockPromptOverride(False)
        log = StringIO()
        config_contents = '1:\n2\n3\n'
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_filename = _write_config_file(temp_dir.path,
                                                 'config.yaml',
                                                 config_contents)

            exit_code = main(config_filename,
                             log=log,
                             prompt_to_override=mock_override.prompt_to_override)

        lines = iter(log.getvalue().split('\n'))
        self.assertEquals(1, exit_code)
        self.assertEquals(False, mock_override.prompt_to_override_was_called)
        self.assertRegex(next(lines), r'===')
        self.assertEqual('config validation: ERROR', next(lines))
        self.assertEqual(('Could not parse this config file [{}]; '
                           'is it valid YAML?').format(config_filename),
                          next(lines))


    def test_main_configWarningStop(self):
        mock_override = MockPromptOverride(False)
        log = StringIO()
        config_contents = \
'''main_factors: yes ^ no
phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
comparisons:
    diet:
    - HD_v_LD
'''
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_filename = _write_config_file(temp_dir.path,
                                                 'config.yaml',
                                                 config_contents)

            exit_code = main(config_filename,
                             log=log,
                             prompt_to_override=mock_override.prompt_to_override)

        lines = iter(log.getvalue().split('\n'))
        self.assertEqual(1, exit_code)
        self.assertEqual(True, mock_override.prompt_to_override_was_called)
        self.assertRegex(next(lines), r'===')
        self.assertRegex(next(lines), r'config validation: WARNING')

    def test_main_configWarningContinue(self):
        mock_override = MockPromptOverride(True)
        log = StringIO()
        config_contents = \
'''main_factors: yes ^ no
phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
comparisons:
    diet:
    - HD_v_LD
'''
        with TempDirectory() as temp_dir:
            temp_dir_path = temp_dir.path
            config_filename = _write_config_file(temp_dir.path,
                                                 'config.yaml',
                                                 config_contents)

            exit_code = main(config_filename,
                             log=log,
                             prompt_to_override=mock_override.prompt_to_override)

        lines = iter(log.getvalue().split('\n'))
        self.assertEqual(0, exit_code)
        self.assertEqual(True, mock_override.prompt_to_override_was_called)
        self.assertRegex(next(lines), r'===')
        self.assertRegex(next(lines), r'config validation: WARNING')
