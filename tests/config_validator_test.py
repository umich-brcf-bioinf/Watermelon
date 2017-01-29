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

from scripts.rnaseq_snakefile_helper import PhenotypeManager
from scripts.config_validator import _ConfigValidator
from scripts.config_validator import _WatermelonConfigFailure
from scripts.config_validator import _WatermelonConfigWarning
from scripts.config_validator import _ValidationCollector
from scripts.config_validator import main

class MockValidation(object):
    def __init__(self, warning_or_error=None):
        self.validate_was_called = False
        self.warning_or_error = warning_or_error

    def validate(self):
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

    def test_validate_ok(self):
        validator = _ConfigValidator({})
        validation1 = MockValidation()
        validation2 = MockValidation()
        validator._validations=[validation1.validate, validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(True, validation_result.passed)
        self.assertEqual([], validation_result.failures)
        self.assertEqual([], validation_result.warnings)

    def test_validate_warning(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        validation1 = MockValidation()
        warning =_WatermelonConfigWarning("I'm angry")
        validation2 = MockValidation(warning)
        validator._validations=[validation1.validate, validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([warning], validation_result.warnings)
        self.assertEqual([], validation_result.failures)


    def test_validate_earlyWarningContinuesProcessingValidations(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        warning =_WatermelonConfigWarning("I'm angry")
        validation1 = MockValidation(warning)
        validation2 = MockValidation()
        validator._validations=[validation1.validate, validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([warning], validation_result.warnings)
        self.assertEqual([], validation_result.failures)

    def test_validate_error(self):
        log = StringIO()
        validator = _ConfigValidator({}, log)
        validation1 = MockValidation()
        error =_WatermelonConfigFailure("I'm angry")
        validation2 = MockValidation(error)
        validator._validations=[validation1.validate, validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([], validation_result.warnings)
        self.assertEqual([error], validation_result.failures)


    def test_validate_mixOfWarningsAndErrors(self):
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
        validator._validations=[validation1.validate,
                                validation2.validate,
                                validation3.validate,
                                validation4.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([problem2, problem3], validation_result.warnings)
        self.assertEqual([problem1, problem4], validation_result.failures)



    def test_comparison_missing_phenotype_value_singleValue(self):
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
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        
        expected_message = (r'\(diet:NA\) are not present in comparisons;'
                            r' some samples \(diet:s4\) will be excluded from comparisons .*')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._comparison_missing_phenotype_value)

    def test_comparison_missing_phenotype_value_muiltipleValues(self):
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
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        
        expected_message = (r'\(pLabelA:pVal3,pVal4,pVal5\)'
                            r' are not present in comparisons; '
                            r'some samples \(pLabelA:s3,s4,s5\) will be excluded from comparisons .*')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._comparison_missing_phenotype_value)

    def test_comparison_missing_phenotype_value_muiltiplePhenotypes(self):
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
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        
        expected_message = (r'\(pLabelA:A3,A4;pLabelB:B3,B4\)'
                            r' are not present in comparisons; '
                            r'some samples \(pLabelA:s3,s4;pLabelB:s1,s2\) will be excluded from comparisons .*')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._comparison_missing_phenotype_value)

    def test_comparison_missing_phenotype_label_ok(self):
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
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        validator._comparison_missing_phenotype_label()
        self.ok()

    def test_comparison_missing_phenotype_label_singleValueRaises(self):
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
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        expected_message = (r'\(gender\) are not present in comparisons.')
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._comparison_missing_phenotype_label)

    def test_samples_excluded_from_comparison_ok(self):
        config_string = \
'''phenotypes: diet ^ gender
samples:
    s1: HD ^ M
    s2: LD ^ F
comparisons:
    diet:
    - HD_v_LD'''
        config = yaml.load(config_string)
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        validator._samples_excluded_from_comparison()
        self.ok()

    def test_samples_excluded_from_comparison_missingBecausePhenoLabelRaises(self):
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
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        expected_message = r'Some samples \(s3,s4\) will not be compared.'
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._samples_excluded_from_comparison)

    def test_samples_excluded_from_comparison_missingBecausePhenoValueRaises(self):
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
        phenotype_manager = PhenotypeManager(config)
        log = StringIO()
        validator = _ConfigValidator(phenotype_manager, log)
        expected_message = r'Some samples \(s2,s4\) will not be compared.'
        self.assertRaisesRegex(_WatermelonConfigWarning,
                               expected_message,
                               validator._samples_excluded_from_comparison)


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
        
        self.assertEqual('config validation: OK\n', log.getvalue())

    def test_log_validation_results_warning(self):
        log = StringIO()
        collector = _ValidationCollector(log)
        validation = MockValidation(_WatermelonConfigWarning('warn1'))
        collector.check(validation.validate)
        collector.log_results()

        lines = iter(log.getvalue().split('\n'))
        self.assertEqual('config validation: WARNING (1 warnings):', next(lines))
        self.assertEqual('warning: warn1', next(lines))

    def test_log_validation_results_failure(self):
        log = StringIO()
        collector = _ValidationCollector(log)
        validation = MockValidation(_WatermelonConfigFailure('fail1'))
        collector.check(validation.validate)
        collector.log_results()

        lines = iter(log.getvalue().split('\n'))
        self.assertEqual('config validation: FAILED (1 failures):', next(lines))
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
        self.assertEqual('config validation: FAILED (2 failures, 2 warnings):', next(lines))
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
        self.assertEqual('', log.getvalue())

    def test_ok_to_proceed_FalseIfFailed(self):
        log = StringIO()
        prompt = MockPromptOverride(False)
        collector = _ValidationCollector(log)
        collector.failures = [1]
        self.assertEqual(False, collector.ok_to_proceed(prompt.prompt_to_override))
        self.assertEqual(False, prompt.prompt_to_override_was_called)
        lines = iter(log.getvalue().split('\n'))
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
        self.assertRegex(next(lines), 'There were.*warnings')
        self.assertRegex(next(lines), 'Watermelon stopped')

    def test_main_configOk(self):
        mock_override = MockPromptOverride(False)
        log = StringIO()
        config_contents = \
'''phenotypes: diet ^ gender
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
        self.assertEquals('config validation: OK', next(lines))

    def test_main_configWarningStop(self):
        mock_override = MockPromptOverride(False)
        log = StringIO()
        config_contents = \
'''phenotypes: diet ^ gender
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
        self.assertEquals(1, exit_code)
        self.assertEquals(True, mock_override.prompt_to_override_was_called)
        self.assertRegex(next(lines), r'config validation: WARNING')

    def test_main_configWarningContinue(self):
        mock_override = MockPromptOverride(True)
        log = StringIO()
        config_contents = \
'''phenotypes: diet ^ gender
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
        self.assertEquals(0, exit_code)
        self.assertEquals(True, mock_override.prompt_to_override_was_called)
        self.assertRegex(next(lines), r'config validation: WARNING')
