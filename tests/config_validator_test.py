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

import yaml

from scripts.rnaseq_snakefile_helper import  PhenotypeManager
from scripts.config_validator import ConfigValidator
from scripts.config_validator import WatermelonConfigWarning

class ConfigValidatorTest(unittest.TestCase):
    class MockValidation(object):
        def __init__(self, warning_or_error=None):
            self.validate_was_called = False
            self.warning_or_error = warning_or_error

        def validate(self):
            self.validate_was_called = True
            if self.warning_or_error:
                raise self.warning_or_error

    def ok(self):
        self.assertEqual(1, 1)

    def test_validate_ok(self):
        validator = ConfigValidator({})
        validation1 = ConfigValidatorTest.MockValidation()
        validation2 = ConfigValidatorTest.MockValidation()
        validator._validations=[validation1.validate, validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(True, validation_result.passed)
        self.assertEqual([], validation_result.errors)
        self.assertEqual([], validation_result.warnings)

    def test_validate_error(self):
        stderr = StringIO()
        validator = ConfigValidator({}, stderr=stderr)
        validation1 = ConfigValidatorTest.MockValidation()
        warning = WatermelonConfigWarning("I'm angry")
        validation2 = ConfigValidatorTest.MockValidation(warning)
        validator._validations=[validation1.validate, validation2.validate]

        validation_result = validator.validate()

        self.assertEqual(True, validation1.validate_was_called)
        self.assertEqual(True, validation2.validate_was_called)
        self.assertEqual(False, validation_result.passed)
        self.assertEqual([warning], validation_result.warnings)
        self.assertEqual([], validation_result.errors)


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
        validator = ConfigValidator(phenotype_manager)
        
        expected_message = (r'\(diet:NA\) are not present in comparisons;'
                            r' some samples \(diet:s4\) will be excluded from comparisons .*')
        self.assertRaisesRegex(WatermelonConfigWarning,
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
        validator = ConfigValidator(phenotype_manager)
        
        expected_message = (r'\(pLabelA:pVal3,pVal4,pVal5\)'
                            r' are not present in comparisons; '
                            r'some samples \(pLabelA:s3,s4,s5\) will be excluded from comparisons .*')
        self.assertRaisesRegex(WatermelonConfigWarning,
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
        validator = ConfigValidator(phenotype_manager)
        
        expected_message = (r'\(pLabelA:A3,A4;pLabelB:B3,B4\)'
                            r' are not present in comparisons; '
                            r'some samples \(pLabelA:s3,s4;pLabelB:s1,s2\) will be excluded from comparisons .*')
        self.assertRaisesRegex(WatermelonConfigWarning,
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
        validator = ConfigValidator(phenotype_manager)
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
        validator = ConfigValidator(phenotype_manager)
        expected_message = (r'\(gender\) are not present in comparisons.')
        self.assertRaisesRegex(WatermelonConfigWarning,
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
        validator = ConfigValidator(phenotype_manager)
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
        validator = ConfigValidator(phenotype_manager)
        expected_message = r'Some samples \(s3,s4\) will not be compared.'
        self.assertRaisesRegex(WatermelonConfigWarning,
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
        validator = ConfigValidator(phenotype_manager)
        expected_message = r'Some samples \(s2,s4\) will not be compared.'
        self.assertRaisesRegex(WatermelonConfigWarning,
                               expected_message,
                               validator._samples_excluded_from_comparison)