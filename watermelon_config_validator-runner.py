#!/usr/bin/env python
"""Convenience wrapper for running config_validator directly from source tree."""

import os
import sys

from scripts.config_validator import main

if __name__ == '__main__':
    config_filepath = os.path.realpath(sys.argv[1])
    schema_filepath = os.path.realpath(sys.argv[2])
    exit(main(config_filepath, schema_filepath))
