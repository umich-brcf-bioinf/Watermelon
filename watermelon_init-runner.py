#!/usr/bin/env python

"""Convenience wrapper for running watermelon_init directly from source tree."""

from scripts.watermelon_init import main
import sys

if __name__ == '__main__':
    main(sys.argv[1:])
