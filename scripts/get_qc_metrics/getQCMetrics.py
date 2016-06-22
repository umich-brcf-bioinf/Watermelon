#!/usr/bin/env python
#Intended to replace getQCMetrics.pl in the RNA-Seq pipeline.
#Jessica Bene (jebene) 03/2016
from __future__ import print_function, absolute_import

import argparse
import os
import re
import sys

STEP_NAME = "align"

def get_arguments(args):
    parser = argparse.ArgumentParser(description="obtain alignment rates for a given sample")
    parser.add_argument("align_dir", help="path to alignment directory")
    parser.add_argument("sample", help="sample name")
    parser.add_argument("output_file", help="path to output file")

    return parser.parse_args(args)

def check_direction(line, left, right):
    if line.startswith("Left reads"):
        left = True
        right = False
    elif line.startswith("Right reads"):
        left = False
        right = True
    elif line.startswith("Aligned pairs"):
        left = False
        right = False
    return left, right

def get_key(line, left, right):
    if left or right:
        if "Input" in line:
            return "{}_input"
        elif "Mapped" in line:
            return "{}_mapped"
        elif "of these" in line:
            return "{}_multiple"
        else:
            return False
    return False

def form_line(sample, left, right, key, line):
    direction = ""

    if left:
        direction = "left"
    elif right:
        direction = "right"

    if direction and key:
        value = line.split(":")[1].strip().split()[0]
        return "\t".join([sample, STEP_NAME, key.format(direction), value]) + "\n"
    return False

def get_alignment_rate(sample, align_summary):
    left = False
    right = False
    lines = []

    for line in align_summary:
        left, right = check_direction(line, left, right)

        key = get_key(line, left, right)
        line = form_line(sample, left, right, key, line)

        if line:
            lines.append(line)

    return lines

def main():
    args = get_arguments(sys.argv[1:])
    align_summary_path = os.path.join(args.align_dir, args.sample, "align_summary.txt")
    if os.path.isfile(align_summary_path):
        align_summary = open(align_summary_path, "r")
        lines = get_alignment_rate(args.sample, align_summary)
        align_summary.close()

    output_file = open(args.output_file, "w")
    for line in lines:
        output_file.write(line)
    output_file.close()

if __name__ == "__main__":
    main()
