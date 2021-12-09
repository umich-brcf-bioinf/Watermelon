#!/usr/bin/env python3

import argparse
import re
import warnings


def _parse_lines(lines):
    '''Parse the lines of rseqc's split_bam.py output
    Return a dictionary with simplified keys 'in', 'ex', 'junk' and integer vals'''
    summary = {}
    for line in lines:
        dkey = line.split(":")[0].strip()
        bamtest = re.match(r".*\.(.*)\.bam .*", dkey) #Capture which type of bam
        if bamtest:
            dkey = bamtest.group(1) #Use bam type as dict key
        dval = int(line.split(":")[1].strip())
        summary[dkey] = dval
    return summary


def _calculate_percentages(data_dict):
    '''Takes the counts derived from '.in.bam', '.ex.bam', and '.junk.bam'
    Returns a tuple of percentages in that order
    '''
    try:
        p_in = round(data_dict['in'] / data_dict['Total records'] * 100, 2)
        p_ex = round(data_dict['ex'] / data_dict['Total records'] * 100, 2)
        p_junk = round(data_dict['junk'] / data_dict['Total records'] * 100, 2)
    except ZeroDivisionError:
        msg = "Warning: Zero value prevents division - Returning 0 % for all calculations."
        warnings.warn(msg)
        (p_in, p_ex, p_junk) = (0,0,0)
    return (p_in, p_ex, p_junk)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="split_bam_percentage.py parses results from rseqc's split_bam.py and returns the percentages in each bam")
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Path to input file - results from rseqc's split_bam.py")
    parser.add_argument("-n", "--samplename", type=str, required=True, help="Sample name to report in-line with results")
    parser.add_argument("-l", "--labels", nargs=3, default=["Incl", "Excl", "QCFail/Unmapped"], help="Labels for printing the results. They correspond to '.in.bam', '.ex.bam', '.junk.bam', in that order.")
    args = parser.parse_args()

    with open(args.input_file, "r") as ifh:
        info = _parse_lines(ifh.readlines())
        percentages = _calculate_percentages(info)
        labels = ["Samplename"] + args.labels
        values = [args.samplename] + [str(x) for x in percentages]
        print("\t".join(labels))
        print("\t".join(values))
