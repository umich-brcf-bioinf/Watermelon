#Usage: python.script.py(--sample_id sample_id, --input_dir /path/to/input, --output_dir /path/to/output)

'''Calculates and compiles read and insert statistics from bbmap output files.

Specifically this does three things:
1) read_stats_bbmap accepts an input directory containing the bbmap output files.
   These files will have the {Sample}'_lhist.txt' and {Sample}'_ihist.txt' naming convention.
   '{Sample}_lhist.txt' contains read-lenght distribution in as text histogram format.
   '{Sample}_ihist.txt' contains insert statistics (mean, s.d., mode, etc.) and distribution.

2) read_stats_bbmap accepts an output directory, to which it will write a tab-delimited
   file (for each file) containing the calculated statistics in a '{Sample}_read_stats.txt' naming convention.
   '{Sample}_read_stats.txt' contains the following fields:
   a) sample_id = sample id based upon the '{Sample}_lhist.txt' file
   b) insert_mean = mean insert value
   c) insert_std_dev = insert standard deviation value
   d) read_mean = mean read length value
   e) inner_mate_dist = inner-mate distance value

3) read_stats_bbmap accepts a sample_id, for which the calculations will be run.
'''
import os
import sys
import argparse
import csv
import numpy as np
import pandas as pd

DESCRIPTION = \
'''Calculates and compiles read and insert statistics from bbmap output files on a per sample basis.'''

_IHIST_LABELS = argparse.Namespace(mean='#Mean', sd='#STDev')
_IHIST_REQUIRED_LABELS = set([_IHIST_LABELS.mean, _IHIST_LABELS.sd])

_LHIST_FIELDS = argparse.Namespace(length='Length', count='Count')

#set definitions

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--sample_id',
        type=str,
        help=('Sample id, from which the {sample_id}_ihist.txt and '
              '{sample_id}_lhist.txt file names will be derived. This id will '
              'also be used to identify the output.'),
        required=True)

    parser.add_argument(
        '--output_dir',
        type=str,
        help='Path to dir to which output files will be written. ',
        required=True)

    parser.add_argument(
        '--input_dir',
        type=str,
        help=('Path to dir containing input files, Sample_ihist.txt and '
              'Sample_lhist.txt. '),
        required=True)

    args = parser.parse_args(sys_argv)

    return args

def _calc_mean(read_data):
    read_data['Tot'] = read_data[_LHIST_FIELDS.length]*read_data[_LHIST_FIELDS.count]
    return read_data['Tot'].sum() / read_data[_LHIST_FIELDS.count].sum()

def _get_mean_stdev_from_ihist(filename, insert_data):
    def is_float(x):
        try:
            float(x.split()[1])
            return True
        except:
            return False

    def is_valid_pair(x):
        return x.startswith('#') and len(x.split()) == 2 and is_float(x)

    header_values = dict([x.split() for x in insert_data if is_valid_pair(x)])

    missing_labels = _IHIST_REQUIRED_LABELS - set(header_values.keys())
    if missing_labels:
        msg_fmt = 'File [{}] is missing label or invalid value for [{}]'
        msg = msg_fmt.format(filename, ', '.join(sorted(missing_labels)))
        raise ValueError(msg)

    return float(header_values[_IHIST_LABELS.mean]), float(header_values[_IHIST_LABELS.sd])

def _build_read_stats(sample_id, ins_mean, ins_sd, read_data_df):
    #calc_mean
    m = _calc_mean(read_data_df)

    #calc_inner_mate
    im = ins_mean - (2 * m)

    d = [['#sample', 'insert_mean', 'insert_std_dev','read_mean', 'inner_mate_dist'],
        [sample_id, ins_mean, ins_sd, m, im]]
    return d

def main(sys_argv):
    #args
    args = _parse_command_line_args(sys_argv)

    #define names and paths
    lhist_file_name = os.path.join(args.input_dir, args.sample_id + '_lhist.txt')
    ihist_file_name = os.path.join(args.input_dir,
                                   args.sample_id + '_ihist.txt')
    out_file_name = os.path.join(args.input_dir,
                                 args.sample_id + '_read_stats.txt')
    # in_path = args.input_dir
    # out_path = args.output_dir

    #read in files
    with open(lhist_file_name) as lhist_file:
        read_data = pd.read_csv(lhist_file,
                                sep='\t',
                                header=0)
    read_data.columns = list(map(lambda x: x.lstrip('#'),
                             read_data.columns.values))

    with open(ihist_file_name, 'r') as insert_file:
        (ins_mean, ins_sd) = _get_mean_stdev_from_ihist(ihist_file_name,
                                                        insert_file)
    read_stats_df = _build_read_stats(args.sample_id, ins_mean, ins_sd, read_data)

    #write out file
    with open(out_file_name, "w") as output:
        wr = csv.writer(output, lineterminator='\n', delimiter='\t')
        wr.writerows(read_stats_df)


if __name__ == '__main__':
    main(sys.argv[1:])
