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

import argparse
import os
import numpy as np
import pandas as pd

DESCRIPTION = \
'''Calculates and compiles read and insert statistics from bbmap output files on a per sample basis.'''

#set definitions

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--sample_id',
        type=str,
        help='Sample id, from which the Sample_ihist.txt and Sample_lhist.txt file '
              'names will be derived. This id will also be used to identify the output.',
        required=True)

    parser.add_argument(
        '--output_dir',
        type=str,
        help='Path to dir to which output files will be written. ',
        required=True)

    parser.add_argument(
        '--input_dir',
        type=str,
        help='Path to dir containing input files, Sample_ihist.txt and Sample_lhist.txt. ',
        required=True)

    args = parser.parse_args(sys_argv)

    return args

_LHIST_LENGTH = 'Length'
_LHIST_COUNT = 'Count'

def _calc_mean(read_data):
    read_data['Tot'] = read_data[_LHIST_LENGTH]*read_data[_LHIST_COUNT]
    return read_data['Tot'].sum() / read_data[_LHIST_COUNT].sum()

def _calc_inner(insert_data, m):
    inner_mate_dist = ins_mean - (2*m)
    return(inner_mate_dist)

def _build_read_stats(read_data_df, insert_data_df):
    #define insert mean and sd
    ins_mean = round(float(insert_data_df[1][0]))
    ins_sd = round(float(insert_data_df[1][3]))

    #calc_mean
    m = _calc_mean(read_data_df)

    #calc_inner_mate
    im = _calc_inner(insert_data_df, m)

    #collect information into dictionary, format to dataframe
    d = {'#sample': args.sample_id,'insert_mean': ins_mean,'insert_std_dev':ins_sd ,'read_mean':m,'inner_mate_dist':im}
    df = pd.DataFrame(list(d.items())).transpose() #convert dictionary to dataframe, transpose
    df.columns = df.iloc[0] #create header row from index 0 values
    df = df.reindex(df.index.drop(0))#drop index 0
    return df


def main(sys_argv):
    #args
    args = _parse_command_line_args(sys_argv)

    #define names and paths
    in_file_name_1 = args.sample_id + '_lhist.txt'
    in_file_name_2 = args.sample_id + '_ihist.txt'
    out_file_name = args.sample_id + '_read_stats.txt'
    in_path = args.input_dir
    out_path = args.output_dir

    #read in files
    read_data = pd.read_csv(os.path.join(in_path,in_file_name_1), sep='\t', header = 0)
    read_data.columns = list(map(lambda x: x.lstrip('#'), read_data.columns.values))

    insert_data = pd.read_csv(os.path.join(in_path,in_file_name_2), sep='\t', header = None)

    read_stats_df = _build_read_stats(read_data, insert_data)

    read_stats_df.to_csv(os.path.join(out_path,out_file_name), sep='\t', encoding='utf-8', header = True, index = False)

if __name__ == '__main__':
    main(sys.argv[1:])
