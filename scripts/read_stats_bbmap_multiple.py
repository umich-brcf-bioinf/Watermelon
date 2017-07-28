#Usage: python.script.py(--input_dir /path/to/input, --output_dir /path/to/output)

'''Calculates and compiles read and insert statistics from bbmap output files.

Specifically this does three things:
1) read_stats_bbmap accepts a source directory containing the bbmap output files.
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

'''

import argparse
import os
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

DESCRIPTION = \
'''Calculates and compiles read and insert statistics from all bbmap output files in a directory.'''

#set definitions

def _parse_command_line_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
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

    args = parser.parse_args()

    if not (args.output_dir or args.input_dir):
        parser.error('You must specify the sample_id, the output_dir, and input_dir')
    return args

def calc_mean(read_data):
    read_data['Tot'] = read_data['Length']*read_data['Count']
    col_list = list(read_data)
    col_sums = read_data[col_list[1:]].sum(axis=0)
    read_mean = col_sums[1]/col_sums[0]
    return(read_mean)

def calc_inner(insert_data, m):
    inner_mate_dist = ins_mean - (2*m)
    return(inner_mate_dist)

#args
args = _parse_command_line_args()

#define names and paths
in_path = args.input_dir
out_path = args.output_dir

#read list of files in directory
onlyfiles = [f for f in listdir(in_path) if isfile(join(in_path, f))]

#get files with _ihist.txt or _lhist.txt
ihist_list = [name for name in onlyfiles if '_ihist' in name]
lhist_list = [name for name in onlyfiles if '_lhist' in name]

#copy list and parse copy list to remove _ihist.txt and lhist.txt for writing out files
sample_list = ihist_list
sample_list = [sample.replace('_ihist.txt', '') for sample in sample_list]

#read in files
for i in range(0,len(sample_list)):
        #read in files
        read_data = pd.read_csv(os.path.join(in_path,lhist_list[i]), sep='\t', header = 0, names = ['Length','Count'])
        insert_data = pd.read_csv(os.path.join(in_path,ihist_list[i]), sep='\t', header = None)

        #define insert mean and sd
        ins_mean = round(float(insert_data[1][0]))
        ins_sd = round(float(insert_data[1][3]))

        #calc_mean
        m = calc_mean(read_data)

        #calc_inner_mate
        im = calc_inner(insert_data, m)

        #collect information into dictionary, format to dataframe
        d = {'#sample': sample_list[i],'insert_mean': ins_mean,'insert_std_dev':ins_sd ,'read_mean':m,'inner_mate_dist':im}
        df = pd.DataFrame(list(d.items())).transpose() #convert dictionary to dataframe, transpose
        df.columns = df.iloc[0] #create header row from index 0 values
        df = df.reindex(df.index.drop(0))#drop index 0

        #write out
        out_file_name = sample_list[i] + '_read_stats.txt'
        df.to_csv(os.path.join(out_path,out_file_name), sep='\t', encoding='utf-8', header = True, index = False)
