'''Combines the count/FPKM/TPM from individual sample outputs into one matrix'''
import argparse
from glob import glob
from os.path import commonprefix, commonpath
import re
import sys

import pandas as pd

__version__ = '0.0.1'
_DESCRIPTION = \
'''Accepts tab-separated sample isoform files and combines into single tab-separated matrix.'''

def _commonsuffix(strings):
    return commonprefix(list(map(lambda s:s[::-1], strings)))[::-1]

def _build_sample_files(file_glob):
    sample_files = []
    file_names = glob(file_glob)
    prefix = commonpath(file_names)
    suffix = _commonsuffix(file_names)
    for file_name in sorted(file_names):
        sample_name = file_name.lstrip(prefix).rstrip(suffix)
        sample_files.append((sample_name, file_name))
    return sample_files

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(
        description=_DESCRIPTION)
    parser.add_argument(
        '-o', '--output_file',
        type=str,
        help='path to combined output file',
        required=True)
    parser.add_argument(
        '-i', '--input_path',
        type=str,
        help='path (including linux wildcards) to sample input files; surround with single quotes when usimg wildcards',
        required=True)
    parser.add_argument(
        '-c', '--column',
        type=str,
        help='full name of column to extract from inputs (e.g. FPKM)',
        required=True)
    parser.add_argument(
        '--id_columns',
        type=str,
        help='gene_id or gene_id,transcript_id',
        required=True)

    parser.add_argument('--version',
                    '-V',
                    action='version',
                    version=__version__)
    args = parser.parse_args(sys_argv)
    args.id_columns=args.id_columns.split(',')
    return args


def main(argv):
    print('combine v{}'.format(__version__))
    print('command line args: {}'.format(' '.join(argv)))
    args = _parse_command_line_args(argv[1:])

    sample_files = _build_sample_files(args.input_path)
    output_filename = args.output_file
    merge_column = args.column

    name, file = sample_files.pop(0)
    name = name + '|' + merge_column
    df=pd.read_csv(file, sep='\t',low_memory=False)

    new=pd.DataFrame(df[args.id_columns+[merge_column]])
    new.rename(columns={merge_column:name},inplace=True)

    for (name, file) in sample_files:
        df=pd.read_csv(file, sep='\t')
        name = name + '|' + merge_column
        new[name]=df[merge_column]

    print('saving {} ({} x {})'.format(output_filename, *new.shape))
    new.to_csv(output_filename,sep='\t',index=False)
    print('done')

if __name__ == '__main__':
    main(sys.argv)
