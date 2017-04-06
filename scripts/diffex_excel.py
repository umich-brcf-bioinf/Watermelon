#!/bin/env python
from __future__ import print_function, absolute_import, division

import argparse
import csv
import datetime
import itertools
import math
import os
import sys
import time

from xlsxwriter.workbook import Workbook

DESCRIPTION = \
'''Combines annotated gene and isoform tab-separated text files into an Excel workbook.
The first line of the file is assumed to be the header. If glossary and info file
paths are specified, the contents of those files will  be copied to new worksheets.'''

#validate required fields present
#validate right number of rows/columns
#validate new field names
#validate correct sheets
#$ scripts/diffex_excel genes.txt isoforms.txt comparison.xlsx 
# files don't exist, output dir doesn't exist
# better error processing for missing newlines
# NaN values

CELL_COMMENT_FORMAT = {'width': 400}
DELIMITER = '\t'
REQUIRED_FIELDS = argparse.Namespace(gene_symbol='gene_symbol', gene_id='gene_id')
NULL = '.'
HEADER_PREFIX = '#'


class _Formatter(object):
    def __init__(self, workbook):
        self.default_format = workbook.add_format()
        self.header_format = workbook.add_format({'bold': True})
        self.url_format = workbook.add_format({'font_color': 'blue', 'underline':  1})
        self.num_format = workbook.add_format({'num_format':'0'})
        self.float_format = workbook.add_format({'num_format':'0.000'})
        self.scientific_format = workbook.add_format({'num_format':'0.00E+00'})
        self.field_formats = {'p_value' : self.scientific_format,
                              'q_value': self.scientific_format}

    @staticmethod
    def _is_number(value):
        if not value:
            return False
        try:
            value = float(value)
            return not (math.isinf(value) or math.isnan(value)) 
        except ValueError:
            return False

    def _guessed_format(self, value):
        if self._is_number(value):
            value = float(value)
            if value == int(value):
                return self.num_format
            else:
                return self.float_format
        else:
            return self.default_format

    def add_field_format(self, field_name, format):
        self.field_formats[field_name] = format

    def get_format(self, field_name, value):
        return self.field_formats.get(field_name, self._guessed_format(value))


class _NcbiGeneHyperlink(object):
    _NCBI_URL_FMT = 'https://www.ncbi.nlm.nih.gov/gene/?term={}' 
    
    def __init__(self, formatter, required_fields):
        self._formatter = formatter
        self._field_name = 'ncbi_gene_link'
        self._gene_id = required_fields.gene_id
        self._gene_symbol = required_fields.gene_symbol
        
    def _get_gene_hyperlink(self, row):
        result = row.get(self._gene_id)
        if result == NULL:
            result = row.get(self._gene_symbol)
        if result == NULL:
            return NULL 
        return self._NCBI_URL_FMT.format(result)

    def _write_value(self, worksheet, formatter, row_index, col_index, row, field_name):
        value = row[self._gene_symbol]
        hyperlink = self._get_gene_hyperlink(row)
        format = formatter.get_format(field_name, value)
        if hyperlink != NULL: 
            worksheet.write_url(row_index, col_index,
                                hyperlink,
                                format,
                                value)
        else:
            worksheet.write(row_index, col_index, value)

    def inject_field(self, header, data_iter, field_writers, delimiter=DELIMITER):
        self._formatter.add_field_format(self._field_name, self._formatter.url_format)
        field_writers = dict(field_writers)
        field_writers[self._field_name] = self._write_value
        field_names = header.split(delimiter)
        new_field_index = field_names.index(self._gene_symbol)
        if self._field_name not in field_names:
            field_names.insert(new_field_index, self._field_name)
        def new_data_iter(data_iter, new_field_index, delimiter):
            for line in data_iter:
                fields = line.split(delimiter)
                fields.insert(new_field_index, '')
                yield delimiter.join(fields)
        return field_names, new_data_iter(data_iter, new_field_index, delimiter), field_writers


def _log(message):
    def _time_stamp():
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        return(st)
    print('{}|diffex_excel|{}'.format(_time_stamp(), message), file=sys.stderr)

def _parse_header_data(reader):
    header = reader.readline()
    return header, reader

def _write_header_row(worksheet, field_names, format):
    for col_index, value in enumerate(field_names):
        worksheet.write(0, col_index, value, format)

def _write_value(worksheet, formatter, row_index, col_index, row, field_name):
    value = row[field_name]
    format = formatter.get_format(field_name, value)
    worksheet.write(row_index, col_index, value, format)

def _write_data_row(worksheet, field_names, formatter, field_writers, row_index, row):
    for col_index, field_name in enumerate(field_names):
        field_writer = field_writers.get(field_name, _write_value)
        field_writer(worksheet, formatter, row_index, col_index, row, field_name)

def _add_worksheet(workbook, formatter, worksheet_name, input_filepath):
    worksheet = workbook.add_worksheet(worksheet_name)
    field_writers = {}
    ncbi_gene_hyperlink = _NcbiGeneHyperlink(formatter, REQUIRED_FIELDS)
    _log('reading {}'.format(input_filepath))
    with open(input_filepath, 'rU') as input_file:
        header, data = _parse_header_data(input_file)
        _validate_required_fields(header.split(DELIMITER), input_filepath)
        field_names, data, field_writers = ncbi_gene_hyperlink.inject_field(header, data, field_writers, DELIMITER)
        reader = csv.DictReader(data, fieldnames=field_names, delimiter=DELIMITER)
        _write_header_row(worksheet, reader.fieldnames, formatter.header_format)
        for row_index, row in enumerate(reader):
            _write_data_row(worksheet, reader.fieldnames, formatter, field_writers, row_index + 1, row)
    worksheet.freeze_panes(1, 0)

def _add_glossary(workbook, formatter, worksheet_name, input_filepath):
    worksheet = workbook.add_worksheet(worksheet_name)
    field_writers = {}
    _log('adding glossary from {}'.format(input_filepath))
    with open(input_filepath, 'rU') as input_file:
        header, data = _parse_header_data(input_file)
        field_names = header.split(DELIMITER)
        reader = csv.DictReader(data, fieldnames=field_names, delimiter=DELIMITER)
        _write_header_row(worksheet, reader.fieldnames, formatter.header_format)
        for row_index, row in enumerate(reader):
            _write_data_row(worksheet, reader.fieldnames, formatter, field_writers, row_index + 1, row)
    worksheet.freeze_panes(1, 0)

def _add_info(workbook, formatter, worksheet_name, input_filepath):
    worksheet = workbook.add_worksheet(worksheet_name)
    field_writers = {}
    _log('adding info from {}'.format(input_filepath))
    with open(input_filepath, 'rU') as input_file:
        for line_num, line in enumerate(input_file):
            worksheet.write(line_num, 0, line)

def _validate_required_fields(field_names, input_filepath):
    required_field_set = set(vars(REQUIRED_FIELDS).values())
    actual_fields = set(field_names)
    missing_fields = sorted(required_field_set - actual_fields)
    if missing_fields:
        msg_fmt = 'Input file [{}] is missing required field(s) [{}].'
        msg = msg_fmt.format(input_filepath, ','.join(missing_fields))
        raise ValueError(msg)

def _parse_command_line_args(sys_argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '-g',
        '--gene_filepath',
        type=str,
        required=True,
        help='path to annotated gene tab-separated text file')
    parser.add_argument(
        '-i',
        '--isoform_filepath',
        type=str,
        required=False,
        help='path to annotated isoform tab-separated text file')
    parser.add_argument(
        '--glossary_filepath',
        type=str,
        required=False,
        help='contents of this tab-separated text file will become glossary sheet')
    parser.add_argument(
        '--info_filepath',
        type=str,
        required=False,
        help='contents of this text file will become info sheet')
    parser.add_argument(
        'output_filepath',
        type=str,
        help='path to Excel file')

    args = parser.parse_args(sys_argv)
    return args 

def main(sys_argv):
    args = _parse_command_line_args(sys_argv)
    workbook = Workbook(args.output_filepath, {'strings_to_numbers': True})
    formatter = _Formatter(workbook)
    _add_worksheet(workbook, formatter, 'genes', args.gene_filepath)
    if args.isoform_filepath:
        _add_worksheet(workbook, formatter, 'isoforms', args.isoform_filepath)
    if args.glossary_filepath:
        _add_glossary(workbook, formatter, 'glossary', args.glossary_filepath)
    if args.info_filepath:
        _add_info(workbook, formatter, 'info', args.info_filepath)

    _log("saving to {}".format(args.output_filepath))
    workbook.close()
    _log("done")

if __name__ == '__main__':
    main(sys.argv[1:])
