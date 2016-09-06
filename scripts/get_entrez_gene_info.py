#!/usr/bin/python

import sys
import os
import csv
import datetime
import time
from collections import defaultdict

def time_stamp():
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)

def log(message):
    print('{}|{}'.format(time_stamp(), message), file=sys.stderr)

#?? How to validate the command line argument; error handling
gene_info = sys.argv[1]         #'NCBI_annotation_2016_06_09/NCBI_gene_info_2016_07_12'
gene_expr = sys.argv[2]         #'gene_diffexp.txt' #C.Plus_v_Six.Plus_gene.foldchange.1.5.txt'
tax_id = sys.argv[3]            #'9606'
outfile_tag = os.path.basename(gene_expr.replace('.txt', '_annot.txt')) # instead of writing to a file write to stdout
outfile_name = os.path.abspath(outfile_tag)

log('reading gene info')
gene_details = defaultdict(list)
with open(gene_info,'r') as f:
    next(f) # skip header
    reader=csv.reader(f,delimiter='\t')
    for row in reader:
        geneinfo_tax_id = row[0]
        geneinfo_gene_id =row[1]
        geneinfo_gene_symbol = row[2]
        if row[0] == tax_id:
            gene_details[row[2]].append(row[1]) 
            gene_details[row[2]].append(row[8])

header_line=[]
with open(gene_expr, 'r') as x, open(outfile_name, 'w') as fh: 
    reader=csv.reader(x,delimiter='\t')

    for row in reader:
        row = [col.strip() for col in row]
        if row[0] == '#test_id':
            row[1] = 'gene_symbol'
            row[2] = 'gene_id'
            row.insert(3, 'gene_desc')
            print('\t'.join(row), file=fh)
        elif row[0].startswith('#'):
            print('\t'.join(row), file=fh)
        else:
            left = row[0:2]
            right = row[3:]
            gene_name = row[1]
            gene_id_symbol = gene_details.get(gene_name, ['.','.']) # If key is not available then returns default value 

            outline = "\t".join(left + gene_id_symbol + right)  
            print(outline, file=fh)
    
#fh.close()
log('done')






# for line in open(gene_expr):
#     li=line.strip()
#     if not li.startswith("#"):
#         print(li)
        
# rdr=csv.reader(fp,delimiter='\t')
# header = csv.DictReader(filter(lambda row: row[0]=='#', fp))
# rdr = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter ='\t')
# print(rdr.fieldnames[3])
# list =[]
# for rows in rdr:
#    print(rows)
# fp.close()


# with open(gene_expr,'r') as f:
#     reader=csv.reader(f,delimiter='\t')
#     for row in reader:
#         print(row[1])
#         name = row[0].rstrip() #gene symbols in this file have a trailing space
#         if name in gene_details:
# #            print(name, "\t", "\t".join(gene_details[name]))
#  