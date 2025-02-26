#!/usr/bin/env python
import os, sys

#1: input -- taxonomy.table file
#2: output -- all_taxon.taxonomy file

inpath = sys.argv[1]
outpath = sys.argv[2]

# kept columns
columns = ['tax_id', 'parent_id', 'rank', 'tax_name',
    'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

lines = []
with open(inpath, 'r') as f:
    col_to_idx = {}
    for line in f:
        parts = list(eval(line.strip()))
        # first line, map idx to column name
        if parts[0] == 'tax_id':
            for i, col in enumerate(parts):
                if col in columns:
                    col_to_idx[col] = i
        else:
            row = [parts[col_to_idx[col]] for i, col in enumerate(columns)]
            lines.append(row)
with open(outpath, 'w') as f:
    # write header
    f.write(','.join(columns) + '\n')
    for row in lines:
        # if empty, write NA
        f.write(','.join([x if x != '' else 'NA' for x in row]) + '\n')
