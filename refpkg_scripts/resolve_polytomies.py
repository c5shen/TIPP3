#!/usr/bin/env python

# randomly resolve polytomies to make sure input tree to raxml-ng is binary
import os, sys
import dendropy

#1: input tree
#2: output tree
inpath = sys.argv[1]
outpath = sys.argv[2]

intree = dendropy.Tree.get_from_path(inpath, schema='newick')
# remove internal node labels
for node in intree.internal_nodes():
    node.label = None

# resolve polytomies
intree.resolve_polytomies()

# write to outpath
intree.write_to_path(outpath, schema='newick')