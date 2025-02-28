#!/usr/bin/env python

'''
Reads in a taxonomy file, outputs a taxonomic tree in newick format
Created on Jan 31, 2012

python2.7 build_taxonomic_tree.py ../data/taxonomy/all_taxon.taxonomy taxonomic.tree
python2.7 build_taxonomic_tree.py ../data/taxonomy/all_taxon.taxonomy taxonomic.tree ../data/taxonomy/species.mapping

Converted for use in Python 3 by Mike Nute some time in 2018.

@author: namphuon

Edited slightly by Chengze Shen @ 2025.02.26
'''
from dendropy import Tree, Node, Taxon
import sys
import os


if __name__ == '__main__':
  taxonomyFile = sys.argv[1]
  speciesList = sys.argv[2]
  taxonomyTree = sys.argv[3]

  species = {}
  lines = open(speciesList,'r')
  for line in lines:
    species[line.strip()] = line.strip()
    
  lines = open(taxonomyFile,'r')
  header = lines.readline()
  nodes_dict = {}

  #Read first line, root node
  line = lines.readline()
  results = line.strip().split(',')
  tree = Tree()  
  root = Node()
  root.__dict__['label'] = results[0].replace("\"","")
  nodes_dict[results[0].replace("\"","")] = root
  
  prune = ['1']
  
  #Add root node to tree
  tree.__dict__['_seed_node'].add_child(root)
  for line in lines:
    results = line.strip().split(',')
    node = Node();
    node.__dict__['label'] = results[0].replace("\"","")
    node.taxon = Taxon(results[0].replace("\"",""))
    nodes_dict[results[0].replace("\"","")] = node
    nodes_dict[results[1].replace("\"","")].add_child(node)
    if results[0].replace("\"","") not in species:
      prune.append(results[0].replace("\"",""))
      
  for taxa in prune:
    nodes_dict[taxa].label=''
  # tree.delete_outdegree_one_nodes()
  tree.suppress_unifurcations()
      
  output = open(taxonomyTree, 'w')
  output.write(str(tree) + ";");
  
  output.close()
  lines.close()

#Usage python build_taxonomic_tree.py taxonomy.table species.updated.txt unrefined.taxonomy
