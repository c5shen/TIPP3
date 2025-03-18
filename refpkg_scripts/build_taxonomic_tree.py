#!/usr/bin/env python

'''
Reads in a taxonomy file, outputs a taxonomic tree in newick format
Created on Jan 31, 2012

python2.7 build_taxonomic_tree.py ../data/taxonomy/all_taxon.taxonomy taxonomic.tree
python2.7 build_taxonomic_tree.py ../data/taxonomy/all_taxon.taxonomy taxonomic.tree ../data/taxonomy/species.mapping

Converted for use in Python 3 by Mike Nute some time in 2018.

@author: namphuon

Edited by Chengze Shen @ 2025.02.26
'''
from dendropy import Tree, Node, Taxon
import sys
import os


if __name__ == '__main__':
    taxonomyFile = sys.argv[1]
    speciesList = sys.argv[2]
    speciesMap = sys.argv[3]
    taxonomyTree = sys.argv[4]

    species = {}
    lines = open(speciesList,'r')
    for line in lines:
        # skip first line which is the header
        if line.startswith('seqname') or line.startswith('tax_id'):
            continue
        items = line.strip().split(',')
        species[items[-1]] = items[-1]

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
    lines.close()

    # find the node with
    _tmp = tree.find_node_with_taxon_label('131567')
    print("Child taxids at root (131567): "
         f"{','.join(cnode.taxon.label for cnode in _tmp.child_nodes())}")
    #print(f"pruning {len(prune)}/{len(nodes_dict)} nodes")
    for taxa in prune:
        nodes_dict[taxa].label = ''

    # find all unifurcating nodes that have taxids (i.e., having some 
    # sequences that would be assigned to them). These nodes, if not dealt
    # with explicitly, will be contracted ("vanished") after running
    # tree.suppress_unifurcations()
    taxid_reassign = {}
    for node in tree.preorder_node_iter():
        #print(str(node.taxon))
        if (len(node.child_nodes()) == 1 and
                node.taxon and node.taxon.label in species):
            print(f"WARNING: node (taxid {node.taxon.label}) unifurcates "
                  "and have assigned sequences!")
            for cnode in node.preorder_iter():
                if (cnode.taxon.label != node.taxon.label and
                        cnode.taxon.label in species):
                    taxid_reassign[cnode.taxon.label] = node.taxon.label
                    print("\tRe-assigning sequences with taxid "
                          f"{cnode.taxon.label} to {node.taxon.label}")
            # remove the only child so that no more unifurcating below
            # this node
            node.remove_child(node.child_nodes()[0])

    # update speciesList and speciesMap, if len(taxid_reassign) != 0
    num_reassigned, num_total = 0, 0
    if len(taxid_reassign) != 0:
        print("Reassigning taxids to sequences, "
                f"backing up {speciesList} and {speciesMap} ...")
        os.system(f'cp {speciesList} {speciesList}.backup')
        os.system(f'cp {speciesMap} {speciesMap}.backup')

        lines = []
        with open(speciesMap, 'r') as f:
            lines = f.read().strip().split('\n')
        num_total = len(lines) - 1

        sl_fptr = open(speciesList, 'w')
        sm_fptr = open(speciesMap, 'w')
        for line in lines:
            if line.startswith('seqname'):
                sl_fptr.write('tax_id\n')
                sm_fptr.write(line.strip() + '\n')
            else:
                items = line.strip().split(',')
                ori_taxid = items[-1]
                # this sequence needs to have a reassigned taxid
                if ori_taxid in taxid_reassign:
                    print(f"{items[0]} taxid reassigned from "
                          f"{ori_taxid} to {taxid_reassign[ori_taxid]}")
                    items[-1] = taxid_reassign[ori_taxid]
                    num_reassigned += 1
                sl_fptr.write(items[-1] + '\n')
                sm_fptr.write(','.join(items) + '\n')
                num_total += 1
        sl_fptr.close()
        sm_fptr.close()
        print(f"Reassigned taxids for {num_reassigned}/{num_total} sequences") 

    # tree.delete_outdegree_one_nodes()
    tree.suppress_unifurcations()

    output = open(taxonomyTree, 'w')
    output.write(str(tree) + ";\n");

    output.close()
