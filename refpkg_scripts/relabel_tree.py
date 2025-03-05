#!/usr/bin/env python
import dendropy
import sys, os

def sanity_check_relabelled_tree(ref_lin_path, relabelled_lin_path,
        outdir, treename):
    ref_map = {}; ref_line_num = 0
    with open(ref_lin_path, 'r') as f:
        for line in f:
            ref_line_num += 1
            ref_map[line.strip('\n').strip()] = 1
    ref_set = set(ref_map.keys())
    assert len(ref_set) == ref_line_num

    relabelled_map = {}; relabelled_line_num = 0
    with open(relabelled_lin_path, 'r') as f:
        for line in f:
            relabelled_line_num += 1
            relabelled_map[line.strip('\n').strip()] = 1
    relabelled_set = set(relabelled_map.keys())
    assert len(relabelled_set) == relabelled_line_num
    
    with open(outdir + '/{}-lineage.sanity_check.out'.format(treename), 'w') as f:
        thesame = len(ref_set) == len(relabelled_set)
        thesame = thesame and (len(ref_set.difference(relabelled_set)) == 0)
        print('sanity check correctness', thesame)
        f.write('{}\n'.format(thesame))

def write_lineage(tree, outpath):
    with open(outpath, 'w') as f:
        for node in tree.leaf_nodes():
            path = [node.taxon.label]
            for pp in node.ancestor_iter():
                if pp.label:
                    path.append(pp.label)
            f.write('{}\n'.format(','.join(path)))
    # sort lineage for easier comparison
    os.system(f"sort {outpath} -o {outpath}")

# max system recursion limit hard encoding to a large number
# a temp fix for dendropy tree recursion issues
sys.setrecursionlimit(20000)

# Arguments: 1: a "reference" tree with internal taxon labels
#               example: 99_otus.tree
#            2: a "new" tree that is a refinement of the refrence tree, but
#               lacks internal nodes
#               example: reference-gg-raxml-bl-rooted.tre
#            3: the output: the new tree with internal nodes added
#               example: reference-gg-raxml-bl-rooted-relabelled.tre
# Assumes there is a single node called k__Archaea to root the output tree on
#
# Added by Chengze Shen
#            4: the root label: should be inputted by the user
#               (not necesarily k__Archaea)

taxa = dendropy.TaxonNamespace()

t = dendropy.Tree.get_from_path(
    src=sys.argv[1], schema='newick', taxon_namespace=taxa)
t2 = dendropy.Tree.get_from_path(
    src=sys.argv[2], schema='newick', taxon_namespace=taxa)
treename = sys.argv[2].split('/')[-1]
outdir = '/'.join(os.path.abspath(sys.argv[3]).split('/')[:-1])
# root_path=sys.argv[4]

# # read in the tree label - the first line in <root_path>
# root_label = ''
# with open(root_path, 'r') as f:
#     line = f.readline().strip('\n').strip()
#     root_label = line
root_label = sys.argv[4]

a = t.find_node_with_label(label=root_label)
t.encode_bipartitions(suppress_unifurcations=False,
        collapse_unrooted_basal_bifurcation=False)
#print(t.seed_node)
#print(t.seed_node.child_nodes())
#print(len(t.internal_nodes()), len(t.leaf_nodes()))
#print(len(t2.internal_nodes()), len(t2.leaf_nodes()))

# ONE SPECIAL scenario: the seed node is NOT labeled 131567, meaning that 
# there is only one superkingdom group in the tree (e.g., only Bacteria {2}
# and no Archaea {2157}).
# To correctly seed the refined tree in this case, we need to find just one
# bipartition from the seed_node.child_nodes() (e.g., the first one), and seed
# at that edge by creating a new node using the original edge length
seed_node_children = {}
if t.seed_node.label != '131567':
    #for cnode in t.seed_node.child_nodes():
    #    print(cnode.label)
    #print("")
    #for cnode in t2.seed_node.child_nodes():
    #    print(len(cnode.leaf_nodes()))
    #print(t.seed_node.label)
    root_label = t.seed_node.label
    samp_node = t.seed_node.child_nodes()[0]
    #t_samp_num_leaf = len(samp_node.bipartition.leafset_taxa(
    #    taxon_namespace=t.taxon_namespace))
    #print(samp_node.label)
    #exit()

    print('reseeding {} (only one superkingdom group)'.format(sys.argv[2]))
    t2.seed_node.label = root_label
    #newroot = dendropy.Node(label=root_label)
    #cnode = samp_node
    #bip = cnode.edge.bipartition; bip.is_mutable = False

    #if bip in t2.bipartition_edge_map:
    #    t2edge = t2.bipartition_edge_map[bip]
    #    t2edge_len = t2edge.length
    #    t2h, t2t = t2edge.head_node, t2edge.tail_node
    #    t2t.remove_child(t2h)
    #    newroot.add_child(t2h); newroot.parent_node = t2t

    #    t2.reseed_at(newroot, update_bipartitions=True,
    #            suppress_unifurcations=False,
    #            collapse_unrooted_basal_bifurcation=False)
    #else:
    #    raise ValueError('bip <{}> in unrefined tree not in refined tree'.format(
    #        cnode))
else:
    for cnode in t.seed_node.child_nodes():
        # map number of leaves to the corresponding cnode label
        # mostly, mostly, will not be an issue since the number of
        # leaves for bacteria, archaea, eukaryota should be quite different
        # under normal conditions
        _num_leaf = len(cnode.bipartition.leafset_taxa(
                taxon_namespace=t.taxon_namespace))
        seed_node_children[_num_leaf] = cnode.label

    # find the bipartition in t2 that's the same as the bac or arch bipartition in t
    # add a new node that splits the original edge in t2. The new node will be the
    # new seed node 131567
    print('reseeding {}'.format(sys.argv[2]))
    newroot = dendropy.Node(label=root_label)
    for cnode in t.seed_node.child_nodes():
        bip = cnode.edge.bipartition; bip.is_mutable = False
        leafset_size = len(bip.leafset_taxa(taxon_namespace=t.taxon_namespace))
        if bip in t2.bipartition_edge_map:
            t2edge = t2.bipartition_edge_map[bip]
            t2h, t2t = t2edge.head_node, t2edge.tail_node

            # - remove head_node from tail_node child
            # - add newroot as the intermediate node between the original head_node
            #   and tail_node
            t2t.remove_child(t2h)
            newroot.add_child(t2h); newroot.parent_node = t2t
            
            # then, reseed ('soft reroot') at the newroot
            t2.reseed_at(newroot, update_bipartitions=True,
                    suppress_unifurcations=False,
                    collapse_unrooted_basal_bifurcation=False)
            #t2.encode_bipartitions(suppress_unifurcations=False,
            #        collapse_unrooted_basal_bifurcation=False,
            #        is_bipartitions_mutable=False)
            #t2t.bipartition._split_bitmask = ~t2h.bipartition._split_bitmask
            #t2t.bipartition.compile_bipartition()
            break
        else:
            raise ValueError('bip <{}> in unrefined tree not in refined tree'.format(
                cnode))
for edge in t2.seed_node.child_edges():
    #print(edge.length)
    if not edge.length:
        edge.length = 0

missing = []
mapped = 0
#mapped_node = {}
#overwritten = False
print('relabeling {}'.format(sys.argv[2]))
for n in t.preorder_node_iter():
    if n.is_internal() and n.label is not None:
        n.edge.bipartition.is_mutable = False
        if n.edge.bipartition not in t2.bipartition_edge_map:
            missing.append(n.edge.bipartition)
        else:
            t2n = t2.bipartition_edge_map[n.edge.bipartition].head_node
            #rpr = hex(id(t2n))
            ## overwriting a node that has already mapped a bipartition
            #if rpr in mapped_node:
            #    overwritten = True
            #mapped_node[rpr] = n.edge.bipartition
            t2n.label = n.label
            mapped = mapped + 1
print('\tmapped', mapped)
print('\tmissing', len(missing))

# somehow child nodes at the root are not labelled correctly
# due to overwriting
if t.seed_node.label != '2':
    for cnode in t2.seed_node.child_nodes():
        _num_leaf = len(cnode.leaf_nodes())
        if _num_leaf in seed_node_children:
            cnode.label = seed_node_children[_num_leaf]
        #if len(cnode.leaf_nodes()) == t_bac_num_leaf:
        #    cnode.label = bac
        #elif len(cnode.leaf_nodes()) == t_arch_num_leaf:
        #    cnode.label = arch
        else:
            raise ValueError("there must be something wrong with the bipartitions in the refined tree.")
else:
    # in this case, the bipartition won't be labelled incorrectly since the
    # other label is always going to be None (and not considered)
    pass

t2.write_to_path(dest=sys.argv[3], schema='newick', suppress_rooting=False,
                 suppress_internal_node_labels=False)

# write lineage of the ref (unrefined) tree and relabelled tree to local
# then do a sanity check to make sure they are the same
ref_lin_path = outdir + '/{}-lineage.unrefined.out'.format(treename)
relabelled_lin_path = outdir + '/{}-lineage.relabelled.out'.format(treename)
write_lineage(t, ref_lin_path)
write_lineage(t2, relabelled_lin_path)
sanity_check_relabelled_tree(ref_lin_path, relabelled_lin_path, outdir, treename)
