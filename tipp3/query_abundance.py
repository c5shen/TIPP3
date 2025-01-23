import os, sys, json
from collections import defaultdict

from tipp3.configs import Configs
from tipp3 import get_logger
from tipp3.jobs import TIPPJsonMergerJob

#from tqdm import tqdm
#from tipp3.helpers.general_tools import tqdm_styles

import concurrent.futures

_LOG = get_logger(__name__)

# taxonomic ranks
ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

'''
Function to extract taxonomic information from all query read placements
and aggregate to form an abundance profile
'''
def queryAbundance(refpkg, query_placement_paths, pool, lock):
    # (1) obtain classification
    _LOG.info("Obtaining read classification from each marker gene")
    classification_paths = {} 

    # failsafe when no placement file are given
    if len(query_placement_paths) == 0:
        _LOG.warning("No placements are found for classification or " \
                "abundance profile, returning...")
        return

    futures = []
    for marker, query_placement_path in query_placement_paths.items():
        clas_outdir = os.path.join(Configs.outdir, 'query_classifications',
                marker)
        if not os.path.isdir(clas_outdir):
            os.makedirs(clas_outdir)

        # necessary files
        taxonomy_path = refpkg[marker]['taxonomy']
        # temp fix mapping issue with taxonomy 
        if not os.path.exists(taxonomy_path):
            taxonomy_path = os.path.join(os.path.dirname(taxonomy_path),
                    'all_taxon.taxonomy')
        mapping_path = refpkg[marker]['seq-to-taxid-map']
        classification_path = os.path.join(clas_outdir, 'placement.classification')
        reordered_placement_path = os.path.join(clas_outdir,
                'placement.reordered.jplace')
        futures.append(pool.submit(getClassification, marker, taxonomy_path,
            mapping_path, query_placement_path, reordered_placement_path,
            clas_outdir, classification_path, lock))

    for future in concurrent.futures.as_completed(futures):
        marker, classification_path = future.result()
        _LOG.info(f"Classification completed on {marker}")
        classification_paths[marker] = classification_path

    # (2) filter classification based on given support value
    support_value = '0.95'
    try:
        support_value = str(
                getattr(Configs, Configs.placement_method).support_value)
    except (AttributeError, ValueError) as e:
        pass
    _LOG.info(f"Filtering with support value={support_value}")

    filtered_paths = {}
    for marker, classification_path in classification_paths.items():
        clas_outdir = os.path.join(Configs.outdir, 'query_classifications',
                marker)
        taxonomy_path = refpkg[marker]['taxonomy']
        if not os.path.exists(taxonomy_path):
            taxonomy_path = os.path.join(os.path.dirname(taxonomy_path),
                    'all_taxon.taxonomy')
        filtered_path = os.path.join(clas_outdir,
                f"placement.classification.{support_value.split('.')[-1]}")
        filterClassification(taxonomy_path,
                classification_path, filtered_path, float(support_value))
        filtered_paths[marker] = filtered_path

    # Updated @ 1.22.2025 - Chengze Shen
    #   - aggregate all classifications and output to a file as well
    # (2.5) aggregate all classifications
    all_classification_path = os.path.join(Configs.outdir,
            'query_classifications.tsv')
    header = False
    with open(all_classification_path, 'w') as f:
        for marker, filtered_path in filtered_paths.items():
            with open(filtered_path, 'r') as fptr:
                lines = fptr.read().strip().split('\n')
            if not header:
                f.write(lines[0] + '\n')
                header = True
            f.write('\n'.join(lines[1:]) + '\n')

    # (3) aggregate classifications to form an abundance profile
    abundance_profile = getAbundanceProfile(refpkg, filtered_paths)

'''
Function to parse filtered classification files and aggregate their abundances
to produce the abundance profile on all taxonomic levels
'''
def getAbundanceProfile(refpkg, filtered_paths):
    global ranks
    _LOG.info("Aggregating abundances for a profile")

    # initializing abundance profile at each taxonomic level
    abundance_profile = {}
    for rank in ranks:
        abundance_profile[rank] = defaultdict(float)
    
    # read in taxid_map from refpkg
    taxid_map = parseTaxonomy(refpkg['taxonomy']['taxonomy'])

    # aggregating results from all markers
    for marker, filtered_path in filtered_paths.items():
        updateAbundanceProfile(filtered_path, abundance_profile)

    # compute percentage at each rank
    _LOG.info("Writing abundance profiles at each taxonomic level to "
            f"{Configs.outdir}") 
    for rank in ranks:
        rank_sum = sum(abundance_profile[rank].values())
        for taxid in abundance_profile[rank].keys():
            abundance_profile[rank][taxid] /= rank_sum

        # write to Configs.outdir as abundance.<rank>.tsv
        outpath = os.path.join(Configs.outdir, f"abundance.{rank}.tsv")
        slist = sorted([(k, v) for k, v in abundance_profile[rank].items()],
                key=lambda x: x[1], reverse=True)
        with open(outpath, 'w') as f:
            f.write('taxa\ttaxid\tabundance\n')
            for (taxid, abundance) in slist:
                # unclassified
                if taxid == 0:
                    f.write('unclassified\t0\t{}\n'.format(abundance))
                else:
                    taxname = taxid_map[taxid][0]
                    f.write('{}\t{}\t{}\n'.format(taxname, taxid, abundance))
        _LOG.info(f"Finished written rank: {rank}") 

'''
Function to update a given abundance profile with the new results
read from the given inpath.
NOTE: updates are "counts", proportions will be calculated later
'''
def updateAbundanceProfile(inpath, abundance_profile):
    global ranks
    with open(inpath, 'r') as f:
        lines = f.read().strip().split('\n')
        # skip header
        for i in range(1, len(lines)):
            parts = lines[i].split('\t')
            taxon = parts[0]
            # fields 1-7 should be what we want
            fields = [int(x) if x != 'NA' else 0 for x in parts[1:]]
            # not classified at all, ignore
            if sum(fields) == 0:
                continue
            
            for j in range(len(ranks)):
                abundance_profile[ranks[j]][fields[j]] += 1

'''
Function to parse a given taxonomy file (all_taxon.taxonomy)
'''
def parseTaxonomy(taxonomy_path):
    taxid_map = {}
    with open(taxonomy_path, 'r') as f:
        for line in f:
            # skip header
            if line.startswith('tax_id'):
                continue
            parts = line.strip().split(',')
            taxid, rank, parent_id, taxname = \
                int(parts[0]), parts[2], int(parts[1]), parts[3]
            taxid_map[taxid] = (taxname, parent_id, rank)
    return taxid_map

'''
Function to reorder a .jplace file to a standard format for reading
'''
def reorderJson(inpath, outpath):
    fh = open(inpath, 'r')
    s_tree = ""
    obj = json.load(fh)
    newobj = {}

    # obtain old order
    old_order = obj['fields']

    new_plc_order = ['edge_num', 'likelihood', 'like_weight_ratio',
            'distal_length', 'pendant_length']
    
    # (1) tree
    newobj['tree'] = obj['tree']
    s_tree = obj['tree'].strip()

    # (2) placements in different order
    newobj['placements'] = [] 
    for placement in obj['placements']:
        if 'nm' in placement:
            ori_p, ori_n = placement['p'], placement['nm']
            ori_n = [x[0] for x in ori_n]
        else:
            ori_p, ori_n = placement['p'], placement['n']
        new_p = []
        for pp in ori_p:
            #new_p.append([pp[1], pp[3], pp[2], pp[0], pp[-1]])
            pp_map = {old_order[i]: pp[i] for i in range(len(pp))}
            tmp = [pp_map[new_plc_order[i]] for i in range(len(new_plc_order))]
            new_p.append(tmp)
        _dict = {'p': new_p, 'n': ori_n}
        newobj['placements'].append(_dict)
    # (3) remaining fields 
    newobj['metadata'] = obj['metadata']
    newobj['version'] = obj['version']
    newobj['fields'] = new_plc_order

    # write to outpath
    ofh = open(outpath, 'w')
    json.dump(newobj, ofh, indent = 4)
    ofh.close()
    fh.close()

    # create stdindata
    s_tree = s_tree.replace("{", "[")
    s_tree = s_tree.replace("}", "]")
    mergeinput = [s_tree]
    mergeinput.append('{}\n{}'.format(s_tree, outpath))
    mergeinput.append("")
    mergeinput.append("")
    s_mergeinput = '\n'.join(mergeinput)
    return s_mergeinput

'''
Function to obtain the taxonomic classification based on placement result,
for the given marker gene taxonomic tree
'''
def getClassification(marker, taxonomy_path, mapping_path,
        ori_jplace_path, jplace_path,
        outdir, classification_path, lock):
    # reorder output .jplace file from the placement method
    stdindata = reorderJson(ori_jplace_path, jplace_path)
    
    # run TIPPJsonMerger software; skip if classification is already there
    if not (os.path.exists(classification_path) \
            and os.stat(classification_path).st_size > 0):
        tippjsonmerger_job = TIPPJsonMergerJob(path=Configs.tippjsonmerger_path,
                taxonomy_path=taxonomy_path, mapping_path=mapping_path,
                outdir=outdir, classification_path=classification_path)
        _ = tippjsonmerger_job.run(stdin=stdindata)

    return marker, classification_path

'''
Function to load in taxonomy
'''
def load_taxonomy(taxonomy_file, lower=True):
    global ranks
    f = open(taxonomy_file, 'r')

    # First line is the keywords for the taxonomy, need to map the keyword to
    # the positional index of each keyword
    results = f.readline().lower().replace('"', '').strip().split(',')
    key_map = dict([(results[i], i) for i in range(0, len(results))])

    # Now fill up taxonomy, level maps keep track of what taxa exist at each
    # level, taxon_map keep track of entire taxonomy
    taxon_map = {}
    level_map = {"species": {}, "genus": {}, "family": {}, "order": {},
            "class": {}, "phylum": {}, "superkingdom": {}}

    for line in f:
        results = line.replace('"', '').strip()
        if (lower):
            results.lower()
        results = results.split(',')
        # insert into taxon map
        taxon_map[results[0]] = results

        # insert into level map
        for level in ranks:
            if (results[key_map[level]] == ''):
                continue
            else:
                if (results[key_map[level]] not in level_map[level]):
                    level_map[level][results[key_map[level]]] = {}
                level_map[level][results[key_map[level]]][results[0]] = \
                    results[0]
    return (taxon_map, level_map, key_map)

'''
Function to filter classification output from a marker gene by a given
support value
'''
def filterClassification(taxonomy_path, classification_path, filtered_path,
        threshold):
    global ranks
    taxon_map, level_map, key_map = load_taxonomy(taxonomy_path)
    
    class_in = open(classification_path, 'r')
    level_map_hierarchy = {"species": 0, "genus": 1, "family": 2, "order": 3,
                "class": 4, "phylum": 5, "superkingdom": 6,
                "root": 7}
    # Need to keep track of last line so we can determine when we switch to
    # new classification
    old_name, old_id, old_rank, old_probability = "", "", "", 1

    # keep track of all fragment names
    names = {}
    classification = {}
    for line in class_in:
        results = line.strip().split(',')
        if (len(results) > 5):
            results = [results[0], results[1], results[2],
                       results[-2], results[-1]]
        (name, id, taxname, rank, probability) = (
            results[0], results[1], results[2], results[3], float(results[4]))
        names[name] = name
        if (name != old_name):
            # when we switch to new fragment, output last classification for
            # old fragment
            if (old_name != ""):
                lineage = taxon_map[old_id]
                output_line = [old_name]
                for level in ranks:
                    clade = lineage[key_map[level]]
                    if (clade == ""):
                        clade = "NA"
                    output_line.append(clade)
                classification[old_name] = output_line
            old_name = name
            old_rank = "root"
            old_probability = 1
            old_id = '1'

        # Switch to new rank if the new probability is higher than threshold
        # and our rank is more specific than our original rank
        if (
            rank in level_map_hierarchy and
            (level_map_hierarchy[old_rank] > level_map_hierarchy[rank]) and
            (probability > threshold)
           ):
            old_rank = rank
            old_probability = probability
            old_id = id
        # Switch to new rank if the new rank matches old rank but has higher
        # probability
        elif (
              rank in level_map_hierarchy and
              (level_map_hierarchy[old_rank] == level_map_hierarchy[rank]) and
              (probability > old_probability)
             ):
            old_rank = rank
            old_probability = probability
            old_id = id

    if old_id in taxon_map:
        lineage = taxon_map[old_id]
        output_line = [old_name]
        for level in ranks:
            clade = lineage[key_map[level]]
            if (clade == ""):
                clade = "NA"
            output_line.append(clade)
        classification[name] = output_line
    class_in.close()

    # write to filtered_path    
    with open(filtered_path, 'w') as class_out:
        class_out.write("fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n")
        keys = list(classification.keys())
        keys.sort()
        for frag in keys:
            class_out.write("%s\n" % "\t".join(classification[frag]))
