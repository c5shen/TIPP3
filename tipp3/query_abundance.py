import csv
import os, json
from collections import defaultdict
import concurrent.futures

from tipp3.configs import Configs
from tipp3 import get_logger

_LOG = get_logger(__name__)

RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

def getSpeciesDetection(detection_thresholds, refpkg, classification_paths,
                        rank='species'):
    """Detect species (or other taxa) based on read classification confidences."""
    _LOG.info(f"Detecting at taxonomic level: {rank.upper()}")
    species_to_marker = parseSpeciesToMarker(
        refpkg['taxonomy']['species-to-marker-map'])
    taxid_map = parseTaxonomy(refpkg['taxonomy']['taxonomy'])

    try:
        rank_idx = RANKS.index(rank)
    except ValueError:
        rank_idx = 0

    # Aggregate likelihood-weight ratios per marker, then average across markers
    detected = defaultdict(float)
    for marker, classification_path in classification_paths.items():
        marker_detected = defaultdict(float)
        marker_cnt = defaultdict(int)

        with open(classification_path, 'r', newline='') as f:
            reader = csv.reader(f)
            for parts in reader:
                if len(parts) < 5:
                    continue
                if parts[-2] != rank:
                    continue
                try:
                    name, taxid, supp = parts[0], int(parts[1]), float(parts[-1])
                except (ValueError, IndexError):
                    continue
                marker_detected[taxid] += supp
                marker_cnt[taxid] += 1

        for taxid in marker_detected:
            detected[taxid] += marker_detected[taxid] / marker_cnt[taxid]

    detected_taxids = sorted(
        detected.keys(),
        key=lambda x: detected[x] / max(len(species_to_marker.get(x, [1])), 1),
        reverse=True)

    outpath = os.path.join(Configs.outdir, f"detected_{rank}_UNFILTERED.tsv")
    _LOG.info(f"Writing UNFILTERED {rank} detection to {outpath}")
    with open(outpath, 'w') as f:
        f.write("taxa\ttaxid\tmarker_confidence\n")
        for taxid in detected_taxids:
            taxname = taxid_map.get(taxid, (str(taxid),))[0]
            n_markers = max(len(species_to_marker.get(taxid, [1])), 1)
            marker_confidence = detected[taxid] / n_markers
            f.write(f"{taxname}\t{taxid}\t{marker_confidence}\n")

    for suffix, thres in detection_thresholds.items():
        outpath = os.path.join(Configs.outdir, f"detected_{rank}_{suffix}.tsv")
        _LOG.info(f"Writing {suffix} (B={thres}) {rank} detection to {outpath}")
        with open(outpath, 'w') as f:
            f.write("taxa\ttaxid\tmarker_confidence\n")
            for taxid in detected_taxids:
                taxname = taxid_map.get(taxid, (str(taxid),))[0]
                n_markers = max(len(species_to_marker.get(taxid, [1])), 1)
                mc = detected[taxid] / n_markers
                if mc < thres:
                    break
                f.write(f"{taxname}\t{taxid}\t{mc}\n")

def getAbundanceProfile(refpkg, filtered_paths):
    """Aggregate filtered classifications into an abundance profile per rank."""
    _LOG.info("Aggregating abundances for a profile")

    abundance_profile = {rank: defaultdict(float) for rank in RANKS}
    taxid_map = parseTaxonomy(refpkg['taxonomy']['taxonomy'])

    for marker, filtered_path in filtered_paths.items():
        updateAbundanceProfile(filtered_path, abundance_profile)

    _LOG.info(f"Writing abundance profiles at each taxonomic level to {Configs.outdir}")
    for rank in RANKS:
        rank_sum = sum(abundance_profile[rank].values())
        if rank_sum == 0:
            _LOG.warning(f"No classifications found at rank '{rank}', skipping.")
            continue
        for taxid in abundance_profile[rank]:
            abundance_profile[rank][taxid] /= rank_sum

        outpath = os.path.join(Configs.outdir, f"abundance.{rank}.tsv")
        slist = sorted(abundance_profile[rank].items(),
                       key=lambda x: x[1], reverse=True)
        with open(outpath, 'w') as f:
            f.write('taxa\ttaxid\tabundance\n')
            for taxid, abundance in slist:
                if taxid == 0:
                    f.write(f'unclassified\t0\t{abundance}\n')
                else:
                    taxname = taxid_map.get(taxid, (str(taxid),))[0]
                    f.write(f'{taxname}\t{taxid}\t{abundance}\n')
        _LOG.info(f"Finished writing rank: {rank}")

def getAllClassification(refpkg, query_placement_paths, pool, lock):
    """Obtain all read classifications from placement results."""
    _LOG.info("Obtaining read classification from all marker genes")
    classification_paths = {}

    if not query_placement_paths:
        _LOG.warning("No placements found for classification. "
                     "This may indicate no reads mapped to marker genes.")
        return {}, {}

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

    filtered_paths = {}
    if Configs.command == 'abundance':
        _LOG.info(f"Filtering with support value={support_value}")
        for marker, classification_path in classification_paths.items():
            clas_outdir = os.path.join(Configs.outdir, 'query_classifications',
                                       marker)
            taxonomy_path = refpkg[marker]['taxonomy']
            if not os.path.exists(taxonomy_path):
                taxonomy_path = os.path.join(
                    os.path.dirname(taxonomy_path), 'all_taxon.taxonomy')
            filtered_path = os.path.join(
                clas_outdir,
                f"placement.classification.{support_value.split('.')[-1]}")
            filterClassification(taxonomy_path, classification_path,
                                 filtered_path, float(support_value))
            filtered_paths[marker] = filtered_path

        # Aggregate all classifications into a single output file
        all_classification_path = os.path.join(
            Configs.outdir, 'query_classifications.tsv')
        header_written = False
        with open(all_classification_path, 'w') as f:
            for marker, filtered_path in filtered_paths.items():
                with open(filtered_path, 'r') as fptr:
                    lines = fptr.read().strip().split('\n')
                if not lines:
                    continue
                if not header_written:
                    f.write(lines[0] + '\n')
                    header_written = True
                if len(lines) > 1:
                    f.write('\n'.join(lines[1:]) + '\n')

    return classification_paths, filtered_paths

def updateAbundanceProfile(inpath, abundance_profile):
    """Update abundance profile counts from a filtered classification file."""
    with open(inpath, 'r') as f:
        lines = f.read().strip().split('\n')
        for i in range(1, len(lines)):
            parts = lines[i].split('\t')
            if len(parts) < 2:
                continue
            fields = [int(x) if x.upper() != 'NA' else 0 for x in parts[1:]]
            if sum(fields) == 0:
                continue
            for j in range(min(len(RANKS), len(fields))):
                abundance_profile[RANKS[j]][fields[j]] += 1

def parseTaxonomy(inpath):
    """Parse taxonomy file (all_taxon.taxonomy) into a taxid->(name, parent, rank) map.

    Uses csv.reader to correctly handle quoted fields (some taxon names
    contain commas).
    """
    taxid_map = {}
    with open(inpath, 'r', newline='') as f:
        reader = csv.reader(f)
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 4:
                continue
            taxid = int(row[0])
            parent_id = int(row[1])
            rank = row[2]
            taxname = row[3]
            taxid_map[taxid] = (taxname, parent_id, rank)
    return taxid_map


def parseSpeciesToMarker(inpath):
    """Parse species-to-marker map file (species_to_marker.tsv)."""
    species_to_marker = defaultdict(list)
    with open(inpath, 'r') as f:
        for line in f:
            if line.startswith('tax_id'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            species = int(parts[0])
            markers = parts[2].split(',')
            species_to_marker[species] = markers
    return species_to_marker

def reorderJson(inpath, outpath):
    """Reorder a .jplace file to a standard field order for classification."""
    with open(inpath, 'r') as fh:
        obj = json.load(fh)

    old_order = obj['fields']
    new_plc_order = ['edge_num', 'likelihood', 'like_weight_ratio',
                     'distal_length', 'pendant_length']

    newobj = {
        'tree': obj['tree'],
        'placements': [],
    }

    for placement in obj['placements']:
        if 'nm' in placement:
            ori_p, ori_n = placement['p'], [x[0] for x in placement['nm']]
        else:
            ori_p, ori_n = placement['p'], placement['n']
        new_p = []
        for pp in ori_p:
            pp_map = {old_order[i]: pp[i] for i in range(len(pp))}
            tmp = [pp_map[field] for field in new_plc_order]
            new_p.append(tmp)
        newobj['placements'].append({'p': new_p, 'n': ori_n})

    newobj['metadata'] = obj['metadata']
    newobj['version'] = obj['version']
    newobj['fields'] = new_plc_order

    with open(outpath, 'w') as ofh:
        json.dump(newobj, ofh, indent=4)

def getClassification(marker, taxonomy_path, mapping_path,
                      ori_jplace_path, jplace_path,
                      outdir, classification_path, lock):
    """Obtain taxonomic classification for a marker gene from placement results."""
    from tipp3.jplace_classifier import classify_jplace

    reorderJson(ori_jplace_path, jplace_path)

    if not (os.path.exists(classification_path)
            and os.stat(classification_path).st_size > 0):
        classify_jplace(
            jplace_path=jplace_path,
            taxonomy_path=taxonomy_path,
            mapping_path=mapping_path,
            classification_path=classification_path,
            threshold=0.0,
            cutoff=0.0,
        )

    if os.path.exists(jplace_path):
        os.remove(jplace_path)

    return marker, classification_path

def loadTaxonomy(taxonomy_file, lower=True):
    """Load full taxonomy with level maps for classification filtering.

    Uses csv.reader to correctly handle quoted fields (some taxon names
    contain commas).
    """
    with open(taxonomy_file, 'r', newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        header = [h.lower().strip() for h in header]
        key_map = {header[i]: i for i in range(len(header))}

        taxon_map = {}
        level_map = {level: {} for level in RANKS}

        for row in reader:
            if lower:
                row = [cell.lower() for cell in row]
            taxon_map[row[0]] = row

            for level in RANKS:
                if key_map.get(level) is None:
                    continue
                idx = key_map[level]
                if idx >= len(row):
                    continue
                val = row[idx]
                if val == '':
                    continue
                if val not in level_map[level]:
                    level_map[level][val] = {}
                level_map[level][val][row[0]] = row[0]

    return taxon_map, level_map, key_map

def filterClassification(taxonomy_path, classification_path, filtered_path,
                         threshold):
    """Filter classification output by a support value threshold."""
    taxon_map, level_map, key_map = loadTaxonomy(taxonomy_path)

    level_map_hierarchy = {
        "species": 0, "genus": 1, "family": 2, "order": 3,
        "class": 4, "phylum": 5, "superkingdom": 6, "root": 7
    }
    old_name, old_id, old_rank, old_probability = "", "", "", 1

    classification = {}
    with open(classification_path, 'r', newline='') as class_in:
        reader = csv.reader(class_in)
        for results in reader:
            if len(results) < 5:
                continue
            name, id_, taxname, rank, probability = (
                results[0], results[1], results[2],
                results[-2], float(results[-1]))

            if name != old_name:
                if old_name != "" and old_id in taxon_map:
                    lineage = taxon_map[old_id]
                    output_line = [old_name]
                    for level in RANKS:
                        clade = lineage[key_map[level]]
                        output_line.append(clade if clade != "" else "NA")
                    classification[old_name] = output_line
                old_name = name
                old_rank = "root"
                old_probability = 1
                old_id = '1'

            if (rank in level_map_hierarchy
                    and level_map_hierarchy[old_rank] > level_map_hierarchy[rank]
                    and probability > threshold):
                old_rank = rank
                old_probability = probability
                old_id = id_
            elif (rank in level_map_hierarchy
                  and level_map_hierarchy[old_rank] == level_map_hierarchy[rank]
                  and probability > old_probability):
                old_rank = rank
                old_probability = probability
                old_id = id_

    if old_name != "" and old_id in taxon_map:
        lineage = taxon_map[old_id]
        output_line = [old_name]
        for level in RANKS:
            clade = lineage[key_map[level]]
            output_line.append(clade if clade != "" else "NA")
        classification[old_name] = output_line

    with open(filtered_path, 'w') as class_out:
        class_out.write("fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n")
        for frag in sorted(classification.keys()):
            class_out.write("\t".join(classification[frag]) + "\n")
