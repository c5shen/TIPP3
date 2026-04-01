"""Native Python replacement for the Java TIPPJsonMerger classification tool.

The Java tool used O(n^2) bipartition-string comparisons to map edges between
trees; for large reference trees (50K+ leaves) this made the classification step
take minutes or hours. This module maps placement edges to taxonomy nodes via
bottom-up LCA propagation in O(n) time, completing in under a second.
"""

import csv
import json
from collections import defaultdict

from tipp3 import get_logger

_LOG = get_logger(__name__)


def classify_jplace(jplace_path, taxonomy_path, mapping_path,
                    classification_path, threshold=0.0, cutoff=0.0):
    """Classify placed reads from a .jplace file using a taxonomy.

    Parameters
    ----------
    jplace_path : str
        Path to the (reordered) .jplace file.
    taxonomy_path : str
        CSV file: tax_id,parent_id,rank,tax_name,...
    mapping_path : str
        CSV file: seqname,tax_id
    classification_path : str
        Output file for classification results.
    threshold : float
        Minimum probability for a lineage to be reported.
    cutoff : float
        Minimum LWR for a single placement record to be included.
    """
    with open(jplace_path) as fh:
        jplace = json.load(fh)

    fields = jplace["fields"]
    placements = jplace["placements"]

    edge_idx = fields.index("edge_num")
    lwr_idx = fields.index("like_weight_ratio")

    seq_to_taxid = _read_mapping(mapping_path)
    taxonomy = _read_taxonomy(taxonomy_path)
    lca_fn = _build_lca_function(taxonomy)
    edge_to_taxnode = _parse_tree_and_map(jplace["tree"], seq_to_taxid,
                                          taxonomy, lca_fn)

    _LOG.info(f"Mapped {len(edge_to_taxnode)} edges to taxonomy nodes")

    with open(classification_path, "w") as out:
        for placement in placements:
            records = placement["p"]
            names = placement.get("n") or [nm[0] for nm in placement["nm"]]

            lwr_sum = sum(r[lwr_idx] for r in records if r[lwr_idx] > cutoff)
            if lwr_sum <= 0:
                continue

            lineage_prob = defaultdict(float)
            for record in records:
                lwr = record[lwr_idx]
                if lwr <= cutoff:
                    continue
                edge_num = str(int(record[edge_idx]))
                taxnode = edge_to_taxnode.get(edge_num)
                if taxnode is None:
                    continue

                prob = lwr / lwr_sum
                node_id = taxnode
                while node_id is not None:
                    lineage_prob[node_id] += prob
                    parent = taxonomy[node_id]["parent"]
                    node_id = parent if parent != node_id else None

            for fragment in names:
                for tax_id, prob in lineage_prob.items():
                    if prob >= threshold:
                        node = taxonomy[tax_id]
                        name = node["name"]
                        if "," in name:
                            name = f'"{name}"'
                        out.write(
                            f"{fragment},{tax_id},{name},"
                            f"{node['rank']},{prob:.4f}\n"
                        )


def _read_mapping(path):
    """Read seqname -> tax_id mapping file."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("seqname"):
                continue
            parts = line.strip().split(",")
            if len(parts) >= 2:
                mapping[parts[0]] = int(parts[1])
    return mapping


def _read_taxonomy(path):
    """Read taxonomy CSV into a dict: tax_id -> {parent, rank, name}.

    Uses csv.reader to correctly handle quoted fields (some taxon names
    contain commas, e.g. 'Salmonella enterica serovar 6,7:-:1,5').
    """
    taxonomy = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 4:
                continue
            tax_id = int(row[0])
            parent_id = int(row[1]) if row[1] else tax_id
            rank = row[2]
            name = row[3]
            taxonomy[tax_id] = {
                "parent": parent_id if parent_id != tax_id else None,
                "rank": rank,
                "name": name,
            }
    return taxonomy


def _build_lca_function(taxonomy):
    """Build an efficient LCA function using depth and ancestor lookup.

    Pre-computes depth for each node, then LCA of two nodes is found by
    walking from deeper to shallower until they meet.
    """
    depth = {}

    def _get_depth(tax_id):
        if tax_id in depth:
            return depth[tax_id]
        path = []
        node = tax_id
        while node is not None and node not in depth:
            path.append(node)
            node = taxonomy.get(node, {}).get("parent")

        base = depth.get(node, 0) if node is not None else 0
        for i, n in enumerate(reversed(path)):
            depth[n] = base + i + (1 if node is not None else 0)
        return depth.get(tax_id, 0)

    for tid in taxonomy:
        _get_depth(tid)

    def lca(a, b):
        if a is None:
            return b
        if b is None:
            return a
        da, db = depth.get(a, 0), depth.get(b, 0)
        while da > db:
            a = taxonomy.get(a, {}).get("parent")
            da -= 1
        while db > da:
            b = taxonomy.get(b, {}).get("parent")
            db -= 1
        while a != b and a is not None and b is not None:
            a = taxonomy.get(a, {}).get("parent")
            b = taxonomy.get(b, {}).get("parent")
        return a

    return lca


def _parse_tree_and_map(tree_str, seq_to_taxid, taxonomy, lca_fn):
    """Parse Newick tree and compute edge -> taxonomy LCA in a single pass.

    Instead of collecting all descendant leaves per edge and then computing
    LCA over potentially thousands of leaves, we propagate LCA bottom-up:
    a parent's LCA is lca(child1_lca, child2_lca, ...).  This means exactly
    one LCA call per internal node.
    """
    tree_str = tree_str.strip().rstrip(";")
    edge_to_taxnode = {}
    stack = []  # stack of LCA-so-far for current parenthetical group
    pos = 0
    n = len(tree_str)

    while pos < n:
        ch = tree_str[pos]

        if ch == "(":
            stack.append(None)
            pos += 1

        elif ch == ",":
            pos += 1

        elif ch == ")":
            pos += 1
            current_lca = stack.pop()
            label, edge_num, pos = _parse_node_annotation(tree_str, pos, n)

            if edge_num is not None and current_lca is not None:
                edge_to_taxnode[edge_num] = current_lca

            if stack:
                stack[-1] = lca_fn(stack[-1], current_lca)

        else:
            label, edge_num, pos = _parse_node_annotation(tree_str, pos, n)
            leaf_name = label.split(":")[0].strip("'\"") if label else None

            leaf_lca = None
            if leaf_name:
                tid = seq_to_taxid.get(leaf_name)
                if tid is not None and tid in taxonomy:
                    leaf_lca = tid

            if edge_num is not None and leaf_lca is not None:
                edge_to_taxnode[edge_num] = leaf_lca

            if stack:
                stack[-1] = lca_fn(stack[-1], leaf_lca)

    return edge_to_taxnode


def _parse_node_annotation(tree_str, pos, n):
    """Parse label:length{edge_num} starting at pos.

    Returns (label_str, edge_num_str_or_None, new_pos).
    """
    label_parts = []
    edge_num = None

    while pos < n and tree_str[pos] not in "(),;":
        if tree_str[pos] == "{":
            brace_end = tree_str.index("}", pos + 1)
            edge_num = tree_str[pos + 1:brace_end]
            pos = brace_end + 1
        elif tree_str[pos] == "[":
            bracket_end = tree_str.index("]", pos + 1)
            edge_num = tree_str[pos + 1:bracket_end]
            pos = bracket_end + 1
        else:
            label_parts.append(tree_str[pos])
            pos += 1

    return "".join(label_parts), edge_num, pos
