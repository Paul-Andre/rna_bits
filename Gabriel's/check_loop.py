import pickle
from collections import Counter, defaultdict
import os
import traceback # simply to print stack trace

import numpy as np

import rna_bits.utils.pdb

# From biopython package:
import Bio.PDB.mmtf


def fetch_residues_from_graph(g, pdb_code, chain_id):
    """
    Returns a dict int->residue
    """
    print("Fetching", pdb_code, chain_id)
    c = rna_bits.utils.pdb.fetch_pdb(pdb_code)[0][chain_id]
    # c = Bio.PDB.mmtf.MMTFParser.get_structure_from_url(pdb_code)[0][chain_id]
    print("Finished fetching", pdb_code, chain_id)

    letters = []

    residues_dict = {}

    for i, node in g._node.items():
        rid = node["fr3d"]

        print(rid)
        if not rid[-1].isnumeric():
            pos = int(rid[:-1])
            ins_code = rid[-1]
        else:
            pos = int(rid)
            ins_code = " "

        # TODO: Here I assume that the residue is never a hetero-residue
        r = c[(" ", pos, ins_code)]

        print(r.resname)
        print(node["real_nt"])
        assert(r.resname == node["real_nt"])

        letters.append(r.resname)

        residues_dict[i] = r

    return residues_dict


def compute_backbone_next(residues_dict):
    backbone_next = {}
    for i,ir in residues_dict.items():
        for j,jr in residues_dict.items():
            if i==j: continue

            i_oxygen = None
            if "O3'" in ir:
                i_oxygen = ir["O3'"]
            if "O3*" in ir:
                i_oxygen = ir["O3*"]

            j_phosphate = None
            if "P" in jr:
                j_phosphate = jr["P"]

            if i_oxygen is not None and j_phosphate is not None:
                dist = np.linalg.norm(i_oxygen.coord - j_phosphate.coord)
                if dist <= 2.0:
                    assert(i not in backbone_next), "Multiple bb connections??"
                    backbone_next[i] = j

    return backbone_next

def get_loop_order(g, backbone_next, considered_set=None):
    """
    g is an interaction network given as a networkx.classes.digraph.DiGraph 

    backbone_next is a dict[int,int] indicating the next nucleotide along the backbone

    considered_set, if given, is a set of integers representing the nodes in
    the graph we consider (i.e. pretend that the rest don't exist)

    Checks if the interaction network is a "loop", e.g. a cycle of only B53 and CWW interactions that
    covers all considered graph nodes

    If the graph forms such a cycle, returns a list of lists representing the "strands" of the cycle
    Otherwise, returns None
    """

    if considered_set is None:
        considered_set = g.nodes

    def follow_label(i, label):
        for j,edge_attributes in g.succ[i].items():
            if j in considered_set and edge_attributes["label"] == label:
                return j
        return None
    
    have_incoming = set()
    for i in considered_set:
        if i in backbone_next and backbone_next[i] in considered_set:
            have_incoming.add(backbone_next[i])

    strand_beginnings = set(considered_set) - have_incoming
    print(strand_beginnings, "strand_beginnings")

    if len (strand_beginnings) == 0:
        return None

    first = min(strand_beginnings)
    ret = [[]]

    vis = set()

    current = first
    vis.add(current)
    ret[-1].append(current)

    prev = None
    while True:
        next = backbone_next.get(current)
        if next is None or next not in considered_set:
            next = follow_label(current, "CWW")
            if next is None or next is prev:
                # Couldn't find a B53 or CWW edge
                print("no next")
                return None
            if next == first:
                break;
            ret.append([])
        if next in vis:
            # if we have already visited next, and it's not the first node,
            # what we have is not a cycle that covers the whole graph
            print(next)
            print("returneeedd")
            return None
        prev = current
        current = next
        vis.add(current)
        ret[-1].append(current)

    if set(vis) != set(considered_set):
        print("not full visit")
        return None

    return ret



rins_by_counts = defaultdict(list)

def process(k,v):

    print()
    print()

    print("Processing rin", k)
    g = v["graphs"]
    (p_id, positions) = v["positions"][0]
    (s_id, sequence) = v["seqs"][0]
    assert(p_id == s_id)

    pdb_code, chain_id = p_id.split(".")

    residues_dict = fetch_residues_from_graph(g, pdb_code, chain_id)

    # Sanity check, make sure that the sequence actually matches the pdb file
    letters = []
    for (_, r) in sorted(residues_dict.items()):
        letters.append(r.resname)
    # if (sequence != letters) :
    #     print("Failed", sequence, letters)
    assert(sequence == letters)

    backbone_next = compute_backbone_next(residues_dict)
    print(backbone_next)

    loop_order = get_loop_order(g, backbone_next)
    print(loop_order)

    if loop_order is not None:
        print("RIN", k, "is a proper loop")
        print(loop_order)

        # Save file, using my utility code for ease
        # parent_dir = "motifs_from_graph_but_only_cycles_take_2/"+str(k)
        # os.makedirs(parent_dir, exist_ok=True)
        # save_file = parent_dir+"/"+p_id+".pdb"

        # model = rna_bits.utils.pdb.build_model_from_lists_of_residues([[r for (_,r) in sorted(residues_dict.items())]])
        # rna_bits.utils.pdb.save_model_as_pdb(model, save_file)
        # print("saved", save_file)


        # Display the interactions as reported in the data:
        cww_count =0
        for sk, sv in g.succ.items():
            for svk, svv in sv.items():
                label = svv["label"]
                if label == "b53" or sk<svk:
                    if label == "CWW":
                        cww_count+=1
                    print(sk, svk, label)
        print("cww_count", cww_count)

        rins_by_counts[cww_count].append(k)


    return loop_order


f = open("3drins_by_id.pickle", "rb")
p = pickle.load(f)


failed_list = []
loop_rins = []
for k, v in sorted(p.items()):
    # if (k % 10 != 0):
    #     continue
    proper_loops= [68, 72, 96, 102, 103, 116, 127, 131, 140, 142, 143, 153, 154, 155, 179, 183, 186, 189, 191, 192, 234, 252, 288, 289, 300, 301, 326, 336, 356]

    failed=[178, 282, 295]

    if k not in failed:
        continue

    try:
        loop_order = process(k,v)
        if loop_order is not None:
            loop_rins.append(k)
    except Exception as e:
        print("FAILED", k)
        print(e)
        traceback.print_exc()
        failed_list.append((k, e))


print()
print(rins_by_counts)
print("proper loops:", loop_rins)
print("rins that failed to process", failed_list)
