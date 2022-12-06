import pickle
from collections import Counter
import rna_bits.utils.pdb as rna_pdb
import os

f = open("3drins_by_id.pickle", "rb")
p = pickle.load(f)

already_saved = set()
inteacts = Counter()
nn = 0


def fetch(k,v, parent):
    pi = 0
    (pk, pv) = v["positions"][0]
    (sk, sv) = v["seqs"][0]
    g = v["graphs"]

    
    assert(pk == sk)
    code, chain_id = pk.split(".")
    print("fetching", code)
    c = rna_pdb.fetch_pdb(code)[0][chain_id]
    print("fetched")
    collect = []
    letters = []
    orig_ids = []

    assert(len(pv) == len(sv))
    assert(len(pv) == len(g._node))
    for i in range(len(pv)):
        ii = i+1

        rid = g._node[ii]["fr3d"]

        if not rid[-1].isnumeric():
            pos = int(rid[:-1])
            ins_code = rid[-1]
        else:
            pos = int(rid)
            ins_code = " "

        # TODO: Here I assume that the residue is never a hetero-residue
        r = c[(" ", pos, ins_code)]
        collect.append(r)
        letters.append(r.resname)
        orig_ids.append(r.id)

    
    model = rna_pdb.build_model_from_lists_of_residues([collect])

    save_file = parent+"/"+pk+"_"+str(pi)+".pdb"  # pv[0] = Me being lazy and not using a proper increment
    print(save_file)
    assert(save_file not in already_saved)
    rna_pdb.save_model_as_pdb(model, save_file)
    already_saved.add(save_file)

    if (sv != letters) :
        print("Failed", sv, letters)
        print("      ",pv, orig_ids)
    assert(sv == letters)


def get_loop_order(g, considered_set=None):
    """
    Interprets the motif as a cycle of only B53 and CWW interactions that
    covers all graph nodes (or considered_set if provided).

    Returns None if it cannot.
    Otherwise returns a list of lists representing the "strands" of the cycle
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
        j = follow_label("B53")
        if j is not None:
            have_incoming.add(j)


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

    while True:
        next = follow_label(current, "B53")
        if next is None:
            next = follow_label(current, "CWW")
            if next is None:
                # Couldn't find a B53 or CWW edge
                return None
            if next == first:
                break;
            ret.append([])
        if next in vis:
            # if we have already visited next, and it's not the first node,
            # what we have is not a cycle
            return None
        current = next
        vis.add(current)
        ret[-1].append(current)

    if vis != considered_set:
        return None

    return ret

failed_list = []
for k, v in sorted(p.items()):
    parent = "motifs_from_graph/"+str(k)
    os.makedirs(parent, exist_ok=True)

    g = v["graphs"]
    if is_loop_motif(g, 

    try:
        fetch(k,v,parent)
    except Exception as e:
        
        failed_list.append((k, e))
        print("FAILED", k)
        print(e)

print(failed_list)





    # for pi,((pk, pv), (sk, sv))  in enumerate(zip(v["positions"], v["seqs"])):



#     for gk, gv, in  v["graphs"].adj.items():
#         for ggk, ggv in gv.items():
#             l = ggv["label"]
#             inteacts[l]+=1

# for k,v in inteacts.items():
#     print(k,v)

# rna_pdb.fetch_pdb("1hmh")
# a = rna_pdb.fetch_pdb("1hmh")
# a[0]["A"]
# c = a[0]["A"]
# for i, r in c:
#     print(i, r)
# for i, r in enumerate(c)
# :
#     print(i, r)
# for i, r in enumerate(c):
#     print(i, r)
# for i, r in enumerate(c):
#     print(i+1
#     , r)
# p[24]
# for i, r in enumerate(c):
#     if i in [12,13,28,29]:
#         print(i+1, r)
# for i, r in enumerate(c):
#     if i+1 in [12,13,28,29]:
#         print(i+1, r)
# p[24]["graphs"]
# g=p[24]["graphs"]
# g.adj.items()
# for k,v in g.adj.items():
#     print(k,v)
# %history
