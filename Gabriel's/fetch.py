import pickle
from collections import Counter
import rna_bits.utils.pdb as rna_pdb
import os

f = open("3drins_by_id.pickle", "rb")
p = pickle.load(f)

already_saved = set()
inteacts = Counter()
for k, v in sorted(p.items()):
    parent = "motifs/"+str(k)
    os.makedirs(parent, exist_ok=True)
    for pi,((pk, pv), (sk, sv))  in enumerate(zip(v["positions"], v["seqs"])):
        if pi == 0:
            continue
        
        assert(pk == sk)
        code, chain_id = pk.split(".")
        print("fetching", code)
        c = rna_pdb.fetch_pdb(code)[0][chain_id]
        print("fetched")
        collect = []
        letters = []
        orig_ids = []
        i = 1

        if False:
            for rid in pv:
                r = c[rid]
                collect.append(r)
                letters.append(r.resname)
                orig_ids.append(r.id)
                i+=1

        else:
            for r in c:
                print(i, r)
                if True:
                    # if r.id[0]=="w":  # Not a hetero residue
                    #     continue
                    if r.id[0]=="w":  # Not a hetero residue
                        continue

                if i in pv:
                    collect.append(r)
                    letters.append(r.resname)
                    orig_ids.append(r.id)
                i+=1

        
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
