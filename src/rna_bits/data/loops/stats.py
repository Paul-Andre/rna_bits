import sys, os
import csv
from Bio import PDB
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import numpy as np

STRUCT_DIR = "norm_representative/"
MCA_DIR = "/home/paul/MC-Annotate"


def get_conv(d):
    with open(d) as f:
        return dict(
            l.strip().split()
            for l in f
            if not l.startswith("#") and len(l.strip()) != 0
        )


res_dict = get_conv("data/residues.list")
atom_dict = get_conv("data/atoms.list")

struct_filenames = [a for a in os.listdir(STRUCT_DIR) if a.endswith(".pdb")]
struct_filenames.sort()
# struct_filenames = struct_filenames[:100]


def is_nuc(residue):
    if res_dict.get(residue.resname, "-") in "AUCG":
        return True
    else:
        return False


atom_cnt = Counter()
res_cnt = Counter()

failed = []
distances = []
unrec_res_cnt = Counter()
for sf in struct_filenames:
    try:
        m = PDBParser().get_structure(sf, os.path.join(STRUCT_DIR, sf))[0]
    except Exception as e:
        print("couldn't open", sf)
        print(traceback.format_exc())
        failed.append(sf)
        continue
    prev = None
    for residue in m.get_residues():
        if residue.resname not in res_dict:
            print("Not recognized residue", residue.resname, residue.full_id)
            unrec_res_cnt[residue.resname] += 1
        if is_nuc(residue):
            res_cnt[residue.resname] += 1
            for a in residue:
                if a.element == "H":
                    continue
                if a.name not in atom_dict:
                    print("not recognized", a.name, a.full_id, residue.resname)
                atom_cnt[atom_dict.get(a.name)] += 1

            #            if prev is None:
            #                pass
            #            elif "P" not in residue:
            #                print(residue.full_id, "does not have P")
            #            elif "O3'" not in prev:
            #                print(prev.full_id, "does not have O3'")
            #            else:
            #                dist = np.linalg.norm(prev["O3'"].coord - residue["P"].coord)
            #                distances.append(dist)
            #                if (dist > 1.8 and dist < 2.0) :
            #                    print(prev.full_id, residue.full_id, dist)

            prev = residue
        else:
            prev = None


# print("got", len(distances), "distances")
# print("Plotting histogram")
# import math
# print(distances[:10])
# for d in distances:
#    if math.isnan(d):
#        print("Nan")
# import seaborn as sns
# import matplotlib.pyplot as plt
# plt.hist(np.array(distances, dtype=float), bins=200, range=(1.4,1.8))
# plt.show()

for a, b in atom_cnt.items():
    print(a, b)
print("###################################################")
for a, b in res_cnt.items():
    print(a, b)

print("###################################################")
for a, b in unrec_res_cnt.items():
    print(a, b)
