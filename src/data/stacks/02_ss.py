import utils.pdb
from Bio.PDB.PDBList import PDBList
from utils.data_path import get_path
import subprocess
import os
from collections import Counter
import sys
import warnings

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from utils import mc_annotate

STRUCT_DIR = get_path("interim/stacks/pdb")
MCOUT_DIR = get_path("interim/stacks/mcout")
SS_DIR = get_path("interim/stacks/ss", create=True)

mcout_filenames = [a for a in os.listdir(MCOUT_DIR) if a.endswith(".mcout")]
mcout_filenames.sort()
#struct_filenames = struct_filenames[:10]
#seen_pairings = {p:Counter() for p in possible_pairs}
assert mcout_filenames

def contains_O3p(residue):
    return "O3'" in residue or "O3*" in residue

def get_O3p(residue):
    if "O3'" in residue:
        return residue["O3'"]
    else:
        return residue["O3*"]

for i,fn in enumerate(mcout_filenames):
    print(fn, str(i+1)+"/"+str(len(mcout_filenames)))
    with open(os.path.join(MCOUT_DIR, fn)) as f:
        mc_out = mc_annotate.MCOut(f.read())

    base_code = fn[:-len(".mcout")]

    struct_fn = base_code +".pdb"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        m = PDBParser().get_structure(struct_fn, os.path.join(STRUCT_DIR, struct_fn))[0]

    all_chains = []
    current_chain = []

    prev = None
    for residue in m.get_residues():
        name = mc_annotate.get_mc_style_name(residue)
        if (name in mc_out.proper_nucs):

            continuing = False
            if (prev is not None) and contains_O3p(prev) and ("P" in residue):
                dist = np.linalg.norm(get_O3p(prev).coord - residue["P"].coord)
                if (dist <= 2.0):
                    continuing = True

            if continuing:
                current_chain.append(name)
            else:
                if len(current_chain) > 0:
                    all_chains.append(current_chain)
                    current_chain = []

                current_chain.append(name)

            prev = residue
        else:
            if len(current_chain) > 0:
                all_chains.append(current_chain)
                current_chain = []

            prev = None

    if len(current_chain) > 0:
        all_chains.append(current_chain)
        current_chain = []


    # Sanity checks:
    cnt = Counter()
    for c in all_chains:
        cnt.update(c)
    for a,b in cnt.items():
        assert(b==1)
    assert(set(cnt) == mc_out.proper_nucs)
    
    out_fn = os.path.join(SS_DIR, base_code+".ss")
    with open(out_fn, "w", encoding="ascii") as f:
        f.write("chains:\n")
        for c in all_chains:
            f.write(" ".join(c))
            f.write("\n")
        f.write("pairs:\n")
        pairings = sorted(set(tuple(sorted((mc_annotate.mc_name_to_tuple(a),mc_annotate.mc_name_to_tuple(b)))) for a,b in mc_out.proper_pairings.items() if b is not None))
        for a,b in pairings:
            a = mc_annotate.tuple_to_mc_name(a)
            b = mc_annotate.tuple_to_mc_name(b)
            f.write(a)
            f.write(" ")
            f.write(b)
            f.write("\n")


