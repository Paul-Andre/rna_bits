import sys,os
import csv
from Bio import PDB
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
from collections import defaultdict
import numpy as np
from utils import *

import json

STRUCT_DIR = "norm_representative/"
SEGMENTATION_DIR = "segmentations/"
MODEL_DIR = "loop_models/"
JSON_DIR = "json/"

if not os.path.isdir(MODEL_DIR):
    os.mkdir(MODEL_DIR)

if not os.path.isdir(JSON_DIR):
    os.mkdir(JSON_DIR)

filenames = [a for a in os.listdir(SEGMENTATION_DIR) if a.endswith(".txt")]
filenames.sort()
#filenames= filenames[:10]

def query_mc_name(model, name):
    (chain_id, res_id, insertion_code) = mc_name_to_tuple(name)
    return model[chain_id][(" ",res_id, insertion_code)]

pattern_cnt = Counter()

all_json_objs = defaultdict(list)

for i,seg_fn in enumerate(filenames):
    print(seg_fn, str(i+1)+"/"+str(len(filenames)))

    loop_segments = []

    with open(os.path.join(SEGMENTATION_DIR, seg_fn)) as f:
        i = 0
        for line in f:
            a = line.strip().split()
            pattern = a[0]
            nucs = a[1:]
            loop_segments.append(( nucs, pattern))

    if len(loop_segments) == 0:
        continue

    struct_fn = remove_string_end(seg_fn, ".txt")
    m = PDBParser().get_structure(struct_fn, os.path.join(STRUCT_DIR, struct_fn))[0]

    out_dir = os.path.join(MODEL_DIR, struct_fn)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)


    for loop_id, (nucs, pattern) in enumerate(loop_segments):
        loop_id += 1


        out_struct = PDB.Structure.Structure(loop_id)

        out_model = PDB.Model.Model(0)
        out_struct.add(out_model)

        # print(nucs, pattern)
        assert(len(nucs) == len(pattern))
        new_chain_id = "A"
        new_residue_id = 1
        out_chain = PDB.Chain.Chain(new_chain_id)
        out_model.add(out_chain)

        seqs = []
        current_seq = ""
        for i, (nuc, c) in enumerate(zip(nucs, pattern)):
            if c == ")" and i!=len(pattern)-1:
                new_chain_id = chr(ord(new_chain_id)+1)
                new_residue_id = 1
                out_chain = PDB.Chain.Chain(new_chain_id)
                out_model.add(out_chain)
                seqs.append(current_seq)
                current_seq=""
            out_residue = query_mc_name(m, nuc).copy()
            out_residue.id = (" ", new_residue_id, " ")
            new_residue_id+=1
            out_chain.add(out_residue)
            resname = out_residue.resname
            assert(resname in "CGAU")
            current_seq+=resname

        seqs.append(current_seq)

        simplified_pattern = ""
        for c in pattern:
            if c in "()":
                simplified_pattern+=c
            else:
                simplified_pattern+="."

        pattern_cnt[tuple(map(len, seqs))]+=1



        io = PDB.PDBIO()
        io.set_structure(out_struct)
        io.save(os.path.join(out_dir, str(loop_id)+".pdb"))

        # TODO: add "original source" in BGSU 3f3d|1|Xx|A|123|||A format
        json_obj = {
                "file": os.path.join(out_dir, str(loop_id)+".pdb"),
                "full_ss": pattern,
                "seqs": seqs,
                "original_nucs": nucs
                };

        all_json_objs["_".join(str(len(a)) for a in seqs)].append(json_obj)


for _,comps, freq in reversed(sorted((len(comps), comps, freq) for (comps, freq) in pattern_cnt.items())):
    pass
    #print(comps, freq)

cnt = Counter()
for (comps, freq) in pattern_cnt.items():
    cnt[len(comps)]+=freq

#dict("simple_ss":, "full_ss":, "strand_lengths", "file", "origin" I
for a,b in sorted(cnt.items()):
    print(a,b)

for k,v in all_json_objs.items():
    with open(os.path.join(JSON_DIR, k+".json"),"w") as f:
        f.write(json.dumps(v, indent=2))
