""" info such as num_nucs and original pdb_chains for the test cases
"""
import sys,os
import csv
import logging
import traceback
from collections import Counter
from os.path import join as pjoin
import json

import numpy as np
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from rna_bits.utils.data_path import get_path
from rna_bits.utils.ss import parse_ss_file
from rna_bits.utils.ss import Segmenter
from rna_bits.utils.misc import remove_string_end
from rna_bits.utils.mc_annotate import query_mc_name, mc_name_to_tuple


SS_DIR = get_path("interim/loops/ss_annotation/")
STRUCT_DIR = get_path("interim/loops/norm_representative/")

LIST_FILE = pjoin(get_path("provided"), "nrlist_3.267_4.0A.csv")

def generate_info(DIR):
    # TODO: hardcoded LIST_FILE

    with open(LIST_FILE, newline="") as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)

    pdb_by_nrlist = {}
    for (name, rep_s, _) in rows:
        pdb_by_nrlist[name] = rep_s

    PROVIDED_DIR = get_path(pjoin(DIR, "provided_ss"))
    OUT_DIR = get_path(pjoin(DIR, "info"), create=True)

    filenames = [a for a in os.listdir(PROVIDED_DIR) if a.endswith(".ref.index_noloose")]
    filenames.sort()
    assert(filenames)

    for i, fn in enumerate(filenames):
        print(fn, str(i+1) + "/" + str(len(filenames)))

        nat = remove_string_end(fn, ".ref.index_noloose")

        with open(pjoin(PROVIDED_DIR, nat+".cons.index_noloose")) as f:
            _,_,num = f.read().strip().split(":")
            num = int(num)

        out = {}
        out["pdb_chains"] = pdb_by_nrlist[nat]
        out["num_nucs"] = num

        with open(os.path.join(OUT_DIR, nat+".json"), "w") as f:
            f.write(json.dumps(out))


if __name__ == "__main__":
    DIR = "benchmark/auto_from_pdb/all/"
    generate_info(DIR);
    DIR = "benchmark/auto_from_pdb/bp2_limited/"
    generate_info(DIR);

