"""
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
import pandas as pd

from rna_bits.utils.data_path import get_path
from rna_bits.utils.ss import parse_ss_file
from rna_bits.utils.ss import Segmenter
from rna_bits.utils.misc import remove_string_end
from rna_bits.utils.mc_annotate import query_mc_name, mc_name_to_tuple

from rna_bits.builder import builder

REFERENCE_DIR = get_path("benchmark/auto_from_pdb/all/insert_loops_force_native/rass")

def get(DIR):
    """
    DIR = "benchmark/auto_from_pdb/all/insert_loops_force_native"
    """
    RASS_DIR = get_path(pjoin(DIR, "rass"))
    OUT_DIR = get_path(pjoin(DIR, "rass_info"), create=True)

    print(RASS_DIR)
    filenames = [a for a in os.listdir(RASS_DIR) if os.path.isdir(pjoin(RASS_DIR, a))]
    filenames.sort()
    assert(filenames)

    for i, fn in enumerate(filenames):
        instance_fns = [a for a in os.listdir(os.path.join(RASS_DIR,fn)) if a.endswith(".rass")]
        instance_fns.sort()

        # TODO: this is a hack in order to avoid counting the amount of loops in the rass file
        try:
            with open(pjoin(REFERENCE_DIR, fn, "1.rass")) as f:
                ref_motifs = f.read().count("motif:")
        except FileNotFoundError:
            ref_motifs = float("NaN")

        for ii, fnn in enumerate(instance_fns):
            print(fnn, fn, str(ii + 1) + "/" + str(len(instance_fns)), str(i+1) + "/" + str(len(filenames)))

            #print(fn, str(i+1) + "/" + str(len(filenames)))
            nat = remove_string_end(fnn, ".rass")

            with open(pjoin(RASS_DIR, fn, fnn)) as f:
                targ_motifs = f.read().count("motif:")
    
            out = {}
            out["num_loops"] = ref_motifs
            out["num_loops_found"] = targ_motifs

            os.makedirs(pjoin(OUT_DIR, fn), exist_ok=True)

            with open(os.path.join(OUT_DIR, fn, nat+".json"), "w") as f:
                f.write(json.dumps(out))

if __name__ == "__main__":
    DIR = "benchmark/auto_from_pdb/all/insert_loops_top_1_exclude_native"
    get(DIR)
    DIR = "benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native"
    get(DIR)
    DIR = "benchmark/auto_from_pdb/all/insert_loops_force_native"
    get(DIR)
    DIR = "benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_all"
    get(DIR)
    DIR = "benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_reliable"
    get(DIR)


