import sys, os
import csv
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import subprocess

import numpy as np
from Bio import PDB

from rna_bits.utils.data_path import get_path
from rna_bits.utils.misc import remove_string_end

RASS_DIR = get_path("benchmark/auto_from_pdb/all/provided_ss")


def run_rna_insert_loop(RASS_DIR, OUT_DIR, command):
    struct_filenames = [a for a in os.listdir(RASS_DIR) if a.endswith(".rass")]
    struct_filenames.sort()

    for i, fn in enumerate(struct_filenames):
        print(fn, str(i + 1) + "/" + str(len(struct_filenames)))
        fp = subprocess.run(["rna_insert_loops"] + command + ["--out_dir", os.path.join(OUT_DIR, remove_string_end(fn, ".rass")),
            os.path.join(RASS_DIR, fn)])
        
    #OUT_DIR_PB2 = get_path("benchmark/auto_from_pdb/bp2_limited/provided_ss", create=True)

if __name__ == "__main__":
    OUT_DIR = get_path("benchmark/auto_from_pdb/all/insert_loops_top_1_exclude_native/rass", create=True)
    run_rna_insert_loop(RASS_DIR, OUT_DIR, [ "--top", "1", "--samples", "1", "--exclude_native"])

    OUT_DIR = get_path("benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/rass", create=True)
    run_rna_insert_loop(RASS_DIR, OUT_DIR, [ "--top", "10", "--samples", "50", "--exclude_native"])
