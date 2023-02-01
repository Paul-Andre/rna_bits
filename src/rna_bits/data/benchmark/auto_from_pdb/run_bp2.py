import sys, os
import csv
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import subprocess
from os.path import join  as pjoin
import shutil

import numpy as np
from Bio import PDB

from rna_bits.utils.data_path import get_path
from rna_bits.utils.misc import remove_string_end



def run_bp2(RASS_DIR, OUT_DIR):
    struct_filenames = [a for a in os.listdir(RASS_DIR) if a.endswith(".rass")]
    struct_filenames.sort()


    for i, fn in enumerate(struct_filenames):
        print(fn, str(i + 1) + "/" + str(len(struct_filenames)))
        nat = remove_string_end(fn, ".rass")

        with open(pjoin(RASS_DIR, fn)) as in_f:
            seq = in_f.readline().strip()
            dot_bracket = in_f.readline().strip()
            rest = in_f.read()

        # fp = subprocess.run(["rna_insert_loops"] + command + ["--out_dir", os.path.join(OUT_DIR, remove_string_end(fn, ".rass")),
        #     os.path.join(RASS_DIR, fn)])

        fp = subprocess.run(["BayesPairing", "-seq", seq, "-ss", dot_bracket, "-d", "RELIABLE"])

        shutil.move("output.json", pjoin(OUT_DIR, nat +".json"))
        

if __name__ == "__main__":
    RASS_DIR = get_path("benchmark/auto_from_pdb/all/provided_ss")
    OUT_DIR = get_path("benchmark/auto_from_pdb/all/bp2_with_ss_reliable/bp2_output", create=True)
    run_bp2(RASS_DIR, OUT_DIR);

    # OUT_DIR = get_path("benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/rass", create=True)
    # run_rna_insert_loop(RASS_DIR, OUT_DIR, [ "--top", "10", "--samples", "50", "--exclude_native"])
