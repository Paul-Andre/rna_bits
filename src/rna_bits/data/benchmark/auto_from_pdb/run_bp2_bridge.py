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

# TODO: somehow, either by changing BP2 code myself, or by doing some kind of
# hack, make it so that this is found automatically
# XXX: Need to set this to the right directory on each computer.
BP2_DIR = "/home/paul/Masters/installation_try/rnabayespairing2/"

BP2_MODELS_DIR = pjoin(BP2_DIR, "bayespairing/models/")

def run_bp2_bridge(RASS_DIR, BP2_RESULT_DIR, OUT_RASS_DIR, DATASET_NAME, extra_args=[]):
    struct_filenames = [a for a in os.listdir(RASS_DIR) if a.endswith(".rass")]
    struct_filenames.sort()


    for i, fn in enumerate(struct_filenames):
        print(fn, str(i + 1) + "/" + str(len(struct_filenames)))
        nat = remove_string_end(fn, ".rass")

        with open(pjoin(RASS_DIR, fn)) as in_f:
            seq = in_f.readline().strip()
            dot_bracket = in_f.readline().strip()
            #rest = in_f.read()

        # fp = subprocess.run(["rna_insert_loops"] + command + ["--out_dir", os.path.join(OUT_DIR, remove_string_end(fn, ".rass")),
        #     os.path.join(RASS_DIR, fn)])


        os.makedirs(pjoin(OUT_RASS_DIR, nat), exist_ok=True)

        fp = subprocess.run(["bp2_bridge",
            # "-seq", seq,
            "-ss", dot_bracket,
            "-d", pjoin(BP2_MODELS_DIR, DATASET_NAME + ".json"), "-r", pjoin(BP2_RESULT_DIR, nat +".json"),
            "-o", pjoin(OUT_RASS_DIR, nat, "out"),
            "-m", get_path("benchmark/bp2b_motifs/" + DATASET_NAME, create=True),
            "--output_absolute_path" 
            ] + extra_args)

        

if __name__ == "__main__":
    # RASS_DIR = get_path("benchmark/auto_from_pdb/bp2_limited/provided_ss")
    # BP2_RESULT_DIR = get_path("benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_reliable/bp2_output")
    # OUT_RASS_DIR = get_path("benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_reliable/rass", create=True)
    # run_bp2_bridge(RASS_DIR, BP2_RESULT_DIR, OUT_RASS_DIR, "RELIABLE");

    RASS_DIR = get_path("benchmark/auto_from_pdb/bp2_limited/provided_ss")
    BP2_RESULT_DIR = get_path("benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_all/bp2_output")
    OUT_RASS_DIR = get_path("benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_all/rass", create=True)
    run_bp2_bridge(RASS_DIR, BP2_RESULT_DIR, OUT_RASS_DIR, "ALL");

# if __name__ == "__main__":
#     OUT_DIR = get_path("benchmark/auto_from_pdb/all/insert_loops_top_1_exclude_native/rass", create=True)
#     run_rna_insert_loop(RASS_DIR, OUT_DIR, [ "--top", "1", "--samples", "1", "--exclude_native"])

#     OUT_DIR = get_path("benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/rass", create=True)
#     run_rna_insert_loop(RASS_DIR, OUT_DIR, [ "--top", "10", "--samples", "50", "--exclude_native"])

    # OUT_DIR = get_path("benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/rass", create=True)
    # run_rna_insert_loop(RASS_DIR, OUT_DIR, [ "--top", "10", "--samples", "50", "--exclude_native"])
