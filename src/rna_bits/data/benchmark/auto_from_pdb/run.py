from sys import argv
import sys, os
import csv
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import subprocess
from os.path import join  as pjoin
import shutil

from rna_bits.utils.data_path import get_path
from rna_bits.data.benchmark.auto_from_pdb.select_simple_strutures import  generate_structure_files
from rna_bits.data.benchmark.auto_from_pdb.insert_loops import  run_rna_insert_loop
from rna_bits.data.benchmark.auto_from_pdb.insert_loops_force_native import  run_insert_loops_force_native
from rna_bits.data.benchmark.auto_from_pdb.build import run_builder
from rna_bits.data.benchmark.auto_from_pdb.evaluate import run_evaluate
from rna_bits.data.benchmark.auto_from_pdb.run_bp2_bridge import run_bp2_bridge

# OUT_DIR_ALL = get_path("benchmark/auto_from_pdb/all/provided_ss", create=True)
# generate_structure_files(OUT_DIR_ALL, allow_multi_chain=False, allow_junctions=True)

# OUT_DIR_PB2 = get_path("benchmark/auto_from_pdb/bp2_limited/provided_ss", create=True)
# generate_structure_files(OUT_DIR_PB2, allow_multi_chain=False, allow_junctions=False)


# Some one-time functions that did things or something idk
def build_and_evaluate(DIR):
    """
    For example,
    DIR = "benchmark/auto_from_pdb/all/insert_loops_top_1_exclude_native/"
    """
    run_builder(DIR)
    run_evaluate(DIR)

def bp2b_something(DATASET_NAME):
    """
    For example,
    DATASET_NAME = "ALL"
    """
    DIR = "benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_"+DATASET_NAME.lower() + "/"
    RASS_DIR =  get_path(pjoin(DIR, "../provided_ss/"))
    BP2_RESULT_DIR = get_path(pjoin(DIR, "bp2_output/"))
    OUT_RASS_DIR =   get_path(pjoin(DIR, "rass/"), create=True)
    run_bp2_bridge(RASS_DIR, BP2_RESULT_DIR, OUT_RASS_DIR, DATASET_NAME, extra_args=["--samples", "50"]);

    build_and_evaluate(DIR)


if __name__ == "__main__":
    if len(argv) >= 3:
        if argv[1] == "build_and_evaluate":
            build_and_evaluate(argv[2])
        if argv[1] == "bp2b_something":
            bp2b_something(argv[2])
