import sys, os
import csv
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import subprocess
import os
import shutil

from rna_bits.utils.data_path import get_path
from rna_bits.utils.misc import remove_string_end


def run_builder(DIR):

    RASS_DIR = get_path(DIR+"rass/")
    PDB_DIR = get_path(DIR+"/generated/", create=True)

    rass_dirs = [a for a in os.listdir(RASS_DIR) if os.path.dirname(os.path.join(RASS_DIR,a))]
    print(RASS_DIR)
    print(rass_dirs)
    assert(len(rass_dirs))
    rass_dirs.sort()

    for i, fn in enumerate(rass_dirs):
        struct_filenames = [a for a in os.listdir(os.path.join(RASS_DIR,fn)) if a.endswith(".rass")]
        struct_filenames.sort()

        for ii, fnn in enumerate(struct_filenames):
            os.makedirs(os.path.join(PDB_DIR, fn), exist_ok=True)

            print(fnn, fn, str(ii + 1) + "/" + str(len(struct_filenames)), str(i+1) + "/" + str(len(rass_dirs)))

            # TODO: add an --out parameter to rna_builder so don't need to move
            command = ["rna_builder", os.path.join(RASS_DIR, fn, fnn)]
            print(command)
            fp = subprocess.run(command)
            shutil.move(os.path.join(RASS_DIR, fn, fnn +".pdb"), os.path.join(PDB_DIR, fn, fnn+".pdb"))


if __name__ == "__main__":
    DIR = "benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/"
    run_builder(DIR)

    DIR = "benchmark/auto_from_pdb/all/insert_loops_top_1_exclude_native/"
    run_builder(DIR)



