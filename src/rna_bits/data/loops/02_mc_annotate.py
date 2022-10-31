import sys, os
import csv
import logging
import traceback
from collections import Counter
import numpy as np
import subprocess

from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from rna_bits.utils.data_path import get_path

STRUCT_DIR = get_path("interim/loops/norm_representative/")
OUT_DIR = get_path("interim/loops/mca_out/", create=True)

# TODO: Put this into rna_bits.utils
MCA_DIR = "/home/paul/MC-Annotate"

struct_filenames = [a for a in os.listdir(STRUCT_DIR) if a.endswith(".pdb")]
struct_filenames.sort()
assert struct_filenames

for i, in_fn in enumerate(struct_filenames):
    print(in_fn, str(i + 1) + "/" + str(len(struct_filenames)))
    fp = subprocess.run(
        [MCA_DIR, "-f", "0", os.path.join(STRUCT_DIR, in_fn)], capture_output=True
    )

    with open(os.path.join(OUT_DIR, in_fn + ".txt"), "wb") as f:
        f.write(fp.stdout)
