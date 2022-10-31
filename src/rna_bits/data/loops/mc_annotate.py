import sys,os
import csv
from Bio import PDB
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import numpy as np
#from utils import *

import subprocess

STRUCT_DIR = "norm_representative/"
OUT_DIR = "mca_out/"
MCA_DIR = "/home/paul/MC-Annotate"

if not os.path.isdir(OUT_DIR):
    os.mkdir(OUT_DIR)

struct_filenames = [a for a in os.listdir(STRUCT_DIR) if a.endswith(".pdb")]
struct_filenames.sort()

for i,in_fn in enumerate(struct_filenames):
    print(in_fn, str(i+1)+"/"+str(len(struct_filenames)))
    fp = subprocess.run([MCA_DIR, "-f", "0", os.path.join(STRUCT_DIR, in_fn)], capture_output=True)

    with open(os.path.join(OUT_DIR, in_fn+".txt"), "wb") as f:
        f.write(fp.stdout)

