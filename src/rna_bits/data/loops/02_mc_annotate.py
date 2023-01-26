import sys, os
import csv
import logging
import traceback
from collections import Counter
import numpy as np

import shutil

from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from rna_bits.utils.data_path import get_path
from rna_bits.utils.mc_annotate import run_mc_annotate

ORIG_STRUCT_DIR = get_path("interim/loops/representative/")
STRUCT_DIR = get_path("interim/loops/norm_representative/")
OUT_DIR = get_path("interim/loops/mca_out/", create=True)

with open(ORIG_STRUCT_DIR + "/version.txt") as f:
    version = f.read().strip()

PROVIDED_DIR = get_path("provided/mca_out");

# Not actually an archive, just a directory
archive_fname = PROVIDED_DIR + "/" + version
if os.path.isdir(archive_fname):
    print("Instead of running MC-Annotate using provided data from",
            archive_fname)
    for fname in os.listdir(archive_fname):
        shutil.copy(archive_fname+"/"+fname, OUT_DIR)
else:
    struct_filenames = [a for a in os.listdir(STRUCT_DIR) if a.endswith(".pdb")]
    struct_filenames.sort()
    assert struct_filenames

    for i, in_fn in enumerate(struct_filenames):
        print("running MC-Annotate", in_fn, str(i + 1) + "/" + str(len(struct_filenames)))

        mcout = run_mc_annotate( os.path.join(STRUCT_DIR, in_fn))

        with open(os.path.join(OUT_DIR, in_fn + ".txt"), "wb") as f:
            f.write(mcout)

    
