import sys, os
import csv
import logging
import traceback
from collections import Counter
import numpy as np

from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from rna_bits.utils.data_path import get_path
from rna_bits.utils.normalize import Normalizer

STRUCT_DIR = get_path("interim/loops/representative/")
SAVE_DIR = get_path("interim/loops/norm_representative/", create=True)


struct_filenames = [a for a in os.listdir(STRUCT_DIR) if a.endswith(".pdb")]
struct_filenames.sort()
# struct_filenames = struct_filenames[:10]
assert struct_filenames

failed = []

for i, sf in enumerate(struct_filenames):
    print(sf, str(i + 1) + "/" + str(len(struct_filenames)))
    s = PDBParser().get_structure(sf, os.path.join(STRUCT_DIR, sf))
    out_s = Normalizer(convert_hetero_nucs=True).normalize_structure(s)
    io = PDB.PDBIO()
    io.set_structure(out_s)
    io.save(os.path.join(SAVE_DIR, sf))
