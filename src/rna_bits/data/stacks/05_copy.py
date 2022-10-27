import subprocess
import os
from collections import Counter
import sys
import warnings
import shutil

from Bio.PDB.PDBList import PDBList

from rna_bits.utils.data_path import get_path
import rna_bits.utils.pdb
from rna_bits import utils

MODEL_DIR = get_path("interim/stacks/extracted/2_2/")
DATABASE_DIR = get_path("database/rna_bits/canonical", create=True)
DEST_DIR = os.path.join(DATABASE_DIR, "2_2")

if os.path.exists(DEST_DIR):
    shutil.rmtree(DEST_DIR)
shutil.copytree(MODEL_DIR, DEST_DIR)
