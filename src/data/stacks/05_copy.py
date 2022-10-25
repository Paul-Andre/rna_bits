import utils.pdb
from Bio.PDB.PDBList import PDBList
from utils.data_path import get_path
import subprocess
import os
from collections import Counter
import sys
import warnings

import shutil

MODEL_DIR = get_path("interim/stacks/extracted/2_2/")
DATABASE_DIR = get_path("database/rna_bits/canonical", create=True)
DEST_DIR = os.path.join(DATABASE_DIR, "2_2")

if os.path.exists(DEST_DIR):
    shutil.rmtree(DEST_DIR)
shutil.copytree(MODEL_DIR, DEST_DIR)

