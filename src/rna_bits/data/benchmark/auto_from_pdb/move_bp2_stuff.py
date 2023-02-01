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
from rna_bits.utils.misc import remove_string_end

PROVIDED = get_path("benchmark/auto_from_pdb/bp2_limited/provided_ss")
OUTS = get_path("benchmark/auto_from_pdb/all/bp2_with_ss_reliable/bp2_output")
OUTS_NEW = get_path("benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_reliable/bp2_output")

struct_filenames = [a for a in os.listdir(PROVIDED) if a.endswith(".rass")]
struct_filenames.sort()

for i,fn in enumerate(struct_filenames):
    nat = remove_string_end(fn, ".rass")
    ii = pjoin(OUTS, nat+".json")
    oo = pjoin(OUTS_NEW, nat+".json")
    shutil.move(ii, oo)

