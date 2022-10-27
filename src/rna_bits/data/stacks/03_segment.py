import subprocess
import os
from collections import Counter
import sys
import warnings

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBList import PDBList

from rna_bits.utils import mc_annotate
from rna_bits.utils import ss
from rna_bits.utils.data_path import get_path
from rna_bits.utils.misc import remove_string_end

STRUCT_DIR = get_path("interim/stacks/pdb")
# MCOUT_DIR = get_path("interim/stacks/mcout")
SS_DIR = get_path("interim/stacks/ss")
SEG_DIR = get_path("interim/stacks/segmentations", create=True)

filenames = [a for a in os.listdir(SS_DIR) if a.endswith(".ss")]
filenames.sort()
assert filenames

for i, fn in enumerate(filenames):
    print(fn, str(i + 1) + "/" + str(len(filenames)))
    with open(os.path.join(SS_DIR, fn)) as f:
        chains, pairs = ss.parse_ss_file(f)

    segmenter = ss.Segmenter(chains, pairs)

    base_code = remove_string_end(fn, ".ss")

    stacks = []
    for helix in segmenter.helices:
        st_a = helix.nucs[0]
        st_b = list(reversed(helix.nucs[1]))
        assert len(st_a) == len(st_b)
        for i in range(1, len(st_a)):
            a = st_a[i - 1]
            b = st_a[i]
            c = st_b[i]
            d = st_b[i - 1]
            stacks.append((a, b, c, d))
            stacks.append((c, d, a, b))

    seg_fn = os.path.join(SEG_DIR, base_code + ".segmentation")
    with open(seg_fn, "w", encoding="ascii") as f:
        for a, b, c, d in stacks:
            f.write(f"(()) {a} {b} {c} {d}\n")

    # if not stacks:
    #     continue

    # m = PDBParser().get_structure(struct_fn, os.path.join(STRUCT_DIR, struct_fn))[0]
