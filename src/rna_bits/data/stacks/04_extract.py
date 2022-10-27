import subprocess
import os
from collections import Counter
import sys
import warnings

import numpy as np
from Bio.PDB.PDBList import PDBList
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from rna_bits.utils import mc_annotate
from rna_bits.utils import ss
from rna_bits.utils.data_path import get_path
from rna_bits.utils.misc import remove_string_end
from rna_bits.utils.pdb import query_segmented
from rna_bits.utils.pdb import build_model_from_lists_of_residues
from rna_bits.utils.pdb import save_model_as_pdb

STRUCT_DIR = get_path("interim/stacks/pdb")
# MCOUT_DIR = get_path("interim/stacks/mcout")
SS_DIR = get_path("interim/stacks/ss")
SEG_DIR = get_path("interim/stacks/segmentations")
MODEL_DIR = get_path("interim/stacks/extracted/2_2/", create=True)


def get_model_dir(sequence):
    directory = get_path("interim/stacks/extracted/2_2/" + sequence, create=True)
    return directory


letters = "AUGC"

filenames = [a for a in os.listdir(SEG_DIR) if a.endswith(".segmentation")]
filenames.sort()
assert len(filenames)

for i, fn in enumerate(filenames):
    print(fn, str(i + 1) + "/" + str(len(filenames)))
    stacks = []
    with open(os.path.join(SEG_DIR, fn)) as f:
        for l in f:
            stuff = l.strip().split()
            assert stuff[0] == "(())"
            stacks.append(stuff[1:])

    base_code = remove_string_end(fn, ".segmentation")

    struct_fn = base_code + ".pdb"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        m = PDBParser().get_structure(struct_fn, os.path.join(STRUCT_DIR, struct_fn))[0]

    for stack_i, stack in enumerate(stacks):
        chains = query_segmented(
            "(())", stack, lambda id: mc_annotate.query_mc_name(m, id)
        )

        sequence = ""
        for l in chains:
            for r in l:
                resname = r.get_resname()
                assert resname in letters
                sequence += resname

        out_model = build_model_from_lists_of_residues(chains)

        out_name = os.path.join(get_model_dir(sequence), f"{base_code}_{stack_i+1}.pdb")

        save_model_as_pdb(out_model, out_name)
