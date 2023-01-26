import sys, os
import csv
import logging
import traceback
from collections import Counter
import numpy as np

from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from rna_bits.utils.data_path import get_path

from utils import *

STRUCT_DIR = get_path("interim/loops/norm_representative/")
MCOUT_DIR = get_path("interim/loops/mca_out/")
OUT_DIR = get_path("interim/loops/ss_annotation/", create=True)

weird_unknown = []
full_unknown = []


# TODO: Use rna_bits.utils versions of stuff that's here

# Takes an input that looks like A149-A151.A or A-44-A-33 or A12-B12-C2324
# Returns ["A149", "A151.A"]
# I assume the chain id is a single alpha character
# I do very little error checking. The purpose of this function is to handle cases with "-" signs
def split_nuc_name_pair(s):
    current = ""
    everything = []
    for c in s:
        if current == "":
            assert c.isalpha()
            current += c
        elif len(current) == 1:
            assert c == "-" or c.isnumeric()
            current += c
        else:
            if c == "-":
                everything.append(current)
                current = ""
            else:
                current += c
    everything.append(current)
    return everything


class MCOut:
    def __init__(self, s):
        global fn
        self.resconf_lines = []
        self.adjstack_lines = []
        self.nonadjstack_lines = []
        self.basepair_lines = []

        current_dest = None
        for line in s.split("\n"):
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("Residue conformations -"):
                current_dest = self.resconf_lines
            elif line.startswith("Adjacent stackings -"):
                current_dest = self.adjstack_lines
            elif line.startswith("Non-Adjacent stackings -"):
                current_dest = self.nonadjstack_lines
            elif line.startswith("Base-pairs -"):
                current_dest = self.basepair_lines
            elif line.startswith("Number of "):
                pass
            else:
                current_dest.append(line)

        self.proper_nucs = set()
        for line in self.resconf_lines:
            a, b = line.split(":")
            nuc_name = a.strip()
            b = b.strip()
            if (
                (
                    b.startswith("G")
                    or b.startswith("C")
                    or b.startswith("A")
                    or b.startswith("U")
                )
                and "unknown" in b
                and "unknown unknown" not in b
            ):
                weird_unknown.append((fn, line))
            if (
                b.startswith("G")
                or b.startswith("C")
                or b.startswith("A")
                or b.startswith("U")
            ) and "unknown unknown" in b:
                full_unknown.append((fn, line))

            if is_valid_nucleotide_conformation(b):
                self.proper_nucs.add(nuc_name)

        self.proper_pairings = {a: None for a in self.proper_nucs}
        for line in self.basepair_lines:
            a, b = line.split(":")
            nuc_a, nuc_b = split_nuc_name_pair(a.strip())
            if (nuc_a not in self.proper_nucs) or (nuc_b not in self.proper_nucs):
                continue
            pairing = b.strip()
            if is_valid_pairing(pairing):
                is_double = False
                if self.proper_pairings[nuc_a] is not None:
                    print(
                        nuc_a,
                        "is paired with both",
                        self.proper_pairings[nuc_a],
                        "and",
                        nuc_b,
                    )
                    is_double = True

                if self.proper_pairings[nuc_b] is not None:
                    print(
                        nuc_b,
                        "is paired with both",
                        self.proper_pairings[nuc_b],
                        "and",
                        nuc_a,
                    )
                    is_double = True

                # if not is_double:
                self.proper_pairings[nuc_b] = nuc_a
                self.proper_pairings[nuc_a] = nuc_b


possible_pairs = ["G-C", "C-G", "A-U", "U-A", "G-U", "U-G"]


def is_valid_pairing(s):
    # Note: this is actually really sketchy. It might be a good idea to print
    # out a copy of every single basepair MC-Annotate reports and stare at them.
    # Or, read the MC-Annotate paper
    a = s.split()
    return (
        (a[0] in possible_pairs)
        and (a[1] in ["Ww/Ww", "Ws/Ww", "Ww/Ws"])
        and ("pairing" in a)
        and ("cis" in a)
        and ("adjacent_5p" not in a)
    )


def is_valid_nucleotide_conformation(s):
    return (
        s.startswith("G ")
        or s.startswith("C ")
        or s.startswith("A ")
        or s.startswith("U ")
    ) and "unknown" not in s


def get_mc_style_name(residue):
    (
        _struct_id,
        _model_id,
        chain_id,
        (_hetero, res_id, insertion_code),
    ) = residue.get_full_id()
    if insertion_code == " ":
        return chain_id + str(res_id)
    else:
        return chain_id + str(res_id) + "." + insertion_code


def mc_name_to_tuple(s):
    chain_id = s[:1]
    s = s[1:]
    if "." in s:
        a, b = s.split(".")
        res_id = int(a)
        insertion_code = b
    else:
        res_id = int(s)
        insertion_code = " "

    return (chain_id, res_id, insertion_code)


def tuple_to_mc_name(t):
    (chain_id, res_id, insertion_code) = t
    if insertion_code in (" ", ""):
        return chain_id + str(res_id)
    else:
        return chain_id + str(res_id) + "." + insertion_code


mcout_filenames = [a for a in os.listdir(MCOUT_DIR) if a.endswith(".txt")]
mcout_filenames.sort()
# struct_filenames = struct_filenames[:10]
assert mcout_filenames

seen_pairings = {p: Counter() for p in possible_pairs}

for i, fn in enumerate(mcout_filenames):
    print(fn, str(i + 1) + "/" + str(len(mcout_filenames)))
    with open(os.path.join(MCOUT_DIR, fn)) as f:
        mc_out = MCOut(f.read())

    struct_fn = remove_string_end(fn, ".txt")
    m = PDBParser().get_structure(struct_fn, os.path.join(STRUCT_DIR, struct_fn))[0]

    all_chains = []
    current_chain = []

    prev = None
    for residue in m.get_residues():
        name = get_mc_style_name(residue)
        if name in mc_out.proper_nucs:

            continuing = False
            if (prev is not None) and ("O3'" in prev) and ("P" in residue):
                dist = np.linalg.norm(prev["O3'"].coord - residue["P"].coord)
                if dist <= 2.0:
                    continuing = True

            if continuing:
                current_chain.append(name)
            else:
                if len(current_chain) > 0:
                    all_chains.append(current_chain)
                    current_chain = []

                current_chain.append(name)

            prev = residue
        else:
            if len(current_chain) > 0:
                all_chains.append(current_chain)
                current_chain = []

            prev = None

    if len(current_chain) > 0:
        all_chains.append(current_chain)
        current_chain = []

    # Sanity checks: each nucleotide should be in at most 1 chain
    cnt = Counter()
    for c in all_chains:
        cnt.update(c)
    for a, b in cnt.items():
        assert b == 1
    assert set(cnt) == mc_out.proper_nucs

    out_fn = os.path.join(OUT_DIR, struct_fn + ".txt")
    with open(out_fn, "w") as f:
        f.write("chains:\n")
        for c in all_chains:
            f.write(" ".join(c))
            f.write("\n")
        f.write("pairs:\n")
        pairings = sorted(
            set(
                tuple(sorted((mc_name_to_tuple(a), mc_name_to_tuple(b))))
                for a, b in mc_out.proper_pairings.items()
                if b is not None
            )
        )
        for a, b in pairings:
            a = tuple_to_mc_name(a)
            b = tuple_to_mc_name(b)
            f.write(a)
            f.write(" ")
            f.write(b)
            f.write("\n")
