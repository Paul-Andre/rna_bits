
import sys,os
import csv
from Bio import PDB
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import numpy as np


def split_nuc_name_pair(s):
    """
    Takes an input that looks like A149-A151.A or A-44-A-33 or A12-B12-C2324
    or '0'123-'0'124
    Returns ["A149", "A151.A"], ["A-44", "A-33"], ["A12", "B12", "C2324"],
    ["'0'123, '0'124]
    I assume the chain id is a single alpha character
    The reason this function exists instead
    of just using "split" is to handle the cases where residue number are
    negative.
    I do very little error checking for cases that aren't correctly formatted.
    """
    current = ""
    everything = []
    i = 0
    while i < len(s):
        c = s[i]
        if current == "":
            if c.isalpha():
                current+=c
            elif c == "'":
                current+=c
                i+=1
                c=s[i]
                current+=c
                i+=1
                c=s[i]
                current+=c
        elif len(current) == 1:
            assert c == "-" or c.isnumeric()
            current+=c
        else:
            if c == "-":
                everything.append(current)
                current = ""
            else:
                current+=c

        i+=1

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
            elif (line.startswith("Number of ")):
                pass
            else:
                current_dest.append(line)


        self.proper_nucs = set()
        for line in self.resconf_lines:
            a,b = line.split(":")
            nuc_name = a.strip()
            b = b.strip()
            if (
                b.startswith("G") or
                b.startswith("C") or
                b.startswith("A") or
                b.startswith("U")
            ) and "unknown" in b and "unknown unknown" not in b:
                weird_unknown.append((fn, line))
            if (
                b.startswith("G") or
                b.startswith("C") or
                b.startswith("A") or
                b.startswith("U")
            ) and "unknown unknown" in b:
                full_unknown.append((fn, line))

            if is_valid_nucleotide_conformation(b):
                self.proper_nucs.add(nuc_name)


        self.proper_pairings = {a:None for a in self.proper_nucs}
        for line in self.basepair_lines:
            a,b = line.split(":")
            nuc_a, nuc_b = split_nuc_name_pair(a.strip())
            if (nuc_a not in self.proper_nucs) or (nuc_b not in self.proper_nucs):
                continue
            pairing = b.strip()
            if is_valid_pairing(pairing):
                is_double = False
                if self.proper_pairings[nuc_a] is not None:
                    print(nuc_a, "is paired with both", self.proper_pairings[nuc_a], "and", nuc_b)
                    is_double = True

                if self.proper_pairings[nuc_b] is not None:
                    print(nuc_b, "is paired with both", self.proper_pairings[nuc_b], "and", nuc_a)
                    is_double = True

                #if not is_double:
                self.proper_pairings[nuc_b] = nuc_a
                self.proper_pairings[nuc_a] = nuc_b



possible_pairs = ["G-C", "C-G", "A-U", "U-A", "G-U", "U-G"]

def is_valid_pairing(s):
    # Note: this is actually really sketchy. It might be a good idea to print
    # out a copy of every single basepair MC-Annotate reports and stare at them.
    # Or, read the MC-Annotate paper
    a = s.split()
    return (
        (a[0] in possible_pairs) and
        (a[1] in ["Ww/Ww", "Ws/Ww", "Ww/Ws"]) and
        ("pairing" in a) and
        ("cis" in a) and
        ("adjacent_5p" not in a)
    )

def is_valid_nucleotide_conformation(s):
    return (
            s.startswith("G ") or
            s.startswith("C ") or
            s.startswith("A ") or
            s.startswith("U ")
        ) and "unknown" not in s


def mc_name_to_tuple(s):
    if s[0] != "'":
        chain_id = s[:1]
        s = s[1:]
    else:
        pos1 = s[1:].find("'") + 1
        chain_id = s[1:pos1]
        s = s[pos1+1:]
    if "." in s:
        a,b = s.split(".")
        res_id = int(a)
        insertion_code = b
    else:
        res_id = int(s)
        insertion_code = " "

    return (chain_id, res_id, insertion_code)


def tuple_to_mc_name(t):
    (chain_id, res_id, insertion_code) = t
    if any(c.isdigit() for c in chain_id):
        chain_id = "'" + chain_id + "'"
    if (insertion_code in (" ", "")):
        return chain_id+str(res_id)
    else:
        return chain_id+str(res_id)+"."+insertion_code

def get_mc_style_name(residue):
    (_struct_id, _model_id, chain_id, (_hetero, res_id, insertion_code)) = residue.get_full_id()
    return tuple_to_mc_name((chain_id, res_id, insertion_code))

def query_mc_name(model, name):
    (chain_id, res_id, insertion_code) = mc_name_to_tuple(name)
    return model[chain_id][(" ",res_id, insertion_code)]
