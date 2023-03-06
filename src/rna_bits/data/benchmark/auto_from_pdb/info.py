""" info such as num_nucs and original pdb_chains for the test cases
"""
import sys,os
import csv
import logging
import traceback
from collections import Counter
from os.path import join as pjoin
import json

import numpy as np
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from rna_bits.utils.data_path import get_path
from rna_bits.utils.ss import parse_ss_file
from rna_bits.utils.ss import Segmenter, parse_parens
from rna_bits.utils.misc import remove_string_end
from rna_bits.utils.mc_annotate import query_mc_name, mc_name_to_tuple


SS_DIR = get_path("interim/loops/ss_annotation/")
STRUCT_DIR = get_path("interim/loops/norm_representative/")

LIST_FILE = pjoin(get_path("provided"), "nrlist_3.267_4.0A.csv")

FORCE_NATIVE = get_path("benchmark/auto_from_pdb/all/insert_loops_force_native")



def ss_stats(ss):
    num_junctions = 0

    current_rep = []

    full_pairing = [None] * len(ss)

    for a, b in parse_parens(ss, start=0):
        full_pairing[a] = b
        full_pairing[b] = a

    main_pairing = [None] * len(ss)

    def has_helix_fragment(x, pairing=main_pairing):
        if x >= len(ss) - 1:
            return False
        xx = x + 1

        y = pairing[x]
        if y is None:
            return False

        yy = pairing[xx]
        if yy is None:
            return False

        if yy >= len(ss) - 1:
            return False

        return yy + 1 == y

    # Remove lonely pairs
    ss_str = ss
    ss = list(ss)
    for (k, v) in enumerate(full_pairing):
        if v is None:
            continue
        if not has_helix_fragment(k, full_pairing) and not has_helix_fragment(
            v, full_pairing
        ):
            full_pairing[k] = None
            full_pairing[v] = None
            ss[k] = "."
            ss[v] = "."

    for k, v in enumerate(full_pairing):
        if ss[k] in "()":
            main_pairing[k] = v

    loops = []
    for x in range(len(ss)):
        if main_pairing[x] and main_pairing[x] > x:
            y = main_pairing[x]
            a = x
            loop_ss = ""
            nucs = []
            loop_ss += ss[a]
            nucs.append(a)
            a += 1
            loop_ss += ss[a]
            nucs.append(a)

            while a != y:
                if main_pairing[a] and main_pairing[a] > a:
                    a = main_pairing[a]
                else:
                    a = a + 1
                loop_ss += ss[a]
                nucs.append(a)
            if loop_ss != "(())":
                typ = loop_ss.count("()") + 1
                if typ >= 3:
                    num_junctions+=1
                    current_rep.append(typ)

    out = {}
    out["num_nucs"] = len(ss_str.strip("."))
    out["num_nucs_with_dangling"] = len(ss)
    out["num_paired"] = ss.count("(")*2
    out["num_junctions"] = num_junctions
    if len(current_rep):
        out["junction_profile"] = "+".join(map(str,sorted(current_rep)))
    else:
        out["junction_profile"] = "No"
    return out


def generate_info(DIR):
    # TODO: hardcoded LIST_FILE

    with open(LIST_FILE, newline="") as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)

    pdb_by_nrlist = {}
    all_by_nrlist = {}
    for (name, rep_s, all_s) in rows:
        pdb_by_nrlist[name] = rep_s
        all_by_nrlist[name] = all_s

    PROVIDED_DIR = get_path(pjoin(DIR, "provided_ss"))
    OUT_DIR = get_path(pjoin(DIR, "info"), create=True)

    filenames = [a for a in os.listdir(PROVIDED_DIR) if a.endswith(".ref.index_noloose")]
    filenames.sort()
    assert(filenames)

    for i, fn in enumerate(filenames):
        print(fn, str(i+1) + "/" + str(len(filenames)))

        nat = remove_string_end(fn, ".ref.index_noloose")

        with open(pjoin(PROVIDED_DIR, nat+".cons.index_noloose")) as f:
            _,_,num = f.read().strip().split(":")
            num = int(num)

        with open(pjoin(PROVIDED_DIR, nat+".rass")) as f:
            seq = f.readline()
            ss = f.readline().strip()
            print(ss)

        out = {}
        out["pdb_chains"] = pdb_by_nrlist[nat]
        out["equivalent_nrlist_entries"] = all_by_nrlist[nat]
        out["num_nucs"] = num

        out_ss = ss_stats(ss)
        assert(out_ss["num_nucs"] == out["num_nucs"])

        out.update(out_ss)

        with open(os.path.join(OUT_DIR, nat+".json"), "w") as f:
            f.write(json.dumps(out))


if __name__ == "__main__":
    DIR = "benchmark/auto_from_pdb/all/"
    generate_info(DIR);
    DIR = "benchmark/auto_from_pdb/bp2_limited/"
    generate_info(DIR);

