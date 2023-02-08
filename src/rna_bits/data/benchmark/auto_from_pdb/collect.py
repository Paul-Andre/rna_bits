""" Used to open a bunch of json files and to put the info into dataframes
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
import pandas as pd

from rna_bits.utils.data_path import get_path
from rna_bits.utils.ss import parse_ss_file
from rna_bits.utils.ss import Segmenter
from rna_bits.utils.misc import remove_string_end
from rna_bits.utils.mc_annotate import query_mc_name, mc_name_to_tuple

# TODO: put this stuff in a more general directory in the library. Oh, well

def collect(DIR):
    """
    DIR = "benchmark/auto_from_pdb/all/info/";
    """
    REAL_DIR = get_path(DIR)
    filenames = [a for a in os.listdir(REAL_DIR) if a.endswith(".json")]
    filenames.sort()
    assert(filenames)

    data_list = []

    for i, fn in enumerate(filenames):
        nat = remove_string_end(fn, ".json")

        with open(pjoin(REAL_DIR, fn)) as f:
            datum = json.load(f)

        datum["DIR"] = DIR
        datum["REAL_DIR"] = REAL_DIR
        datum["name"] = nat
        data_list.append(datum)

    df = pd.DataFrame(data_list)
    df.set_index(["name"], inplace=True)
    return df

def collect_multi(DIR):
    """
    DIR = "benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/measures"
    """
    REAL_DIR = get_path(DIR)
    filenames = [a for a in os.listdir(REAL_DIR) if os.path.isdir(pjoin(REAL_DIR, a))]
    filenames.sort()
    assert(filenames)

    data_list = []

    for i, fn in enumerate(filenames):
        instance_fns = [a for a in os.listdir(pjoin(REAL_DIR, fn)) if a.endswith(".json")]
        instance_fns.sort()

        for ii, fnn in enumerate(instance_fns):
            instance_nat = remove_string_end(fnn, ".json")

            with open(pjoin(REAL_DIR, fn, fnn)) as f:
                datum = json.load(f)

            datum["DIR"] = DIR
            datum["REAL_DIR"] = REAL_DIR
            datum["name"] = fn
            datum["instance"] = instance_nat
            data_list.append(datum)

    df = pd.DataFrame(data_list)
    df.set_index(["name", "instance"], inplace=True)
    return df


if __name__ == "__main__":
    # DIR = ("benchmark/auto_from_pdb/all/info/")
    # df = collect(DIR)
    # print(df)

    DIR = ("benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/measures")
    df = collect_multi(DIR)
    print(df)

