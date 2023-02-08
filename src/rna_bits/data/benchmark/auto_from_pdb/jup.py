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

from rna_bits.data.benchmark.auto_from_pdb.collect import collect, collect_multi


info = collect("benchmark/auto_from_pdb/all/info/")

singles = {}
singles["from_native"] = collect_multi("benchmark/auto_from_pdb/all/insert_loops_force_native/measures")
singles["loops_top"] = collect_multi("benchmark/auto_from_pdb/all/insert_loops_top_1_exclude_native/measures")

multis =  {}
multis["loops_sample"] = collect_multi("benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/measures")
multis["bp2_reliable"] = collect_multi("benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_reliable/measures")
multis["bp2_all"] = collect_multi("benchmark/auto_from_pdb/bp2_limited/bp2_with_ss_all/measures")


def draw_table(mesure, best_is_max:bool):
    df = info[.copy()
    for (k,v) in singles.items():




