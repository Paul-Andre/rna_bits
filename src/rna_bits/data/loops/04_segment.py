import sys, os
import csv
import logging
import traceback
from collections import Counter

import numpy as np

from rna_bits.utils.data_path import get_path
from rna_bits.utils.ss import parse_ss_file
from rna_bits.utils.ss import Segmenter

SS_DIR = get_path("interim/loops/ss_annotation/")
OUT_DIR = get_path("interim/loops/segmentations/", create=True)

SEGMENT_EXTERNAL_LOOPS = False

filenames = [a for a in os.listdir(SS_DIR) if a.endswith(".txt")]
filenames.sort()
# filenames = filenames[:10]
assert filenames

for i, fn in enumerate(filenames):
    print(fn, str(i + 1) + "/" + str(len(filenames)))
    with open(os.path.join(SS_DIR, fn)) as f:
        chains, pairs = parse_ss_file(f)

    segmenter = Segmenter(chains, pairs, remove_lonely_pairs=False)

    loops = segmenter.segment_loops()
    if SEGMENT_EXTERNAL_LOOPS:
        loops += segmenter.segment_external_loops()

    out_fn = fn

    if len(loops) == 0:
        continue

    with open(os.path.join(OUT_DIR, out_fn), "w") as f:
        for ss, nucs in loops:
            # print(" ".join(nucs))
            # print(ss)

            f.write(ss)
            f.write(" ")
            f.write(" ".join(nucs))
            f.write("\n")
    print(len(loops))
