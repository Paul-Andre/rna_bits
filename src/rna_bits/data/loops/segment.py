import sys, os
import csv
import logging
import traceback
from collections import Counter
import numpy as np

from Segmenter import Segmenter

# from utils import *

SS_DIR = "ss_annotation/"
OUT_DIR = "segmentations/"

SEGMENT_EXTERNAL_LOOPS = False

filenames = [a for a in os.listdir(SS_DIR) if a.endswith(".txt")]
filenames.sort()
# filenames = filenames[:10]

if not os.path.isdir(OUT_DIR):
    os.mkdir(OUT_DIR)


for i, fn in enumerate(filenames):
    print(fn, str(i + 1) + "/" + str(len(filenames)))
    with open(os.path.join(SS_DIR, fn)) as f:
        segmenter = Segmenter(f)

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
