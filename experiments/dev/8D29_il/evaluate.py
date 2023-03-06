import sys, os
import csv
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import subprocess
import os
import shutil
import json

from rna_bits.utils.data_path import get_path
from rna_bits.utils.misc import remove_string_end




import sys, os

from os.path import abspath

# https://github.com/RNA-Puzzles/RNA_assessment
import RNA_normalizer

from operator import attrgetter


def get_mid_file(start, end):
    return (a[: -len(end)] for a in os.listdir(start) if a.endswith(end))


def file_pattern_intersection(*all_patterns):
    sets_of_mids = [set(get_mid_file(a, b)) for (a, b) in all_patterns]
    union_mids = set.union(*sets_of_mids)
    intersection_mids = set.intersection(*sets_of_mids)

    # # Print what files don't coincide
    # intersection_mids = set.intersection(*sets_of_mids)
    # for mid in union_mids - intersection_mids:
    #     has = []
    #     hasnt = []
    #     for (a,b),mids in zip(all_patterns,sets_of_mids):
    #         if mid in mids:
    #             has.append(a+mid+b)
    #         else:
    #             hasnt.append(a+mid+b)
    #     print("The file(s)", has, "exist, but the file(s)", hasnt, "do not.")

    print("We have", len(intersection_mids))

    return [
        (mid,) + tuple(a + mid + b for (a, b) in all_patterns)
        for mid in intersection_mids
    ]


def InteractionNetworkFidelity(
    prediction_file, prediction_index, native_file, native_index
):
    res_struct = RNA_normalizer.PDBStruct()
    res_struct.load(native_file, native_index)
    res_raw_seq = res_struct.raw_sequence()

    sol_struct = RNA_normalizer.PDBStruct()
    sol_struct.load(prediction_file, prediction_index)
    sol_raw_seq = sol_struct.raw_sequence()

    if sol_raw_seq != res_raw_seq:
        sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
        sys.stderr.write("DATA Solution sequence --> '%s'\n" % sol_raw_seq)
        sys.stderr.write("DATA Result sequence   --> '%s'\n" % res_raw_seq)
        raise Exception("sequence don't match")

    # computes the RMSD
    comparer = RNA_normalizer.PDBComparer()
    rmsd = comparer.rmsd(sol_struct, res_struct)
    INF_ALL = comparer.INF(sol_struct, res_struct, type="ALL")
    DI_ALL = rmsd / INF_ALL
    pvalue = comparer.pvalue(rmsd, len(sol_raw_seq), "-")

    INF_WC = comparer.INF(sol_struct, res_struct, type="PAIR_2D")
    INF_NWC = comparer.INF(sol_struct, res_struct, type="PAIR_3D")
    INF_STACK = comparer.INF(sol_struct, res_struct, type="STACK")
    return {
        "rmsd": rmsd,
        "pvalue": pvalue,
        "DI_ALL": DI_ALL,
        "INF_ALL": INF_ALL,
        "INF_WC": INF_WC,
        "INF_NWC": INF_NWC,
        "INF_STACK": INF_STACK,
    }


def run_evaluate():

    REF_PDB_DIR = get_path("interim/loops/norm_representative/")
    PDB_DIR = "./"

    PROVIDED_DIR = get_path("benchmark/auto_from_pdb/all/provided_ss/")


    OUT_DIR = get_path("experiments/8D29_insert_loops/measures/", create=True)

    pdb_dirs = [a for a in os.listdir(PDB_DIR)]
    print(pdb_dirs)
    pdb_dirs = [a for a in pdb_dirs if os.path.dirname(os.path.join(PDB_DIR,a)) and a.startswith("loops_sample")]
    print(PDB_DIR)
    print(pdb_dirs)
    assert(len(pdb_dirs))
    pdb_dirs.sort()

    for i, fn in enumerate(pdb_dirs):
        struct_filenames = [a for a in os.listdir(os.path.join(PDB_DIR,fn)) if a.endswith(".rass.pdb")]
        struct_filenames.sort()

        for ii, fnn in enumerate(struct_filenames):
            os.makedirs(os.path.join(OUT_DIR, fn), exist_ok=True)


            print(fnn, fn, str(ii + 1) + "/" + str(len(struct_filenames)), str(i+1) + "/" + str(len(pdb_dirs)))

            clipped = remove_string_end(fnn, ".rass.pdb")

            cons = os.path.join(PDB_DIR, fn, fnn)
            nr = "NR_4.0_56838.1"
            cons_index = os.path.join(PROVIDED_DIR, nr +".cons.index_noloose")

            ref = os.path.join(REF_PDB_DIR, nr + ".pdb")
            ref_index = os.path.join(PROVIDED_DIR, nr +".ref.index_noloose")
            print(cons_index)
            print(ref_index)

            if not os.path.isfile(cons_index):
                # If this is the case, it probably means the file was excluded because of insertion codes
                continue
            if not os.path.isfile(ref_index):
                # If this is the case, it probably means the file was excluded because of insertion codes
                continue

            out = InteractionNetworkFidelity(cons, cons_index, ref, ref_index)
            print(out)

            with open(os.path.join(OUT_DIR, fn, clipped+".json"), "w") as f:
                f.write(json.dumps(out))


# Calculate 


# DIR = "benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/"
# run_builder(DIR)

# def run_evaluate(DIR):
#     results = []

#     cons_pattern = (
#         "/home/paul/LoopLibrary/simple_rass_motifs_possibly_native/",
#         ".rass.pdb",
#     )
#     # cons_pattern = ("/home/paul/LoopLibrary/simple_rass_motifs_force_native/", ".rass.pdb")
#     # cons_pattern = ("/home/paul/LoopLibrary/simple_rass_motifs/", ".rass.pdb")
#     # cons_pattern = ("/home/paul/LoopLibrary/simple_rass_motifs/", ".rass_minimize.pdb")

#     ref_pattern = ("/home/paul/LoopLibrary/norm_representative/", ".pdb")
#     # ref_pattern = ("/home/paul/LoopLibrary/simple_rass_motifs/", ".rass.pdb")

#     # cons_index_pattern = ("/home/paul/LoopLibrary/simple_rass/", ".cons.index")
#     cons_index_pattern = ("/home/paul/LoopLibrary/simple_rass/", ".cons.index_noloose")

#     # ref_index_pattern = ("/home/paul/LoopLibrary/simple_rass/", ".ref.index")
#     ref_index_pattern = ("/home/paul/LoopLibrary/simple_rass/", ".ref.index_noloose")
#     # ref_index_pattern = ("/home/paul/LoopLibrary/simple_rass/", ".cons.index_noloose")

#     all_patterns = (cons_pattern, cons_index_pattern, ref_pattern, ref_index_pattern)

#     print(*all_patterns)

#     cnt = 0
#     rmsd_tot = 0.0

#     failed = []
#     for name, cons, cons_index, ref, ref_index in file_pattern_intersection(
#         *all_patterns
#     ):
#         print(name)
#         try:
#             out = InteractionNetworkFidelity(cons, cons_index, ref, ref_index)
#             cnt += 1
#             rmsd_tot += out["rmsd"]
#         except Exception:
#             out = {
#                 "rmsd": "",
#                 "pvalue": "",
#                 "DI_ALL": "",
#                 "INF_ALL": "",
#                 "INF_WC": "",
#                 "INF_NWC": "",
#                 "INF_STACK": "",
#             }
#             failed.append(name)

#         print(out)
#         out.update(name=name)
#         # to get the number of nucleotides, look at the cons_index file
#         with open(cons_index) as f:
#             l = int(f.readline().strip().split(":")[2])
#         out.update(num_nucs=l)
#         results.append(out)

#     columns = [
#         "name",
#         "num_nucs",
#         "rmsd",
#         "DI_ALL",
#         "INF_ALL",
#         "INF_WC",
#         "INF_NWC",
#         "INF_STACK",
#     ]

#     print("\t".join(columns))
#     for a in results:
#         print("\t".join(str(a.get(k, "")) for k in columns))

#     print(*all_patterns)
#     print("cnt", cnt, "rmsd_tot", rmsd_tot, "average", rmsd_tot / cnt)
#     print("failed", failed)

if __name__ == "__main__":
    # DIR = "benchmark/auto_from_pdb/all/insert_loops_top_1_exclude_native/"
    # run_evaluate(DIR)

    # DIR = "benchmark/auto_from_pdb/all/insert_loops_sample_50_exclude_native/"
    # run_evaluate(DIR)

    run_evaluate()



