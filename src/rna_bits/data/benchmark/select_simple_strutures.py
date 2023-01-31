""" Used to generate testing data
"""
import sys,os
import csv
import logging
import traceback
from collections import Counter
import numpy as np
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from Segmenter import Segmenter
from utils import *

SS_DIR = "ss_annotation/"
OUT_DIR = "simple_rass/"
STRUCT_DIR = "norm_representative/"

SEGMENT_EXTERNAL_LOOPS = False

filenames = [a for a in os.listdir(SS_DIR) if a.endswith(".txt")]
filenames.sort()
#filenames = filenames[:10]

if not os.path.isdir(OUT_DIR):
    os.mkdir(OUT_DIR)


good = []
no_insertion_code=0
for i,fn in enumerate(filenames):
    if fn != "NR_4.0_94683.1.pdb.txt":
        pass
        #continue


    print(fn, str(i+1)+"/"+str(len(filenames)))
    with open(os.path.join(SS_DIR, fn)) as f:
        segmenter = Segmenter(f)

    # only one strand
    if len(segmenter.fp_nucs) != 1:
        continue
    
    external_loops = segmenter.segment_external_loops()
    if len(external_loops) != 1:
        continue

    xl_seq = external_loops[0][0]
    # no complicated external loop
    if xl_seq.count("()") != 1:
        continue
    # ony 5' or 3' loose strand allowed, not both
    if xl_seq[0] == "." and xl_seq[-1] == ".":
        continue

    chain = segmenter.chains[0]
    rev_chain = {c:i for i,c in enumerate(chain)}

    pairing = segmenter.pairing

    def has_crossover(l, r):
        if (l > r):
            return False

        while(l <= r):
            while l<=r and chain[l] not in pairing:
                l+=1

            if l>r:
                break;

            a = l
            b = rev_chain[pairing[chain[a]]]
            if (b<a or b not in range(l,r+1)):
                return True

            if has_crossover(a+1, b-1):
                return True

            l = b+1

        return False

    if has_crossover(0, len(chain)-1):
        print("has_crossover")
        continue

    dot_bracket = ""
    for i in range(len(chain)):
        if chain[i] not in pairing:
            dot_bracket+= "."
            continue
        j = rev_chain[pairing[chain[i]]]
        if (j>i):
            dot_bracket+= "("
            continue
        else:
            assert(j<i)
            dot_bracket+= ")"
            continue

    sequence = ""

    # Open the pdb to get the sequence
    struct_fn = remove_string_end(fn, ".txt")
    m = PDBParser().get_structure(struct_fn, os.path.join(STRUCT_DIR, struct_fn))[0]
    for nuc in chain:
        residue = query_mc_name(m, nuc)
        resname = residue.resname
        assert(resname in "CGAU")
        sequence+=resname

    rass_fn = remove_string_end(fn, ".pdb.txt")+".rass"
    with open(os.path.join(OUT_DIR, rass_fn), "w") as f:
        f.write(sequence)
        f.write("\n")
        f.write(dot_bracket)
        f.write("\n")

    # .index files are used by RNA-puzzles' RNA_assessment to indicate which nucleotides are to be compared
    index_fn = remove_string_end(fn, ".pdb.txt")+".cons.index"
    ref_index_fn = remove_string_end(fn, ".pdb.txt")+".ref.index"
    with open(os.path.join(OUT_DIR, index_fn), "w") as f:
        f.write("A:1:"+str(len(sequence)))

    parsed_chain = list(map(mc_name_to_tuple, chain))
    if all(insertion_code== " " for (_chain_id, _res_id, insertion_code) in parsed_chain):
        # RNA_assessment doesn't take into account insertion codes
        # If the structure contains insertion codes, I'll just not output the file
        # TODO: make the normalization step renumber nucleotides
        no_insertion_code += 1
        with open(os.path.join(OUT_DIR, ref_index_fn), "w") as f:
            f.write(",".join(
                chain_id+":"+str(res_id)+":1" for (chain_id, res_id, _insertion_code) in parsed_chain
                ))

    # .index_noloose are index files with loose 5' or 3' strand excluded
    begin = 0
    while(begin<len(dot_bracket) and dot_bracket[begin]=="."):
        begin+=1
    assert(begin < len(dot_bracket))

    end = len(dot_bracket)-1
    while(end>=0 and dot_bracket[end]=="."):
        end-=1
    assert(end >= 0)
        
    # .index files are used by RNA-puzzles' RNA_assessment to indicate which nucleotides are to be compared
    index_nl_fn = remove_string_end(fn, ".pdb.txt")+".cons.index_noloose"
    ref_index_nl_fn = remove_string_end(fn, ".pdb.txt")+".ref.index_noloose"
    with open(os.path.join(OUT_DIR, index_nl_fn), "w") as f:
        f.write("A:"+str(begin+1)+":"+str(end-begin+1))

    parsed_chain = list(map(mc_name_to_tuple, chain))
    if all(insertion_code== " " for (_chain_id, _res_id, insertion_code) in parsed_chain):
        # RNA_assessment doesn't take into account insertion codes
        # If the structure contains insertion codes, I'll just not output the file
        # TODO: make the normalization step renumber nucleotides
        #no_insertion_code += 1
        with open(os.path.join(OUT_DIR, ref_index_nl_fn), "w") as f:
            f.write(",".join(
                chain_id+":"+str(res_id)+":1" for (chain_id, res_id, _insertion_code) in parsed_chain[begin:end+1]
                ))


    good.append((len(dot_bracket), dot_bracket, fn))
    # TODO: assert that if I parse this dot-brackets string it gives the initial pairing
    # For that I will need to import something from somewhre.

    print(dot_bracket)

good.sort()
for a in good:
    print(*a)

print(len(good))
print("no_insertion_code:", no_insertion_code)

