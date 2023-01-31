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

from rna_bits.utils.data_path import get_path
from rna_bits.utils.ss import parse_ss_file
from rna_bits.utils.ss import Segmenter
from rna_bits.utils.misc import remove_string_end
from rna_bits.utils.mc_annotate import query_mc_name, mc_name_to_tuple


SS_DIR = get_path("interim/loops/ss_annotation/")
STRUCT_DIR = get_path("interim/loops/norm_representative/")


def generate_structure_files(OUT_DIR, allow_multi_chain=True, allow_junctions=True):
    filenames = [a for a in os.listdir(SS_DIR) if a.endswith(".txt")]
    filenames.sort()
    #filenames = filenames[:10]

    good_files = []
    no_insertion_code=0
    for i,fn in enumerate(filenames):
        #if fn != "NR_4.0_94683.1.pdb.txt":
        #    pass
        #    #continue
        # if fn != "NR_4.0_03110.1.pdb.txt":
        #     continue


        print(fn, str(i+1)+"/"+str(len(filenames)))
        with open(os.path.join(SS_DIR, fn)) as f:
            chains, pairs = parse_ss_file(f)

        segmenter = Segmenter(chains, pairs)

        # only one strand
        if not allow_multi_chain and len(segmenter.fp_nucs) != 1:
            continue

        
        internal_loops = segmenter.segment_loops()
        external_loops = segmenter.segment_external_loops()

        if not allow_multi_chain and len(segmenter.chains) != 1: continue
        assert(not allow_multi_chain), "It's not clear if allow_multi_chain works properly"
        # TODO: for multi-chain
        if not (1 <= len(segmenter.chains) <= 2): continue 

        good = True
        for xl_ss, _xl_nucs in external_loops:
            if set(xl_ss).union(set("().")) != set("()."):
                # Loop contains a crossover
                good = False
                break;
            # no complicated external loop
            if xl_ss.count("()") != 1:
                good = False
                break;
            # Only 5' or 3' loose strand allowed, not both:
            # TODO: perhaps a better strategy would be to just clip them off?
            # Or, better yet, allow rna_insert_loops to insert "external loops" and
            # see how well it works.
            if xl_ss[0] == "." and xl_ss[-1] == ".":
                good = False
                break;

        if not good:
            continue

        print(segmenter.fp_nucs)
        # If there's not weird external loops, there should be 1 external loop per 5' nucleotide
        assert len(external_loops) == len(segmenter.fp_nucs)



        pairing = segmenter.pairing


        for loop_ss, _loop_nucs in internal_loops:
            if set(loop_ss).union(set("().")) != set("()."):
                # Loop contains a crossover
                good = False
                break;
            if not allow_junctions and loop_ss.count("()") > 1:
                # A 3-way junction, for example, will look like (...()..().....)
                good = False
                break;
        if not good:
            continue

        # Order the chains so that the dot-bracket works out
        # TODO ^
        # As a hack, if there's <=2 chains, they don't need to be ordered
        if not (1 <= len(segmenter.chains) <= 2):
            print(segmenter.chains, "chains")
            print("Not processing")
            continue;

        print(len(segmenter.chains))
        if (len(segmenter.chains) >= 2): exit()

        Space = object()

        nucs_in_order = []
        for chain in segmenter.chains:
            nucs_in_order.extend(chain)
            nucs_in_order.append(Space)
        nucs_in_order.pop()

        rev_chain = {c:i for i,c in enumerate(nucs_in_order) if c is not Space}

        dot_bracket = ""
        stack = []  # Used for sanity check
        for i,nuc in enumerate(nucs_in_order):
            if nuc is Space:
                dot_bracket+="+";
                continue
            if nuc not in pairing:
                dot_bracket+= "."
                continue
            j = rev_chain[pairing[nuc]]
            if (j>i):
                dot_bracket+= "("
                stack.append(i)
                continue
            else:
                assert(j<i)
                jj = stack.pop()
                assert(j==jj)
                dot_bracket+= ")"
                continue

        sequence = ""

        # Open the pdb to get the sequence
        struct_fn = remove_string_end(fn, ".txt")
        m = PDBParser().get_structure(struct_fn, os.path.join(STRUCT_DIR, struct_fn))[0]
        for nuc in nucs_in_order:
            if nuc is Space:

                sequence+="+"
                continue
            residue = query_mc_name(m, nuc)
            resname = residue.resname
            assert(resname in "CGAU")
            sequence+=resname

        rass_fn = remove_string_end(fn, ".pdb.txt")+".rass"
        with open(os.path.join(OUT_DIR, rass_fn), "w") as f:
            # TODO: was there an utility to create rass files?
            f.write(sequence)
            f.write("\n")
            f.write(dot_bracket)
            f.write("\n")
            f.write("native: ")
            f.write(remove_string_end(fn, ".pdb.txt"))
            f.write("\n")


        # .index files are used by RNA-puzzles' RNA_assessment to indicate which nucleotides are to be compared
        index_fn = remove_string_end(fn, ".pdb.txt")+".cons.index"
        ref_index_fn = remove_string_end(fn, ".pdb.txt")+".ref.index"
        with open(os.path.join(OUT_DIR, index_fn), "w") as f:
            f.write("A:1:"+str(len(sequence)))

        # TODO: Adjust for multiple chains.
        # Note the the "chain" variable is a loop variable from an above loop,
        # but this just happens to work in the case where there's a single
        # chain.
        assert(not allow_multi_chain), "This hasn't been adjusted for multiple chains"
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


        good_files.append((len(dot_bracket), dot_bracket, fn))
        # TODO: assert that if I parse this dot-brackets string it gives the initial pairing
        # For that I will need to import something from somewhre.

        print(dot_bracket)

    good_files.sort()
    for a in good_files:
        print(*a)

    print(len(good_files))
    print("no_insertion_code:", no_insertion_code)



OUT_DIR_ALL = get_path("benchmark/auto_from_pdb/all/provided_ss", create=True)
OUT_DIR_PB2 = get_path("benchmark/auto_from_pdb/bp2_limited/provided_ss", create=True)

generate_structure_files(OUT_DIR_ALL, allow_multi_chain=False, allow_junctions=True)
generate_structure_files(OUT_DIR_PB2, allow_multi_chain=False, allow_junctions=False)
