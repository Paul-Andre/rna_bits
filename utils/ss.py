import sys,os
import csv
import logging
import traceback
from collections import Counter
import numpy as np
from typing import List, Hashable, Tuple, Union

from Project.parser import parseParens as parse_parens
from Project.parser import isSeparator as _is_separator

class Helix:
    def __init__(self):
        self.nucs = ([],[])
        self.out = [None,None]

class Strand:
    def __init__(self):
        self.nucs = []
        self.out_helix = None
        self.out_helix_half = None

def ts(*a):
    return tuple(sorted(a))

def parens_to_chains(parens, skip_separator=True, start=1):
    """
    ".((&))",skip_separator=True -> [[1,2,3],[4,5]]
    ".((&))",skip_separator=False -> [[1,2,3],[5,6]]
    """
    chains = []
    pos = start
    current = []
    for c in parens:
        if _is_separator(c):
            chains.append(current)
            current = []
            if not skip_separator:
                pos+=1
        else:
            current.append(pos)
            pos+=1
    if len(current):
        chains.append(current)
    return chains


NucId = Hashable
class Segmenter:
    @classmethod
    def from_parens(cls, parens, start=0, remove_lonely_pairs=False):
        chains = parens_to_chains(parens, skip_separator=True, start=start)
        pairs = parse_parens(parens, start=start)
        print(chains)
        print(pairs)
        return Segmenter(chains, pairs, remove_lonely_pairs=remove_lonely_pairs)

    def __init__(self, chains: List[List[NucId]], pairs_original: Tuple[NucId, NucId], remove_lonely_pairs=True):
        nucleotides = []
        fp_nucs = []  # 5' nucs
        next = {}
        pairing = {}
        pairs = []

        for a in chains:
            fp_nucs.append(a[0])
            for i in range(1,len(a)):
                assert(a[i-1] not in next)
                next[a[i-1]] = a[i]
            nucleotides += a
        for (x,y) in pairs_original:
            assert(pairing.get(x) != y)
            assert(pairing.get(y) != x)
            if (x not in pairing) and (y not in pairing):
                # (Don't add double basepairs)
                pairing[x] = y
                pairing[y] = x
                pairs.append((x,y))

        pair_to_helix = {}
        def has_helix_fragment(x):
            xx = next.get(x)
            if xx is None:
                return False

            y = pairing.get(x)
            if y is None:
                return False

            yy = pairing.get(xx)
            if yy is None:
                return False

            return next.get(yy) == y

        # Remove lonely base pairs
        if remove_lonely_pairs:
            new_pairs = []
            for x,y in pairs:
                if not has_helix_fragment(x) and not has_helix_fragment(y):
                    del pairing[x]
                    del pairing[y]
                else:
                    new_pairs.append((x,y))
            pairs = new_pairs

        helices = []

        seen_nucs = Counter() # Strictly for sanity checking
        for (x,y) in pairs:
            if ts(x,y) in pair_to_helix:
                continue

            helix = Helix()

            # First go fully "down"
            while has_helix_fragment(y):
                y = next[y]
                x = pairing[y]

            bot_y = y
            bot_x = x

            # Now go "up"
            pair_to_helix[ts(x,y)] = helix
            while has_helix_fragment(x):
                x = next[x]
                y = pairing[x]
                pair_to_helix[ts(x,y)] = helix

            top_y = y
            top_x = x

            #print(bot_x, top_x, top_y, bot_y)

            half_a = []
            curr = bot_x
            half_a.append(curr)
            while has_helix_fragment(curr):
                curr = next[curr]
                half_a.append(curr)

            half_b = []
            curr = top_y
            half_b.append(curr)
            while has_helix_fragment(curr):
                curr = next[curr]
                half_b.append(curr)

            seen_nucs.update(half_a)
            seen_nucs.update(half_b)
            helix.nucs = (half_a, half_b)

            #print(half_a, half_b)
            assert(len(half_a) == len(half_b))
            #assert(len(half_a) >= 2)

            helices.append(helix)


        ignored_nucs = set()
        strands = []

        def follow_strand(x):
            strand = Strand()

            if x is None:
                strand.out_helix = None
                strand.out_helix_half = None

                strands.append(strand)
                return strand

            while (x not in pairing) and (x in next):
                strand.nucs.append(x)
                x = next[x]

            if x in pairing:
                xx = x
                yy = pairing[xx]
                out_h = pair_to_helix[ts(xx,yy)]

                #print(xx)
                #print(strand.nucs)
                #print(out_h.nucs)
                strand.out_helix = out_h
                if xx == out_h.nucs[0][0]:
                    strand.out_helix_half = 0
                else:
                    assert xx == out_h.nucs[1][0]
                    strand.out_helix_half = 1

                strands.append(strand)
                seen_nucs.update(strand.nucs)
            else:
                # Add loose 3' strand
                strand.nucs.append(x)

                strand.out_helix = None
                strand.out_helix_half = None

                strands.append(strand)
                seen_nucs.update(strand.nucs)

            return strand

        for helix in helices:
            helix.out[0] = follow_strand(next.get(helix.nucs[0][-1]))
            helix.out[1] = follow_strand(next.get(helix.nucs[1][-1]))

        fp_strands = []
        # Add loose 5' strands
        for x in fp_nucs:
            fp_strands.append(follow_strand(x))

        #print(tuple(a.nucs for a in fp_strands))

        #print("seen", sorted(seen_nucs))
        #print("ignore", sorted(ignored_nucs))
        #print(set(seen_nucs))
        #print("################################")
        #print(set(nucleotides))
        #print(set(nucleotides)-set(seen_nucs))
        #print(seen_nucs)
        assert(set(seen_nucs) == set(nucleotides))
        print(seen_nucs)
        print([h.nucs for h in helices])
        print([h.out for h in helices])
        print([h.nucs for h in strands])
        assert(all(v == 1 for v in seen_nucs.values()))

        self.pairing = pairing
        self.helices = helices
        self.strands = strands
        self.fp_strands = fp_strands

        self.nucleotides = nucleotides
        self.fp_nucs = fp_nucs
        self.next = next
        self.pairing = pairing
        self.pairs = pairs
        self.chains = chains


    def traverse_strand(self, strand, target, target_half, loops, visited=set(), stack=[]):

        if strand is None:
            return
        visited.add(strand)
        stack.append(strand.nucs)
        try:
            helix = strand.out_helix
            half = strand.out_helix_half

            if helix == target:
                if helix is None:
                    loops.append(list(stack))
                elif half == target_half:
                    stack.append((helix.nucs[half][0], helix.nucs[1-half][-1]))
                    loops.append(list(stack))
                    stack.pop()
                else:
                    pass
                return

            if helix is None:
                return

            if helix in visited:
                return

            visited.add(helix)

            stack.append(helix.nucs[half])
            self.traverse_strand(helix.out[half], target, target_half, loops)
            stack.pop()

            stack.append((helix.nucs[half][0], helix.nucs[1-half][-1]))
            self.traverse_strand(helix.out[1-half], target, target_half, loops)
            stack.pop()

            visited.remove(helix)

        finally:
            stack.pop()
            visited.remove(strand)

    def display_loop(self, l):
        """
        Output like ("(..().)", [8,9,10,11,2,3,4])
        """
        ss = ""
        nucs = []
        if len(l) % 2 == 1:
            external = True
        else:
            external = False

        if not external:
            assert(self.pairing[l[-1][0]] == l[-1][1])
            ss+="("
            nucs.append(l[-1][1])

        pseudobrackets = "][}{><aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ"
        pbi = 0

        i = 0
        while(i<len(l)):
            s_nucs = l[i]
            ss+="."*len(s_nucs)
            nucs+=s_nucs

            if (i+1<len(l)):
                h_nucs = l[i+1]
                if (len(h_nucs) == 2) and (h_nucs[0] in self.pairing) and (self.pairing[h_nucs[0]] == h_nucs[1]):
                    if i+1 == len(l)-1:
                        ss+=")"
                        nucs.append(h_nucs[0])
                    else:
                        ss+="()"
                        nucs+=h_nucs
                else:
                    ss+=pseudobrackets[pbi]*len(h_nucs)
                    pbi+=1
                    nucs+=h_nucs
            i+=2

        return (ss, nucs)

    def segment_loops(self):
        loops = []
        for helix in self.helices:
            self.traverse_strand(helix.out[0], helix, 1, loops)
            self.traverse_strand(helix.out[1], helix, 0, loops)

        return [self.display_loop(l) for l in loops]

    def segment_external_loops(self):
        loops = []
        for strand in self.fp_strands:
            self.traverse_strand(strand, None, None, loops)

        return [self.display_loop(l) for l in loops]

