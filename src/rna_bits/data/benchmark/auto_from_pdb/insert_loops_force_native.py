import sys
import os
import numpy as np
from os.path import join as pjoin

from rna_bits.utils.data_path import get_path
from rna_bits.utils.misc import remove_string_end
from rna_bits.insert_loops import match
from rna_bits.utils.ss import parse_parens



def parse_index(f):
    cons_units = []
    for c, s, l in (a.strip().split(":") for a in f.read().strip().split(",")):
        s = int(s)
        l = int(l)
        for i in range(s, s + l):
            cons_units.append(c + str(i))
            # cons_units.append((i))
    return cons_units


def insert_loops(seq, dot_bracket, out_f, rest_of_file, cons_to_ref, native=None):
    assert len(seq) == len(dot_bracket)

    original_seq = seq
    seq = seq.upper()

    full_pairing = [None] * len(seq)

    for a, b in parse_parens(dot_bracket, start=0):
        full_pairing[a] = b
        full_pairing[b] = a

    main_pairing = [None] * len(seq)

    def has_helix_fragment(x, pairing=main_pairing):
        if x >= len(seq) - 1:
            return False
        xx = x + 1

        y = pairing[x]
        if y is None:
            return False

        yy = pairing[xx]
        if yy is None:
            return False

        if yy >= len(seq) - 1:
            return False

        return yy + 1 == y

    # Remove lonely pairs
    dot_bracket = list(dot_bracket)
    for (k, v) in enumerate(full_pairing):
        if v is None:
            continue
        if not has_helix_fragment(k, full_pairing) and not has_helix_fragment(
            v, full_pairing
        ):
            full_pairing[k] = None
            full_pairing[v] = None
            dot_bracket[k] = "."
            dot_bracket[v] = "."

    for k, v in enumerate(full_pairing):
        if dot_bracket[k] in "()":
            main_pairing[k] = v

    out_f.write(seq)
    out_f.write("\n")
    out_f.write("".join(dot_bracket))
    out_f.write("\n")
    out_f.write(rest_of_file)
    if len(rest_of_file) >= 1 and rest_of_file[-1] != "\n":
        out_f.write("\n")
    out_f.write("\n")
    out_f.write("# Forcing " + repr(native))
    out_f.write("\n")
    out_f.write("\n")
    loops = []
    for x in range(len(seq)):
        if main_pairing[x] and main_pairing[x] > x:
            y = main_pairing[x]
            a = x
            loop_seq = ""
            loop_ss = ""
            nucs = []
            loop_seq += seq[a]
            loop_ss += dot_bracket[a]
            nucs.append(a)
            a += 1
            loop_seq += seq[a]
            loop_ss += dot_bracket[a]
            nucs.append(a)

            while a != y:
                if main_pairing[a] and main_pairing[a] > a:
                    a = main_pairing[a]
                else:
                    a = a + 1
                loop_seq += seq[a]
                loop_ss += dot_bracket[a]
                nucs.append(a)
            matches = match.match(loop_seq, loop_ss)
            matches = [
                (score, os.path.join(match.LOOP_LIBRARY_DIR, d["file"]), d)
                for (score, d) in matches
            ]
            want = cons_to_ref["A" + str(nucs[0] + 1)]
            if native is not None and len(native) > 0:
                matches = [
                    a
                    for a in matches
                    if (
                        any(e in a[2]["file"] for e in native)
                        and a[0] == 0.0
                        and a[2]["original_nucs"][0] == want
                    )
                ]
            if loop_ss != "(())":
                if len(matches) != 1:
                    print(matches)
                assert len(matches) == 1
                # print("#", matches[:10])
                # print(matches[0])
                best_match = matches[0]
                assert best_match[2]["original_nucs"][0] == want
                best_file = best_match[1]
                out_f.write(
                    " ".join(
                        map(
                            str,
                            (
                                "# wanted:",
                                loop_seq,
                                loop_ss,
                                " found:",
                                "".join(best_match[2]["seqs"]),
                                best_match[2]["full_ss"],
                                " score:",
                                best_match[0],
                            ),
                        )
                    )
                    + "\n"
                )
                out_f.write(
                    "motif:" + best_file + ":  " + ",".join(str(n + 1) for n in nucs)
                )
                out_f.write("\n")
                out_f.write("\n")


def in_and_out(in_f, out_f, cons_to_ref, native=None):
    seq = in_f.readline().strip()
    dot_bracket = in_f.readline().strip()
    rest = in_f.read()
    for l in rest.split("\n"):
        if l.startswith("native:") and native is None:
            native = [s.strip() for s in l[len("native:") :].split(",")]
        
    insert_loops(seq, dot_bracket, out_f, rest, cons_to_ref=cons_to_ref, native=native)


def run_insert_loops_force_native(RASS_DIR, OUT_DIR):
    struct_filenames = [a for a in os.listdir(RASS_DIR) if a.endswith(".rass")]
    struct_filenames.sort()
    assert(len(struct_filenames))

    for i, fn in enumerate(struct_filenames):
        print(fn, str(i + 1) + "/" + str(len(struct_filenames)))
        ffn = os.path.join(RASS_DIR, fn)
        nat = remove_string_end(fn, ".rass")

        offn = pjoin(OUT_DIR, nat, "1.rass")
        os.makedirs(pjoin(OUT_DIR, nat), exist_ok=True)

        try:
            with open(os.path.join(RASS_DIR, nat + ".cons.index_noloose")) as f:
                cons_units = parse_index(f)
            with open(os.path.join(RASS_DIR, nat + ".ref.index_noloose")) as f:
                ref_units = parse_index(f)
        except IOError:
            continue

        cons_to_ref = dict(zip(cons_units, ref_units))

        with open(ffn) as inf:
            with open(offn, "w") as outf:
                in_and_out(inf, outf, cons_to_ref)


if __name__ == "__main__":
    DIR = "benchmark/auto_from_pdb/all/"
    RASS_DIR = get_path(pjoin(DIR, "provided_ss"))
    OUT_DIR = get_path(pjoin(DIR, "insert_loops_force_native/rass"), create=True)
    run_insert_loops_force_native(RASS_DIR, OUT_DIR)
