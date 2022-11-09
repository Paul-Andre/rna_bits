import sys
import os
import numpy as np
import random
import io

import rna_bits.insert_loops.match as match
from rna_bits.utils.ss import parse_parens

def output_match(out_f, loop_seq, loop_ss, loop_nuc_ids, match):
    filename = match[1]
    out_f.write(
        " ".join(
            map(
                str,
                (
                    "# wanted:",
                    loop_seq,
                    loop_ss,
                    " selected:",
                    "".join(match[2]["seqs"]),
                    match[2]["full_ss"],
                    " score:",
                    match[0],
                ),
            )
        )
        + "\n"
    )
    out_f.write(
        "motif:" + filename + ":  " + ",".join(str(n + 1) for n in loop_nuc_ids)
    )
    out_f.write("\n")
    out_f.write("\n")


def insert_loops(seq, dot_bracket, out_fs, rest_of_file, exclude=None):
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

    for out_f in out_fs:
        out_f.write(seq)
        out_f.write("\n")
        out_f.write("".join(dot_bracket))
        out_f.write("\n")
        out_f.write(rest_of_file)
        if len(rest_of_file) >= 1 and rest_of_file[-1] != "\n":
            out_f.write("\n")
        out_f.write("\n")
        out_f.write("# Excluding " + repr(exclude))
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
            if exclude is not None and len(exclude) > 0:
                matches = [
                    a for a in matches if all(e not in a[2]["file"] for e in exclude)
                ]
            if len(matches) > 0:
                how_much = 10
                top_ten = matches[:how_much]
                for out_f in out_fs:
                    i, chosen_match = random.choice(list(enumerate(top_ten)))
                    if how_much > 1:
                        out_f.write(
                            f"# top {i+1}, selected uniformly at random from the top {len(top_ten)} matches\n"
                        )
                    output_match(out_f, loop_seq, loop_ss, nucs, chosen_match)

            elif loop_ss != "(())":
                print("Couldn't find a motif for ", loop_seq, loop_ss)


def in_and_out(in_f, out_fs, exclude=None):
    seq = in_f.readline().strip()
    dot_bracket = in_f.readline().strip()
    rest = in_f.read()
    for l in rest.split("\n"):
        if l.startswith("native:") and exclude is None:
            exclude = [s.strip() for s in l[len("native:") :].split(",")]
    insert_loops(seq, dot_bracket, out_fs, rest, exclude=exclude)


def main_cli():
    if len(sys.argv) >= 2:
        f = open(sys.argv[1])
    else:
        f = sys.stdin

    out_dir = sys.argv[2]
    num_times = int(sys.argv[3])

    out_fs = []
    for i_ in range(num_times):
        # Use StringIO so that we are not limited by the OS regarding how many
        # "files" we can write to at the same time
        out_fs.append(io.StringIO())

    in_and_out(f, out_fs, exclude=None)

    # Write the strings to files
    for (i_, string_io) in enumerate(out_fs):
        i = i_ + 1
        print(f"Writing to {i}")
        with open(os.path.join(out_dir, str(i) + ".rass"), "w") as out_f:
            string_io.seek(0)
            text = string_io.read()
            out_f.write(text)

if __name__ == "__main__":
    main_cli()
