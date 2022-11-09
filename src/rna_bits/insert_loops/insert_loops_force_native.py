import sys
import os
import numpy as np
import match

sr_dir = "/home/paul/LoopLibrary/simple_rass/"


def parse_index(f):
    cons_units = []
    for c, s, l in (a.strip().split(":") for a in f.read().strip().split(",")):
        s = int(s)
        l = int(l)
        for i in range(s, s + l):
            cons_units.append(c + str(i))
            # cons_units.append((i))
    return cons_units


def insert_loops(seq, dot_bracket, out_f, rest_of_file, native=None):
    assert len(seq) == len(dot_bracket)
    nat = native[0]
    with open(os.path.join(sr_dir, nat + ".cons.index_noloose")) as f:
        cons_units = parse_index(f)
    with open(os.path.join(sr_dir, nat + ".ref.index_noloose")) as f:
        ref_units = parse_index(f)

    cons_to_ref = dict(zip(cons_units, ref_units))

    original_seq = seq
    seq = seq.upper()

    full_pairing = [None] * len(seq)

    import Project.parser

    for a, b in Project.parser.parseParens(dot_bracket):
        a -= 1
        b -= 1
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
            print(want)
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


def in_and_out(in_f, out_f, native=None):
    seq = in_f.readline().strip()
    dot_bracket = in_f.readline().strip()
    rest = in_f.read()
    for l in rest.split("\n"):
        if l.startswith("exclude:") and native is None:
            exclude = [s.strip() for s in l[len("exclude:") :].split(",")]
    insert_loops(seq, dot_bracket, out_f, rest, native=native)


# if __name__ == "__main__":
if False:
    if len(sys.argv) >= 2:
        f = open(sys.argv[1])
    else:
        f = sys.stdin

    if len(sys.argv) >= 3 and sys.argv[2] == "-n":
        in_and_out(f, sys.stdout, exclude=[])
    else:
        in_and_out(f, sys.stdout, exclude=None)
else:
    OUT_DIR = "/home/paul/LoopLibrary/simple_rass_motifs_force_native"
    if not os.path.isdir(OUT_DIR):
        os.mkdir(OUT_DIR)

    for fn in (
        fn
        for fn in os.listdir("/home/paul/LoopLibrary/simple_rass")
        if fn.endswith("rass")
    ):
        print(fn[: -len(".rass")])
        ffn = os.path.join("/home/paul/LoopLibrary/simple_rass", fn)
        offn = os.path.join(
            "/home/paul/LoopLibrary/simple_rass_motifs_force_native", fn
        )
        with open(ffn) as inf:
            with open(offn, "w") as outf:
                # in_and_out(inf, outf, native = [fn[:-len(".rass")]])
                nat = fn[: -len(".rass")]
                try:
                    f = open(os.path.join(sr_dir, nat + ".ref.index_noloose"))
                    f.close()
                except IOError:
                    continue

                in_and_out(inf, outf, native=[fn[: -len(".rass")]])
