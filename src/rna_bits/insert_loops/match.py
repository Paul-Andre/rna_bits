# Matches a loop
import json
import os


from rna_bits.utils.data_path import get_path

# TODO: these can fail depending if the libraries have been generated or not so
# putting these in the global scope is sketchy.
LOOP_LIBRARY_DIR = get_path("out/loops/")
JSON_DIR = get_path("out/loops/json/")


def get_simple_ss(ss):
    return "".join((a if a in "()" else ".") for a in ss)


def letter_score(a, b):
    if a == b:
        return 0
    elif (a in "CU" and b in "CU") or (a in "AG" and b in "AG"):
        return -1
    else:
        return -2


# t=target, l=library
# l=letter as in ACGU c=character as in ".()"
def letter_and_seq_score(tl, tc, ll, lc):
    if tc in "()":
        assert lc in "()"
        return letter_score(tl, ll) * 0.25

    if (tc != ".") and (lc != "."):
        return letter_score(tl, ll) * 0.25

    if tc == "." and lc != ".":
        return letter_score(tl, ll) * 0.5 - 1
    elif tc != "." and lc == ".":
        return letter_score(tl, ll) * 0.5 - 2

    return letter_score(tl, ll)


def calculate_score(tar_seq, tar_ss, lib_seq, lib_ss):
    # the ss looks like (...()...]]]..)
    # anything character that isn't in "()." is assumed to represent a pseudoknotted
    # region
    assert len(tar_seq) == len(tar_ss)
    assert len(tar_ss) == len(lib_seq)
    assert len(lib_seq) == len(lib_ss)
    assert get_simple_ss(tar_ss) == get_simple_ss(lib_ss)
    score = 0
    for tl, tc, ll, lc in zip(tar_seq, tar_ss, lib_seq, lib_ss):
        score += letter_and_seq_score(tl, tc, ll, lc)

    return score


def match(seq, ss):
    assert len(seq) == len(ss)

    comp_sizes = []
    cur_size = 0

    for i, c in enumerate(ss):
        if i != len(ss) - 1 and c == ")":
            comp_sizes.append(cur_size)
            cur_size = 0
        cur_size += 1
    comp_sizes.append(cur_size)

    assert sum(comp_sizes) == len(seq)

    fn = "_".join(map(str, comp_sizes)) + ".json"
    if not os.path.isfile(os.path.join(JSON_DIR, fn)):
        return []

    with open(os.path.join(JSON_DIR, fn)) as f:
        info = json.load(f)
    scored = []
    for entry in info:
        lib_seq = "".join(entry["seqs"])
        lib_ss = entry["full_ss"]
        score = calculate_score(seq, ss, lib_seq, lib_ss)
        scored.append((score, entry))

    scored.sort(key=lambda a: (-a[0], a[1]["file"]))

    return scored


if __name__ == "__main__":
    seq = input()
    ss = input()
    assert len(seq) == len(ss)

    split_seq = []
    split_ss = []
    cur_seq = ""
    cur_ss = ""
    for i, (l, c) in enumerate(zip(seq, ss)):
        if i != len(seq) - 1 and c == ")":
            split_seq.append(cur_seq)
            split_ss.append(cur_ss)
            cur_seq = ""
            cur_ss = ""
        cur_seq += l
        cur_ss += c
    split_seq.append(cur_seq)
    split_ss.append(cur_ss)

    comp_sizes = tuple(len(a) for a in split_seq)
    assert tuple(len(a) for a in split_ss) == comp_sizes

    scored = match(seq, ss)

    for a in reversed(scored):
        print(a)
