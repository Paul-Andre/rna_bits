import sys, os
import csv
import logging
import traceback

from Bio import PDB
import Bio.PDB.mmtf

from rna_bits.utils.data_path import get_path
from rna_bits.utils.pdb import fetch_pdb

# TODO: put this file somewhere else?
# TODO: I think I might have edited the file by hand to fix a model numbering issue, don't do that
LIST_FILE = "nrlist_3.267_4.0A.csv"

SAVE_DIR = get_path("interim/loops/representative/", create=True)

with open(SAVE_DIR + "/version.txt", "w") as f:
    f.write(LIST_FILE[:-4])


def parse_rep_s(a):
    # input ooks like '6ZQD|1|D4+6ZQD|1|D2'
    ret = list((c, int(d), e) for c, d, e in (b.split("|") for b in a.split("+")))
    return ret


with open(LIST_FILE, newline="") as csvfile:
    reader = csv.reader(csvfile)
    rows = list(reader)

rows = rows
log = open("log.txt", "a")
num_succ = 0
for i, (name, rep_s, _) in enumerate(rows):

    # This was to re-process files that weren't properly saved because of too many hetero residues:
    # if name not in ['NR_4.0_06650.22', 'NR_4.0_26150.2', 'NR_4.0_50188.1', 'NR_4.0_52597.1']:
    #    continue

    if name == "NR_4.0_28129.1":
        # This file has too many atoms for the .pdb file format!
       continue

    try:
        # each row in the csv file looks like:
        # ['NR_4.0_07561.2', '6ZQD|1|D4+6ZQD|1|D2', '6ZQD|1|D4+6ZQD|1|D2,7AJU|1|D4+7AJU|1|D2']
        print(name, rep_s, str(i + 1) + "/" + str(len(rows)))
        rep = parse_rep_s(rep_s)
        # save_name = os.path.join(SAVE_DIR,name)
        pdb_id = rep[0][0]
        model_id = rep[0][1]
        assert model_id == 1
        accepted_chain_names = []
        for r in rep:
            assert len(r) == 3
            assert r[0] == pdb_id
            assert r[1] == model_id
            accepted_chain_names.append(r[2])

        print("starting download", pdb_id)
        s = fetch_pdb(pdb_id)
        print("finished download")
        if model_id == 0:
            model_id = 1
            print("for some reason", name, "has model_id", 0)
            print("I assume it was supposed to be 1")
        m = s[model_id - 1]

        out_struct = PDB.Structure.Structure(name)

        out_model = PDB.Model.Model(0)
        out_struct.add(out_model)

        # I rename my chains because .pdb need the chains to have only a single character
        new_chain_id = "A"
        for cn in accepted_chain_names:
            c = m[cn].copy()
            c.id = new_chain_id
            new_chain_id = chr(ord(new_chain_id) + 1)

            # There are some cases where there were more hetero residues than
            # the .pdb file format allows. (MC-Annotate only accepts the .pdb format)
            # I could have deleted the hetero residues, but I decided to keep them.
            # To reduce their number, I delete waters and renumber the rest of the
            # hetero residues.
            #
            # To be honest, at this time I don't actually use the hetero
            # residues anywhere, so they could be as well deleted for simplicity

            to_del = []
            to_renum = []
            max_id = 0
            for r in c:
                if r.resname == "HOH":
                    to_del.append(r.id)
                else:
                    if r.id[0].startswith("H_"):
                        to_renum.append(r)
                    else:
                        max_id = max(max_id, r.id[1])
            for rid in to_del:
                c.detach_child(rid)
            for r in to_renum:
                c.detach_child(r.id)
            for i, r in enumerate(to_renum):
                prev_id = r.id
                new_id = (prev_id[0], max_id + i, prev_id[2])
                r.id = new_id
                c.add(r)
            # end of the code for hetero residue reduction

            out_model.add(c)



        io = PDB.PDBIO()
        io.set_structure(out_struct)
        io.save(os.path.join(SAVE_DIR, name + ".pdb"))

        # log.write("Succeeded "+ name + "\n")
        log.flush()
        num_succ += 1
    except Exception as e:
        m = traceback.format_exc()
        print(m)
        log.write("Failed " + name + ":\n")
        log.write(m)
        log.write("\n")
        log.flush()

print("successful", num_succ, "out of", len(rows))
