import sys,os
import csv
from Bio import PDB
import Bio.PDB.mmtf
import logging
import traceback

LIST_FILE = "nrlist_3.226_4.0A.csv"
SAVE_DIR = "representative/"

def parse_rep_s(a):
    # input ooks like '6ZQD|1|D4+6ZQD|1|D2'
    ret = list(
            (c,int(d),e) for
            c,d,e in (b.split("|") for b in a.split("+"))
            )
    return ret

with open(LIST_FILE, newline="") as csvfile:
    reader = csv.reader(csvfile)
    rows = list(reader)
    
rows = rows
log = open("log.txt", "a")
num_succ = 0
for i, (name, rep_s, _) in enumerate(rows):

    # This was to re-process files that weren't properly saved because of too many hetero residues:
    #if name not in ['NR_4.0_06650.22', 'NR_4.0_26150.2', 'NR_4.0_50188.1', 'NR_4.0_52597.1']:
    #    continue

    try:
        # each row in the csv file looks like:
        # ['NR_4.0_07561.2', '6ZQD|1|D4+6ZQD|1|D2', '6ZQD|1|D4+6ZQD|1|D2,7AJU|1|D4+7AJU|1|D2']
        print(name, rep_s, str(i+1)+"/"+str(len(rows)))
        rep = parse_rep_s(rep_s)
        #save_name = os.path.join(SAVE_DIR,name)
        pdb_id = rep[0][0]
        model_id = rep[0][1]
        assert(model_id==1)
        accepted_chain_names = []
        for r in rep:
            assert(len(r) == 3)
            assert(r[0] == pdb_id)
            assert(r[1] == model_id)
            accepted_chain_names.append(r[2])

        print("starting download")
        m = PDB.mmtf.MMTFParser.get_structure_from_url(pdb_id)[model_id-1]
        print("finished download")

        
        out_struct = PDB.Structure.Structure(name)

        out_model = PDB.Model.Model(0)
        out_struct.add(out_model)

        # I rename my chains because .pdb need the chains to have only a single character
        new_chain_id = "A"
        for cn in accepted_chain_names:
            c = m[cn].copy()
            c.id = new_chain_id
            new_chain_id = chr(ord(new_chain_id)+1)
            to_del = []
            # Remove waters (strictly to fit in the .pdb constraints)
            for r in c:
                if r.resname == "HOH":
                    to_del.append(r.id)
            for rid in to_del:
                c.detach_child(rid)

            out_model.add(c)

        io = PDB.PDBIO()
        io.set_structure(out_struct)
        io.save(os.path.join(SAVE_DIR, name+".pdb"))

        #log.write("Succeeded "+ name + "\n")
        log.flush()
        num_succ +=1
    except Exception as e:
        m = traceback.format_exc()
        print(m)
        log.write("Failed "+ name + ":\n")
        log.write(m)
        log.write("\n")
        log.flush()

print("successful", num_succ, "out of", len(rows))

