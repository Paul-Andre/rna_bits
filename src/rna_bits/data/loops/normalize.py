import sys,os
import csv
from Bio import PDB
import logging
import traceback
from Bio.PDB.PDBParser import PDBParser
from collections import Counter
import numpy as np

STRUCT_DIR = "representative/"
SAVE_DIR = "norm_representative/"
MCA_DIR = "/home/paul/MC-Annotate"


DELETE_RESIDUES = True
DELETE_ATOMS = True

def get_conv(d):
    with open(d) as f:
        return dict(l.strip().split() for l in f if not l.startswith("#") and len(l.strip())!=0)

res_dict = get_conv("data/residues.list")
atom_dict = get_conv("data/atoms.list")

def canonicalize_res_name(name):
    if name in res_dict and res_dict[name] != "-":
        return res_dict[name]
    else:
        if DELETE_RESIDUES:
            return None
        else:
            return name

def canonicalize_atom_name(name):
    if name in atom_dict and atom_dict[name] != "-":
        return atom_dict[name]
    else:
        if DELETE_ATOMS:
            return None
        else:
            return name

# Creates a new model that has canonical atom representation
def canonicalize_structure(struct):
    out_struct = PDB.Structure.Structure(struct.id)

    for model in struct:
        out_model = PDB.Model.Model(model.id)
        out_struct.add(out_model)

        for chain in model:
            out_chain = PDB.Chain.Chain(chain.id)
            out_model.add(out_chain)

            for residue in chain:
                if residue.resname == "HOH":
                    continue

                new_res_name = canonicalize_res_name(residue.resname)
                if new_res_name is None:
                    continue
                out_residue = PDB.Residue.Residue(residue.id, new_res_name, residue.segid)
                out_chain.add(out_residue)

                for atom in residue:
                    if atom.element == "H":
                        continue

                    if is_nuc(residue):
                        new_atom_name = canonicalize_atom_name(atom.name)
                        if new_atom_name is None:
                            continue
                    else:
                        new_atom_name = atom.name

                    full_atom_name = new_atom_name
                    out_atom = PDB.Atom.Atom(
                            name = new_atom_name,
                            coord = atom.coord,
                            bfactor = atom.bfactor,
                            occupancy = atom.occupancy,
                            altloc = atom.altloc,
                            fullname = full_atom_name,
                            serial_number = atom.serial_number,
                            element = atom.element,
                            )
                    out_residue.add(out_atom)

    return  out_struct

def is_nuc(residue):
    if res_dict.get(residue.resname,"-") in "AUCG":
        return True
    else:
        return False

struct_filenames = [a for a in os.listdir(STRUCT_DIR) if a.endswith(".pdb")]
struct_filenames.sort()
struct_filenames = struct_filenames#[:10]

failed = []

for i,sf in enumerate(struct_filenames):
    print(sf, str(i+1)+"/"+str(len(struct_filenames)))
    s = PDBParser().get_structure(sf, os.path.join(STRUCT_DIR, sf))
    out_s = canonicalize_structure(s)
    io = PDB.PDBIO()
    io.set_structure(out_s)
    io.save(os.path.join(SAVE_DIR, sf))

