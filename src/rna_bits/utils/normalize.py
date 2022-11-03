from weakref import WeakValueDictionary
from typing import Union, List, Sequence, Tuple, TypeVar, Callable, TextIO, Dict
import warnings
import os
import sys
import pickle
import importlib.resources as importlib_resources

from Bio import PDB
import Bio.PDB.mmtf
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

from Bio.PDB.PDBExceptions import PDBConstructionWarning

from rna_bits.utils.bgsu import UnitId
from rna_bits.utils.data_path import DATA_PATH


def parse_rename_list(s: Union[str, TextIO]):
    if isinstance(s, str):
        s = io.StringIO(s)

    return dict(
        l.strip().split()
        for l in s
        if not l.startswith("#") and len(l.strip()) != 0
    )


from . import resources

with importlib_resources.open_text(resources, "residues.list") as f:
    RES_RENAMES: Dict[str, str] = parse_rename_list(f)

with importlib_resources.open_text(resources, "atoms.list") as f:
    ATOM_RENAMES: Dict[str, str] = parse_rename_list(f)

DELETE_RESIDUES = True
DELETE_ATOMS = True


def canonicalize_res_name(name, res_dict):
    if name in res_dict and res_dict[name] != "-":
        return res_dict[name]
    else:
        if DELETE_RESIDUES:
            return None
        else:
            return name


def canonicalize_atom_name(name, atom_dict):
    if name in atom_dict and atom_dict[name] != "-":
        return atom_dict[name]
    else:
        if DELETE_ATOMS:
            return None
        else:
            return name


# Creates a new model that has canonical atom representation
def canonicalize_structure(struct, res_renames=RES_RENAMES, atom_renames=ATOM_RENAMES):
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

                new_res_name = canonicalize_res_name(residue.resname, res_renames)
                if new_res_name is None:
                    continue
                out_residue = PDB.Residue.Residue(
                    residue.id, new_res_name, residue.segid
                )
                out_chain.add(out_residue)

                for atom in residue:
                    if atom.element == "H":
                        continue

                    if is_nuc(residue, res_renames):
                        new_atom_name = canonicalize_atom_name(atom.name, atom_renames)
                        if new_atom_name is None:
                            continue
                    else:
                        new_atom_name = atom.name

                    full_atom_name = new_atom_name
                    out_atom = PDB.Atom.Atom(
                        name=new_atom_name,
                        coord=atom.coord,
                        bfactor=atom.bfactor,
                        occupancy=atom.occupancy,
                        altloc=atom.altloc,
                        fullname=full_atom_name,
                        serial_number=atom.serial_number,
                        element=atom.element,
                    )
                    out_residue.add(out_atom)

    return out_struct


def normalize(struct):
    return canonicalize_structure(struct)


def is_nuc(residue, res_dict):
    if res_dict.get(residue.resname, "-") in "AUCG":
        return True
    else:
        return False
