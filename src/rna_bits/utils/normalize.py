from weakref import WeakValueDictionary
from typing import (
    Union,
    List,
    Sequence,
    Tuple,
    TypeVar,
    Callable,
    TextIO,
    Dict,
    Optional,
)
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

from . import resources


def load_rename_list(filename):
    with open(filename) as f:
        return parse_rename_list(f)


def parse_rename_list(s: Union[str, TextIO]) -> Dict[str, Optional[str]]:
    """Takes the contents of file of the form:

    # heavy atoms
    C2  C2
    C4  C4
    C5  C5
    # ...
    # backbone
    C1' C1'
    C1* C1'
    # ...
    P1 -

    and return a dict. If the atom is renamed to "-" it is set to None in
    the dict."""
    if isinstance(s, str):
        s = io.StringIO(s)

    ret = {}
    for l in s:
        l = l.strip()
        if len(l) == 0 or l.startswith("#"):
            continue
        a, b = l.split()
        if b == "-":
            ret[a] = None
        else:
            ret[a] = b

    return ret


with importlib_resources.open_text(resources, "residues.list") as f:
    RESIDUE_RENAMES: Dict[str, str] = parse_rename_list(f)

with importlib_resources.open_text(resources, "atoms.list") as f:
    ATOM_RENAMES: Dict[str, str] = parse_rename_list(f)


def choose_rename(
    name: str, renames: Dict[str, Optional[str]], keep_unknown=False
) -> Optional[str]:
    if name in renames:
        return renames[name]
    if keep_unknown:
        return name
    return None


def is_nuc(
    residue: Residue, residue_renames: Dict[str, Optional[str]] = RESIDUE_RENAMES
) -> bool:
    if residue_renames.get(residue.resname, "-") in "AUCG":
        return True
    else:
        return False


class Normalizer:
    def __init__(
        self,
        residue_renames: Dict[str, Optional[str]] = RESIDUE_RENAMES,
        atom_renames: Dict[str, Optional[str]] = ATOM_RENAMES,
        keep_unknown_residues=False,
        keep_unknown_atoms=False,
        delete_hydrogens=True,
    ):
        self.residue_renames = residue_renames
        self.atom_renames = atom_renames
        self.keep_unknown_residues = keep_unknown_residues
        self.keep_unknown_atoms = keep_unknown_atoms
        self.delete_hydrogens = delete_hydrogens


    # TODO: do this
    # def normalize_nucleotide(self, residue: Residue) -> 

    def normalize_chain(self, chain: Chain) -> Chain:
        """Creates a new, normalized, chain."""
        out_chain = PDB.Chain.Chain(chain.id)
        for residue in chain:
            new_res_name = choose_rename(
                residue.resname, self.residue_renames, self.keep_unknown_residues
            )
            if new_res_name is None:
                continue

            out_residue = PDB.Residue.Residue(
                residue.id, new_res_name, residue.segid
            )
            out_chain.add(out_residue)

            for atom in residue:
                if self.delete_hydrogens:
                    if atom.element == "H":
                        continue

                if is_nuc(residue, self.residue_renames):
                    new_atom_name = choose_rename(
                        atom.name, self.atom_renames, self.keep_unknown_atoms
                    )
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

        return out_chain


    def normalize_model(self, model: Model) -> Model:
        """Creates a new, normalized, model."""
        out_model = PDB.Model.Model(model.id)

        for chain in model:
            out_chain = self.normalize_chain(chain)
            out_model.add(out_chain)

        return out_model

    def normalize_structure(
        self,
        structure: Structure,
    ) -> Structure:
        """Creates a new, normalized, structure."""
        out_struct = Structure(structure.id)

        for model in structure:
            out_model = self.normalize_model(model)
            out_struct.add(out_model)

        return out_struct
