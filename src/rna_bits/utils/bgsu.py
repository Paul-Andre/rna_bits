import dataclasses
from dataclasses import dataclass
import urllib
import json
from typing import Optional, Tuple, Sequence, List
import os
import warnings
import time

from Bio import PDB

from .data_path import DATA_PATH

MOTIF_GROUP_PATH = os.path.join(DATA_PATH, "original/rna3dhub/")


@dataclass(frozen=True, eq=True)
class UnitId:
    """
    https://www.bgsu.edu/research/rna/help/rna-3d-hub-help/unit-ids.html
    """

    pdb_code: Optional[str] = None
    model_id: Optional[
        int
    ] = None  # Note: starts at 1, whereas Biopython models start at 0
    chain_id: Optional[str] = None
    residue_name: Optional[str] = None
    residue_id: Optional[int] = None
    atom_name: Optional[str] = None
    atom_altloc: Optional[str] = None
    insertion_code: Optional[str] = None
    symmetry: Optional[str] = None

    # TODO: maybe just make it part of __init__
    @classmethod
    def parse(cls, s: str) -> "UnitId":
        """
        Parse a BGSU unit id as described here:
        https://www.bgsu.edu/research/rna/help/rna-3d-hub-help/unit-ids.html
        Example ids:
        "1ABC|1|A", "1ABC|1|B|U|10", "6TQA|1|B|ARG|188||A", "1J5E|1|A|G|190|||G"

        """
        a = s.split("|")
        if a[-1] == "":
            a.pop()
        assert len(a) <= 9, f"Can't understand BGSU unit id '{id}'"
        while len(a) < 9:
            a.append("")

        def str_or_none(t: str):
            if t != "":
                return t
            return None

        def int_or_none(t: str):
            if t != "":
                return int(t)
            return None

        return UnitId(
            pdb_code=str_or_none(a[0]),
            model_id=int_or_none(a[1]),
            chain_id=str_or_none(a[2]),
            residue_name=str_or_none(a[3]),
            residue_id=int_or_none(a[4]),
            atom_name=str_or_none(a[5]),
            atom_altloc=str_or_none(a[6]),
            insertion_code=str_or_none(a[7]),
            symmetry=str_or_none(a[8]),
        )

    def __str__(self) -> str:
        """Print the id in the BGSU format"""
        # TODO test
        def str_or_empty(a):
            if a is None:
                return ""
            return str(a)

        a = [str_or_empty(x) for x in dataclasses.astuple(self)]
        while a[-1] == "":
            a.pop()
        return "|".join(a)

    def __repr__(self) -> str:
        return f"UnitId.parse('{str(self)}')"

    def get_biopython_residue_id(self) -> Tuple[str, int, str]:
        assert self.residue_id is not None

        if self.insertion_code is not None:
            ic = self.insertion_code
        else:
            ic = " "
        # Note: the first element is used to indicate hetero-residues
        # I assume that BGSU residues are never hetero but I didn't check
        return (" ", self.residue_id, ic)

    def contains(self, other) -> bool:
        if self.pdb_code is None:
            return True  # Uhh.. okay...
        if other.pdb_code != self.pdb_code:
            return False

        if self.model_id is None:
            return True
        if other.model_id != self.model_id:
            return False

        if self.chain_id is None:
            return True
        if other.chain_id != self.chain_id:
            return False

        if self.residue_id is None:
            return True
        if other.residue_id != self.residue_id or other.insertion_code != self.insertion_code:
            return False

        if self.atom_name is None:
            return True
        if other.atom_name != self.atom_name:
            return False

        if self.atom_altloc is None:
            return True
        if other.atom_altloc != self.atom_altloc:
            return False

        # TODO: make this work in situation when atom_altloc is set, but
        # atom_name isn't (representing a mutated version of a nucleotide)



def download_motif_group_info(group_id: str):
    url = f"http://rna.bgsu.edu/rna3dhub/motif/view/{group_id}/json"
    path = MOTIF_GROUP_PATH + "/" + group_id + ".json"
    try:
        with open(path, "rt") as f:
            data = json.load(f)
    except FileNotFoundError:

        # Download while retrying with exponential backoff
        retries = 5
        wait_secs = 1
        while True:
            try:
                # Try downloading
                with urllib.request.urlopen(url) as request:
                    s = request.read()
                    break
            except urllib.error.URLError:
                if retries > 0:
                    print("retrying")
                    warnings.warn(f"Failed to get {url} , retrying in {wait_secs}")
                    time.sleep(wait_secs)
                    wait_secs *= 2
                    retries -= 1
                else:
                    raise

        data = json.loads(s)
        os.makedirs(MOTIF_GROUP_PATH, exist_ok=True)
        with open(path, "wb") as f:
            f.write(s)

    return data


# def get_unit_in_model(model: PDB.Model.Model, id: UnitId):
#     chain = model[id.chain_id]
#     residue =
