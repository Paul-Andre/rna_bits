from weakref import WeakValueDictionary
from typing import List, Sequence, Tuple, TypeVar, Callable
import warnings
import os
import sys
import pickle

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

mmcif_parser = PDB.MMCIFParser()

# "Opportunistic" cache for PDBs that were downloaded
# Objects stay in the cache as long as they haven't been garbage collected.
_pdb_cache = WeakValueDictionary()
_pdb_seen = set()  # For debugging; can probably remove


_PDB_DIR = os.path.join(DATA_PATH, "original/PDB/")


class HiddenPrints:
    """
    Used to block prints (useful when calling some Biopython functions)
    https://stackoverflow.com/a/45669280
    """

    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, "w")
        sys.stderr = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr


def fetch_pdb(
    pdb_code: str,
    cache: MutableMapping[str, Structure] = _pdb_cache,
    pickle_dir=_PDB_DIR,
) -> Structure:
    """
    Dowloads the pdb code.
    Allows both an in-memory cache and an on-disk cache in form of pickle files.
    By default, the returned object might be cached, so copy() it before
    modifying it.
    To avoid caching set cache=None
    """
    pdb_code = pdb_code.lower()
    # Using get() to make sure the object won't be garbage collected between
    # the check and the retrieval (if using WeakValueDictionary)
    if cache is not None:
        res = cache.get(pdb_code)
        if res is not None:
            return res

    pickle_path = os.path.join(pickle_dir, pdb_code + ".pickle")

    try:
        with open(pickle_path, "rb") as f:
            struct = pickle.load(f)
    except (EOFError, OSError) as e:
        if not isinstance(e, FileNotFoundError):
            warnings.warn(f"Failed to read {pickle_path}, redownloading.")

        with HiddenPrints():  # (Silence prints)
            with warnings.catch_warnings():  # (Warning silencing context)
                warnings.simplefilter("ignore", category=PDBConstructionWarning)
                struct = PDB.mmtf.MMTFParser.get_structure_from_url(pdb_code)

        os.makedirs(pickle_dir, exist_ok=True)
        with open(pickle_path, "wb") as f:
            pickle.dump(struct, f)

    if cache is not None:
        cache[pdb_code] = struct
    return struct


def build_model_from_lists_of_residues(a: Sequence[Sequence[Residue]]) -> Model:
    """
    Returns a multi-chain Biopython Model.
    """
    chain_inc = "A"
    model = Model(0)
    for b in a:
        chain = Chain(chain_inc)
        chain_inc = chr(ord(chain_inc) + 1)  # "Increment letter"
        res_inc = 1
        for residue in b:
            new_residue = residue.copy()
            # TODO: sus, I don't know if I'm allowed to assign to .id
            new_residue.id = (" ", res_inc, " ")
            res_inc += 1
            chain.add(new_residue)
        model.add(chain)
    return model


class UnsupportedError(Exception):
    pass


def get_residue_from_chain(chain: Chain, unit: UnitId) -> Residue:
    """Note: there's no validation that the UnitId actually references this
    specific chain
    """
    assert unit.atom_name is None
    # TODO: handle disordered atoms (Will need to copy them over?)
    # Currently I just pretend the altloc selected doesn't exist and I end up
    # selecting all the atoms in the residue.
    # assert x.atom_altloc is None
    if unit.symmetry is not None:
        raise UnsupportedError("Symmetry not supported")

    residue = chain[unit.get_biopython_residue_id()]
    if unit.residue_name is not None:
        # "disordered residues" means that there's a "point mutation", meaning
        # there's two residues that have the same number but a different
        # residue name, hence we need to discern them using residue_name
        if isinstance(residue, PDB.Residue.DisorderedResidue):
            raise UnsupportedError(
                "Disordered residues not supported "
                "(Yet. They shouldn't be too hard to implement)"
            )
        assert isinstance(residue, Residue)
        assert (
            residue.resname == unit.residue_name
        ), "resname mismatch. pointwise mutation?"
    return residue


def get_residue_from_model(model: Model, unit: UnitId) -> Residue:
    """Note: there's no validation that the UnitId actually references this
    specific model
    """
    chain = model[unit.chain_id]
    residue = get_residue_from_chain(chain, unit)
    return residue


def fetch_residue(unit: UnitId) -> Residue:
    structure = fetch_pdb(unit.pdb_code)
    model = structure[unit.model_id - 1]
    residue = get_residue_from_model(model, unit)
    return residue


def fetch_residues(units: Sequence[UnitId], assert_same_model=True) -> List[Residue]:
    """
    Please make a copy of the returned objects before modifying them.
    """
    # TODO (Maybe make a copy myself?) In the case of disordered atoms I will
    # need to do that anyway.
    if assert_same_model:
        pdb_code = units[0].pdb_code
        model_id = units[0].model_id
        for x in units:
            assert x.pdb_code == pdb_code
            assert x.model_id == model_id

    ret = []
    for unit in units:
        ret.append(fetch_residue(unit))

    return ret


# TODO: make path accept a file handle
def save_model_as_pdb(model: Model, path: str) -> None:
    structure = Structure("saving_as_pdb")
    structure.add(model)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(path)


def save_model_as_cif(
    model: Model,
) -> None:
    assert False, "Not Implemented"


T = TypeVar("T")


def query_segmented(
    pattern: str,
    ids: Sequence[T],
    querying_function: Callable[[T], Residue],
    copy=False,
) -> List[List[Residue]]:
    """pattern is something like '(...()..().....)'
    ids are given in whatever format querying_function understands
    """
    assert len(pattern) == len(ids)
    chains = [[]]
    for i, (id, c) in enumerate(zip(ids, pattern)):
        if c == ")" and i != len(pattern) - 1:
            chains.append([])
        residue = querying_function(id)
        if copy:
            residue = residue.copy()
        chains[-1].append(residue)

    return chains
