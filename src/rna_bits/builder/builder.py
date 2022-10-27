import sys
import numpy as np
from dataclasses import dataclass
from typing import Optional, Union, List, Tuple, Dict, TextIO, Literal, Container, Set
import os

# A shim to make relative imports work when executing this file directly as a
# script (python builder.py). Note, the package *still* needs to have been
# installed.
# This shim is for convenience only and I don't guarantee that it works. If it
# fails, then go read the readme to see what the correct way of executing this
# file is.
if __name__ == "__main__" and __package__ is None:
    __package__ = "builder"
    sys.path.remove(os.path.abspath(os.path.dirname(__file__)))


def get_input_filenames(argv: List[str]) -> List[str]:
    return sys.argv[1:]


@dataclass(frozen=True)
class RnamoipCompInfo:
    name: str
    comp_id: int
    first_nuc: int
    last_nuc: int


@dataclass(frozen=True)
class MotifInfo:
    filename: str
    nucs: List[int]


@dataclass(frozen=True)
class RassData:
    seq: str
    # TODO: RassData has no business modifying seq. It should keep data as-is
    # as much as possible
    original_seq: str
    dot_bracket: str
    # TODO: Maybe make these all a single list to keep their priority
    user_motifs: List[MotifInfo]
    stacks: List[Tuple[int, int]]
    rnamoip_components: List[RnamoipCompInfo]


def parse_and_expand_ranges(ranges: str) -> List[int]:
    """Takes in a string like "1, 2,5-8" and returns [1,2,5,6,7,8]"""
    ret = []
    for r in ranges.split(","):
        r = r.strip()
        if "-" in r:
            first, last = r.split("-")
            first = int(first)
            last = int(last)
        else:
            first = int(r)
            last = first
        ret.extend(range(first, last + 1))
    return ret


def parse_rass_file(f: TextIO) -> RassData:
    seq = f.readline().strip()
    original_seq = seq
    seq = seq.upper()
    dot_bracket = f.readline().strip()
    assert len(seq) == len(dot_bracket)
    user_motifs = []
    stacks = []
    rnamoip_components = []
    while True:
        line = f.readline()
        if not line:
            break
        line = line.strip()
        if line.startswith("%"):
            pass
        elif line.startswith("C-"):
            info = line.split("-")
            _, name, first, last, comp = info
            first, last, comp = map(int, (first, last, comp))
            rnamoip_components.append(
                RnamoipCompInfo(name=name, comp_id=comp, first_nuc=first, last_nuc=last)
            )
        elif line.startswith("motif:"):
            _, name, ranges = line.split(":")
            nucs = parse_and_expand_ranges(ranges)
            nucs = [n for n in nucs]
            user_motifs.append(MotifInfo(filename=name, nucs=nucs))
        elif line.startswith("stack"):
            # The same kind of stack that would form in a helix.
            # Used to simulate coaxial stacking
            # Currently, only doing stacking of nucleotids that are joined by the backbone
            # Stacking makes sense in other contexts, but I'm not doing it right now
            rest = line[len("stack") :].strip()
            a, b = rest.split("-")
            a = int(a)
            b = int(b)
            # TODO: warn instead of asserting
            assert a >= 1
            assert b <= len(seq)
            assert b > a
            for i in range(a, b):
                stacks.append((i, i + 1))

    return RassData(
        seq=seq,
        original_seq=original_seq,
        dot_bracket=dot_bracket,
        user_motifs=user_motifs,
        stacks=stacks,
        rnamoip_components=rnamoip_components,
    )


import rna_bits.utils as utils
import rna_bits.utils.ss


def parse_parens(dot_bracket: str) -> Dict[int, int]:
    pairs = utils.ss.parse_parens(dot_bracket, skip_separator=False, start=1)
    pairing = utils.ss.pairs_to_pairing_dict(pairs)
    return pairing


Priority = Optional[int]


def is_lhs_higher_priority(lhs: Priority, rhs: Priority) -> bool:
    # The smaller the number, the higher the priority. None is the lowest priority
    if lhs is None:
        return False
    if rhs is None:
        return True
    return lhs < rhs


NodeKind = Union[
    Tuple[Literal["module"], str],
    Tuple[Literal["ncm"], str, str],
    Tuple[Literal["pair"], str],
    Tuple[Literal["helix_stack"], str],
    Tuple[Literal["nucleotide"], str],
]

ModelSourceInfo = Union[str, Tuple[str, str]]


class Node:
    id_increment = 10000

    def __init__(self, kind: NodeKind, nucs: List[int]) -> None:
        self.kind = kind  # Stuff like: ("module", "/filnemanaf/asdf/adsf/asdf.pdb"), ("ncm", "2_2", "AGCU")
        self.nucs = nucs
        self.nucs_set = set(nucs)

        self.id = Node.id_increment
        Node.id_increment += 1

        self.priority: Priority = None

        self.visited = False

        self.model: Optional[Model] = None
        self.model_filename: Optional[ModelSourceInfo] = None
        self.residue_mapping: Optional[Dict[int, Residue]] = None


class Builder:
    def __init__(self, rass_data: RassData) -> None:

        self.rass_data: RassData = rass_data
        self.seq: str = rass_data.seq

        self.chain_nuc_ids = utils.ss.parens_to_chains(
            rass_data.dot_bracket, skip_separator=False, start=1
        )

        ids_list = [aa for a in self.chain_nuc_ids for aa in a]

        self.next_dict = {}
        self.prev_dict = {}
        for l in self.chain_nuc_ids:
            for a, b in zip(l[:-1], l[1:]):
                self.next_dict[a] = b
                self.prev_dict[b] = a

        seq_mapping_list = [
            aa
            for a in utils.ss.parens_to_chains(
                rass_data.dot_bracket, skip_separator=False, start=0
            )
            for aa in a
        ]

        self.letters: Dict[int, str] = {
            a: self.seq[b] for (a, b) in zip(ids_list, seq_mapping_list)
        }

        self.valid_nuc_ids: Set[int] = set(ids_list)

        # For each nucleotide, a list of (module, residue_in_module)
        self.module_assignment: Dict[int, List[Tuple[Node, int]]] = {
            k: [] for k in self.valid_nuc_ids
        }
        self.nodes: List[Node] = []

        self.rigids: List["Rigid"] = []
        self.chain_of_nucleotides: Dict[int, "Nucleotide"] = {}
        self.rigidCounter: int = 0

        self.models_used: List[Tuple[List[int], ModelSourceInfo]] = []

    # TODO: perhaps add the nucs numbering in this function instead of Node
    # initialization
    def assign_module(self, module: Node) -> None:
        for i, n in enumerate(module.nucs):
            self.module_assignment[n].append((module, i))
        self.nodes.append(module)


def new_node_from_motif_info(info: MotifInfo) -> Node:
    return Node(("module", info.filename), info.nucs)


def generate_motif_nodes(builder: Builder) -> None:
    for info in builder.rass_data.user_motifs:
        n = new_node_from_motif_info(info)
        builder.nodes.append(n)
        for motif_nuc, target_nuc in enumerate(info.nucs):
            builder.module_assignment[target_nuc].append((n, motif_nuc))


def generate_nodes_from_rnamoip_components(builder: Builder) -> None:
    # TODO: Somehow take care of rnamoip components?
    # Thing is, mapping rnamoip style motifs requires to actually look at the
    # pdb file in order to see where there are gaps
    assert (
        len(builder.rass_data.rnamoip_components) == 0
    ), "rnamoip style input temporarily unavailable"
    return
    modules_by_name = {}
    rnamoip_components = rass_data.rnamoip_components.copy()
    for info in rnamoip_components:
        name = info.name
        if name in modules_by_name:
            module = modules_by_name[name]
        else:
            module = Node(("module", name))
            builder.nodes.append(module)
            modules_by_name[name] = module
        # print(name, comp, first, last)
        assert (
            len(module.components) == comp
        ), "Index of module component does not match order"
        module.components.append((first, last))
        for i in range(first, last + 1):
            builder.module_assignment[i].append((module, comp, i - first))
            module.nucls.add(i)


def have_overlapping_module(builder: Builder, nucls: Tuple[int, ...]) -> bool:
    """Checks whether nucls are already covered in their entirety by at least one existing node"""
    b = [{node.id for (node, _) in builder.module_assignment[id]} for id in nucls]
    a = set.intersection(*b)
    return len(a) >= 1


def generate_auto_nodes(builder: Builder, pairing: Dict[int, int]) -> None:
    """Generate and add nodes related to ncm's, loose ends, lonely pairs and
    lonely nucleotides"""
    # TODO: split this function in ncm, pairs, etc...
    # TODO: for ncms, check if they exist. Bulges and inner loops might not exist
    module_assignment = builder.module_assignment
    nodes = builder.nodes
    rass_data = builder.rass_data
    added_ncms = set()
    stacks = rass_data.stacks
    letters = builder.letters

    def make_ncm_node(shape, ids):
        assert all(type(i) is int for i in ids)
        return Node(("ncm", shape, "".join(letters[i] for i in ids)), ids)

    assert None not in pairing
    for i in sorted(pairing.keys()):
        j = pairing[i]
        ii = builder.next_dict.get(i)
        jj = builder.prev_dict.get(j)
        iii = builder.next_dict.get(ii)
        jjj = builder.prev_dict.get(jj)

        # Note the "assert ii is not None ..." is to appease pytype
        # TODO: perhaps there is a better way?
        # "Crazy" idea: make the segmentation generic such that stacks are
        # considered as a type of loop. Then this stuff won't be needed.
        if ii in pairing and pairing[ii] == jj:
            ts = tuple(sorted((i, ii, jj, j)))
            if ts in added_ncms or have_overlapping_module(builder, ts):
                continue
            added_ncms.add(ts)
            assert ii is not None and jj is not None
            module = make_ncm_node("2_2", [i, ii, jj, j])
            builder.assign_module(module)

        elif ii in pairing and pairing[ii] == jjj:
            ts = tuple(sorted((i, ii, jjj, jj, j)))
            if ts in added_ncms or have_overlapping_module(builder, ts):
                continue
            added_ncms.add(ts)
            assert ii is not None and jjj is not None and jj is not None
            module = make_ncm_node("2_2", [i, ii, jjj, jj, j])
            builder.assign_module(module)

        elif iii in pairing and pairing[iii] == jj:
            ts = tuple(sorted((i, iii, iii, jj, j)))
            if ts in added_ncms or have_overlapping_module(builder, ts):
                continue
            added_ncms.add(ts)
            assert ii is not None and jjj is not None and jj is not None
            module = make_ncm_node("2_3", [i, ii, iii, jj, j])

        elif iii in pairing and pairing[iii] == jjj:
            ts = tuple(sorted((i, ii, iii, jjj, jj, j)))
            if ts in added_ncms or have_overlapping_module(builder, ts):
                continue
            added_ncms.add(ts)
            assert (
                ii is not None
                and iii is not None
                and jjj is not None
                and jj is not None
            )
            module = make_ncm_node("2_3", [i, ii, iii, jjj, jj, j])
            builder.assign_module(module)

    for i, j in pairing.items():
        ts = tuple(sorted((i, j)))
        if ts in added_ncms or have_overlapping_module(builder, ts):
            continue
        added_ncms.add(ts)
        module = Node(("pair", letters[i] + letters[j]), [i, j])
        builder.assign_module(module)

    def add_helix_stack_module(i, j):
        assert j == i + 1
        ts = tuple(sorted((i, j)))
        if ts in added_ncms or have_overlapping_module(builder, ts):
            # TODO: check if the overlap is a helix fragment or not <- why?
            return
        added_ncms.add(ts)
        module = Node(("helix_stack", letters[i] + letters[j]), [i, j])
        builder.assign_module(module)

    # Add stacking that was given by the instructions
    for i, j in stacks:
        add_helix_stack_module(i, j)

    # Turn dangling ends into stacks
    fp_start = 0
    for ch in builder.chain_nuc_ids:
        if len(ch) == 0:
            continue
        for i in ch[:-1]:
            if len(module_assignment[i]) == 0:
                fp_start = i + 1
            else:
                break
        # print("fp_start", fp_start)
        for i in range(ch[0], fp_start):
            add_helix_stack_module(i, i + 1)

        tp_start = ch[-1]
        for i in reversed(ch[1:]):
            if len(module_assignment[i]) == 0:
                tp_start = i - 1
            else:
                break
        for i in range(ch[-1], tp_start, -1):
            add_helix_stack_module(i - 1, i)

    # Add single nucleotides
    for i in builder.valid_nuc_ids:
        if len(module_assignment[i]) == 0:
            ts = (i,)
            if ts in added_ncms or have_overlapping_module(builder, ts):
                continue
            added_ncms.add(ts)
            module = Node(("nucleotide", letters[i]), [i])
            builder.assign_module(module)


def assign_priorities(builder: Builder) -> None:
    for i, node in enumerate(builder.nodes):
        if node.kind[0] == "module":
            node.priority = i


import vpython as vp
import Bio.PDB as PDB


import gzip


# TODO: take the correspondence files from RNA-Puzzles assessment
def canonical_atom_name(name: str) -> str:
    if name == "O1P":
        return "OP1"
    if name == "O2P":
        return "OP2"
    return name.replace("*", "'")


def canonical_residue_name(name: str) -> str:
    # TODO
    # Note: in this function we won't rename modified nucleotide names
    return name


def canonical_unmodified_residue_name(name: str) -> str:
    # TODO
    # Note: in this function we *will* rename modified nucleotide names
    return name


from Bio.PDB.Atom import Atom
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

# Creates a new model that has canonical atom representation
def canonicalize_model(model: Model) -> Model:
    """Given a Biopython Model of RNA, returns a new Biopython Model that has
    all its atom names canonicalized"""
    output_model = PDB.Model.Model(0)
    for chain in model:
        output_chain = PDB.Chain.Chain(chain.id)
        output_model.add(output_chain)
        for residue in chain:
            output_residue = PDB.Residue.Residue(
                residue.id, canonical_residue_name(residue.resname), residue.segid
            )
            output_chain.add(output_residue)
            for atom in residue:
                can_name = canonical_atom_name(atom.name)
                full_can_name = can_name
                output_atom = PDB.Atom.Atom(
                    name=can_name,
                    coord=atom.coord,
                    bfactor=atom.bfactor,
                    occupancy=atom.occupancy,
                    altloc=atom.altloc,
                    fullname=full_can_name,
                    serial_number=atom.serial_number,
                    element=atom.element,
                )
                output_residue.add(output_atom)

    return output_model


def is_struct_filename(filename: str) -> bool:
    return (
        filename.endswith(".pdb")
        or filename.endswith(".pdb.gz")
        or filename.endswith(".cif")
        or filename.endswith(".cif.gz")
    )


import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning


def load_struct_file(filename: str, struct_name: Optional[str] = None) -> Structure:
    if struct_name is None:
        struct_name = filename

    with warnings.catch_warnings():  # (Silence specific warnings)
        warnings.simplefilter("ignore", category=PDBConstructionWarning)

        if filename.endswith(".pdb"):
            return PDB.PDBParser().get_structure(struct_name, filename)
        elif filename.endswith(".pdb.gz"):
            with gzip.open(filename, "rt") as f:
                return PDB.PDBParser().get_structure(struct_name, f)
        elif filename.endswith(".cif"):
            return PDB.MMCIFParser().get_structure(struct_name, filename)
        elif filename.endswith(".cif.gz"):
            with gzip.open(filename, "rt") as f:
                return PDB.MMCIFParser().get_structure(struct_name, f)

    assert False, "Not a proper filename"


import random


def choose_or_not(options: List[str]) -> str:
    """A hack to switch between deterministic and randomized fragment choice."""
    # return random.choice(options)
    o = sorted(options)
    return o[0]


from rna_bits.utils.data_path import DATA_PATH
from rna_bits.utils.data_path import get_path

LIBRARY_MOTIFS_DIR = os.path.join(DATA_PATH, "database")


def resolve_filename(filename: str, user_motifs_dir: str) -> str:
    """
    The idea is that if motif files specified in a .rass file start with "./",
    then search around the .rass file
    """
    if filename.startswith("./"):
        return os.path.join(user_motifs_dir, filename)
    if not filename.startswith("/"):
        return os.path.join(LIBRARY_MOTIFS_DIR, filename)
    return filename


def get_struct_paths_from_path(path: str) -> List[str]:
    """If path is a directory, return all the structure files within it.
    If path is a structure file, then return that file.
    """

    def allow_filename(f: str) -> bool:
        return is_struct_filename(f) and not f.startswith(".")

    if os.path.isdir(path):
        ret = [os.path.join(path, f) for f in os.listdir(path) if allow_filename(f)]
        assert ret, f"No struct files found in {path}"
        return ret
    else:
        assert is_struct_filename(path), f"Unknown file type {path}"
        return [path]


def load_model_from_path(resolved_path: str) -> Tuple[Model, str]:
    struct_paths = get_struct_paths_from_path(resolved_path)
    choice = choose_or_not(struct_paths)

    model = load_struct_file(choice)[0]
    model = canonicalize_model(model)
    return model, choice


def load_model_from_mcsym_db(shape: str, letters: str) -> Tuple[Model, str]:
    """Shape is like "2_2" and letters is like "AGCU"."""
    # Some folders in mcsym-db have lowercase letters (I'm not sure why.)
    # Because of that I manually iterate over
    struct_directory = None
    ncm_shape_dir = os.path.join(LIBRARY_MOTIFS_DIR, "mcsym-db/", shape)
    for folder in sorted(os.listdir(ncm_shape_dir)):
        if folder.upper() == letters:
            struct_directory = os.path.join(ncm_shape_dir, folder)
            break
    assert struct_directory is not None, f"Couldn't find files for {str}/{letters}"

    # if the theoretical model exists:
    theoretical_filename = os.path.join(
        struct_directory, "1-" + shape + "-" + letters + "_t.pdb.gz"
    )
    if os.path.isfile(theoretical_filename):
        chosen_path = theoretical_filename
    else:
        chosen_path = struct_directory

    return load_model_from_path(struct_directory)


def load_model_from_db(db: str, shape: str, letters: str) -> Tuple[Model, str]:
    path = os.path.join(LIBRARY_MOTIFS_DIR, db, shape, letters)
    return load_model_from_path(path)


def load_model_from_rosetta_db(shape: str, letters: str) -> Tuple[Model, str]:
    return load_model_from_db("rosetta_canonical", shape, letters)


def load_stack_model(letters: str) -> Tuple[Model, str]:
    return load_model_from_mcsym_db("2_2", letters)


def delete_residues_by_number(model: Model, to_delete: Container[int]) -> None:
    """If to_delete is, say, (0,3), the function will delete the 0th and 3rd
    residues in the model (regardless of the chain numbers and the residues'
    formal positions in the chains)."""
    i = 0
    for chain in model:
        residue_ids_to_delete = []
        for residue in chain:
            if i in to_delete:
                residue_ids_to_delete.append(residue.id)
            i += 1
        for id in residue_ids_to_delete:
            chain.detach_child(id)


def has_wobble_pair(letters: str) -> bool:
    return letters[0] + letters[3] in ("GU", "UG") or letters[1] + letters[2] in (
        "GU",
        "UG",
    )


# returns (model, filename)
# TODO: put the filename selection into a different file
# TODO: add a caching layer
# TODO: instead of passing the rass_filename, pass the directory to search for
# TODO: Try to simplify
def load_model_from_kind(
    kind: NodeKind, user_motifs_dir: str
) -> Tuple[Model, ModelSourceInfo]:
    if kind[0] == "module":
        resolved_path = resolve_filename(kind[1], user_motifs_dir)
        assert os.path.exists(
            resolved_path
        ), f"{repr(resolved_path)}, aka {repr(kind[1])}, does not exist"
        return load_model_from_path(resolved_path)
    if kind[0] == "ncm":
        shape = kind[1]
        letters = kind[2]
        if has_wobble_pair(letters):
            return load_model_from_db("rna_bits/canonical", shape, letters)
        return load_model_from_mcsym_db(shape, letters)
    if kind[0] == "pair":
        # Load a stack and truncate
        model, filename = load_stack_model(kind[1][0] + "GC" + kind[1][1])
        delete_residues_by_number(model, (1, 2))
        return model, (filename, "trimmed")
    if kind[0] == "helix_stack":
        complementary = {"A": "U", "U": "A", "G": "C", "C": "G"}
        a = kind[1][0]
        b = kind[1][1]
        c = complementary[b]
        d = complementary[a]
        # Load a stack and truncate
        model, filename = load_stack_model(a + b + c + d)
        delete_residues_by_number(model, (2, 3))
        return model, (filename, "trimmed")
    if kind[0] == "nucleotide":
        name = kind[1]

        directory = name + ".pdb"
        model = load_resource_structure(directory)[0]
        model = canonicalize_model(model)
        return model, directory
    assert False, "We don't recognize " + repr(kind)


# returns a map target_nucleotide_id -> residue_in_the_model
def map_residues(node: Node) -> Dict[int, Residue]:
    model = node.model
    assert model is not None

    pos_list = node.nucs
    pos_set = node.nucs_set

    res_list = []
    for chain in model:
        for residue in chain:
            res_list.append(residue)

    # print(node.kind)
    # print(pos_list, res_list)

    # TODO: perhaps not enforce this, for example to allow adding ions in the fragments
    assert len(pos_list) == len(res_list)

    ret = {}
    for pos, res in zip(pos_list, res_list):
        ret[pos] = res

    return ret


import importlib.resources as importlib_resources


def load_resource_structure(fname: str) -> Model:
    if __package__:
        with importlib_resources.open_text(__package__, fname) as f:
            return PDB.PDBParser().get_structure(fname, f)
    else:
        path = os.path.join(os.path.dirname(__file__), fname)
        with open(path) as f:
            return PDB.PDBParser().get_structure(fname, f)


def load_substitution_model() -> Model:
    return load_resource_structure("bases.pdb")[0]


# Model contains only the 4 bases, used for base substitution
# Each base is in its own chain, where the chain name is the name of the base.
SUBSTITUTION_MODEL: Model = load_substitution_model()


def get_reference_base(letter: str) -> Residue:
    """Returns (a reference to) a Residue representing the base of the given letter.
    Please don't modify it.
    """
    return SUBSTITUTION_MODEL[letter][1]


BB_SUGAR_ATOMS = {
    "C1'",
    "C2'",
    "O2'",
    "C3'",
    "O3'",
    "C4'",
    "O4'",
    "C5'",
    "O5'",
    "P",
    "OP1",
    "OP2",
}


def substitute_base(residue: Residue, newResname: str) -> None:
    oldResname = canonical_unmodified_residue_name(residue.resname)
    oldBaseAtoms = []
    referenceOldBase = get_reference_base(oldResname)
    fixed_list = []
    moving_list = []
    for atom in referenceOldBase:
        if atom.id in residue:
            fixed_list.append(residue[atom.id])
            moving_list.append(atom)

    superimposer = PDB.Superimposer()
    superimposer.set_atoms(fixed=fixed_list, moving=moving_list)

    # Remove the old base
    for atom in oldBaseAtoms:
        residue.detach_child(atom.id)

    # Remove everything that isn't part of our known atoms.
    # I use remove_list because I got the impression that removing atoms from a
    # Residue while iterating through it causes problems.
    remove_list = []
    for atom in residue:
        if atom.id not in BB_SUGAR_ATOMS:
            remove_list.append(atom.id)
    for id_to_remove in remove_list:
        residue.detach_child(id_to_remove)

    freshBaseAtoms = []
    for atom in get_reference_base(newResname):
        fresh_atom = atom.copy()
        residue.add(fresh_atom)
        freshBaseAtoms.append(fresh_atom)

    residue.resname = newResname

    superimposer.apply(freshBaseAtoms)


def assign_models(builder: Builder, user_motifs_dir: str) -> None:
    nodes = builder.nodes
    models_used = builder.models_used
    letters = builder.letters
    for node in nodes:
        # print(node.kind, node.nucs)
        node.model, model_filename = load_model_from_kind(node.kind, user_motifs_dir)
        node.model_filename = model_filename
        models_used.append((node.nucs, model_filename))
        node.residue_mapping = map_residues(node)
        # print(node.residue_mapping)

        did_substitute = False
        # TODO: instead of passing seq, add the sequence information into the Node object
        for k, v in node.residue_mapping.items():
            if letters[k].isupper() and letters[k] != canonical_residue_name(v.resname):
                did_substitute = True
                substitute_base(v, letters[k])


BB_ATOM_NAMES = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]


# Given two nucleotides of the same type,
# returns two lists of corresponding atoms
# (references to the existing atoms, not copies)
def correspond_nucleotide_atoms(
    a: Residue, b: Residue, only_bb: bool = False
) -> Tuple[List[Atom], List[Atom]]:
    # TODO: Make this nicer
    a_canon = {atom.name: atom for atom in a}
    b_canon = {atom.name: atom for atom in b}
    assert len(a_canon) == len(a)
    assert len(b_canon) == len(b)
    all_atoms = set(a_canon.keys()) | set(b_canon.keys())
    common_atoms = set(a_canon.keys()) & set(b_canon.keys())
    diff = all_atoms - common_atoms
    if diff:
        warnings.warn("Some atoms were not shared {diff}")
    a_ret = []
    b_ret = []
    for atom_name in common_atoms:
        if (not only_bb) or atom_name in BB_ATOM_NAMES:
            a_ret.append(a_canon[atom_name])
            b_ret.append(b_canon[atom_name])
    return (a_ret, b_ret)


def superimpose_nodes(fixed, moving, nucls):
    # TODO: this function isn't currently used; could we use it?
    fixed_list = []
    moving_list = []
    for n in nucls:
        fixed_res = fixed.residue_mapping[n]
        moving_res = moving.residue_mapping[n]
        fixed_atoms, moving_atoms = correspond_nucleotide_atoms(
            fixed_res, moving_res, False
        )
        fixed_list += fixed_atoms
        moving_list += moving_atoms

    superimposer = PDB.Superimposer()
    superimposer.set_atoms(fixed_list, moving_list)

    superimposer.apply(moving.model.get_atoms())


class Nucleotide:
    def __init__(
        self, model: Residue, rigid: "Rigid", pos: int, priority: Priority
    ) -> None:
        self.model = model
        self.rigid = rigid
        self.pos = pos
        self.priority = priority
        rigid.nucleotides.add(self)

    def __repr__(self) -> str:
        return "Nuc{}".format(self.pos + 1)

    def __lt__(self, other: "Nucleotide") -> bool:
        return self.pos < other.pos


class Rigid:
    id_increment = 1000

    def __init__(self) -> None:
        # self.translation = nunmpy.array([0.,0.,0.])
        # self.rotation = nunmpy.array([0.,0.,0.])
        self.nucleotides = set()
        self.parent = None
        self.id = Rigid.id_increment
        Rigid.id_increment += 1

    # Merge other into self
    def merge(self, other: "Rigid") -> None:
        if self is other:
            return
        for n in other.nucleotides:
            n.rigid = self
        self.nucleotides.update(other.nucleotides)
        other.nucleotides = None
        other.parent = self


def traverse_and_stack(builder: Builder, node: Node, currentRigid: Rigid) -> None:
    if node.visited:
        return
    print("Putting down", node.model_filename)
    node.visited = True

    fixed_residues: List[Residue] = []
    aligning_residues: List[Residue] = []

    to_place_ids: List[int] = []
    fixed_ids: List[int] = []

    for nuc_id in node.nucs:
        if nuc_id in builder.chain_of_nucleotides:
            fixed_residues.append(builder.chain_of_nucleotides[nuc_id].model)
            fixed_ids.append(nuc_id)
            # print(node.residue_mapping[nuc_id])
            aligning_residues.append(node.residue_mapping[nuc_id])
        else:
            to_place_ids.append(nuc_id)

    if len(fixed_residues) != 0:

        fixed_atoms = []
        aligning_atoms = []

        for fixed_res, aligning_res in zip(fixed_residues, aligning_residues):
            fixed_a, aligning_a = correspond_nucleotide_atoms(
                fixed_res, aligning_res, False
            )
            fixed_atoms += fixed_a
            aligning_atoms += aligning_a

        superimposer = PDB.Superimposer()
        superimposer.set_atoms(fixed_atoms, aligning_atoms)

        # Apply the superimposer transformation to the fragment
        assert node.residue_mapping is not None
        superimposer.apply(node.residue_mapping.values())

    for k in to_place_ids:
        v = node.residue_mapping[k]
        wrapper = Nucleotide(v.copy(), currentRigid, k, node.priority)
        builder.chain_of_nucleotides[k] = wrapper

    for k in fixed_ids:
        old_nucleotide = builder.chain_of_nucleotides[k]
        if is_lhs_higher_priority(node.priority, old_nucleotide.priority):
            v = node.residue_mapping[k]
            wrapper = Nucleotide(v.copy(), currentRigid, k, node.priority)
            builder.chain_of_nucleotides[k] = wrapper

    for nuc_id in node.nucs:
        for (other, _) in builder.module_assignment[nuc_id]:
            if not other.visited:
                traverse_and_stack(builder, other, currentRigid)


def assemble(builder: Builder) -> None:
    rigids = builder.rigids
    while True:
        # Find biggest node to place first
        central_node = None
        for node in builder.nodes:
            if not node.visited and (
                central_node is None or (len(node.nucs) > len(central_node.nucs))
            ):
                central_node = node
        if central_node is None:
            break

        # center the first model, just because
        # pytype: disable=attribute-error
        acc = sum(atom.coord for atom in central_node.model.get_atoms())
        tot = sum(1 for atom in central_node.model.get_atoms())
        # pytype: enable=attribute-error
        centroid = acc / tot
        builder.rigidCounter += 1

        currentRigid = Rigid()
        rigids.append(currentRigid)

        for atom in central_node.model.get_atoms():  # pytype: disable=attribute-error
            atom.coord -= centroid + builder.rigidCounter * 20.0

        # Actually assemble
        traverse_and_stack(builder, central_node, currentRigid)


def generate_c3prime_nuc(builder: Builder, coord, letter, pos):
    output_residue = PDB.Residue.Residue((" ", 0, " "), letter, "    ")
    # 'name', 'coord', 'bfactor', 'occupancy', 'altloc', 'fullname', 'serial_number', element=None
    output_atom = PDB.Atom.Atom("C3'", coord, 0.0, 1.0, " ", "C3'", 1, "C")
    output_residue.add(output_atom)
    rigid = Rigid()
    builder.rigids.append(rigid)
    nucWrapper = Nucleotide(output_residue, rigid, pos, None)
    return nucWrapper


def generate_unplaced_nucleotides(builder: Builder) -> None:
    chain_of_nucleotides = builder.chain_of_nucleotides
    # Find nucleotides that haven't been placed and turn them into rigids
    for i in builder.valid_nuc_ids:
        # assert i in chain_of_nucleotides
        if i not in chain_of_nucleotides:
            warnings.warn(f"{i} is not in chain of nucleotides for some reason")

            builder.rigidCounter += 1
            coord = np.array([0.0, 0, 0]) - builder.rigidCounter * 20.0
            chain_of_nucleotides[i] = generate_c3prime_nuc(
                builder, coord, builder.letters[i], i
            )


def is_rigid_edge_nuc(builder: Builder, i: int) -> bool:
    chain_of_nucleotides = builder.chain_of_nucleotides
    fp = builder.valid_nuc_ids - builder.prev_dict.keys()
    tp = builder.valid_nuc_ids - builder.next_dict.keys()
    if i in builder.next_dict:
        ii = builder.next_dict[i]
        if chain_of_nucleotides[ii].rigid is not chain_of_nucleotides[i].rigid:
            return True
    if i in builder.prev_dict:
        ii = builder.prev_dict[i]
        if chain_of_nucleotides[ii].rigid is not chain_of_nucleotides[i].rigid:
            return True
    return False


def check_cycle_not_degenerate(
    prevMap: Dict[Nucleotide, Tuple[Nucleotide, Nucleotide]], startNuc: Nucleotide
) -> bool:
    assert startNuc in prevMap
    c, d = prevMap[startNuc]
    a, b = prevMap[c]
    return not ((c is d) and (a is b) and a is startNuc)


# Distance between C3' atoms.
# Took from NAST
NUCLEOTIDE_DISTANCE = 5.78

# Some arbitrary distance
# NUCLEOTIDE_DISTANCE = 7.5

import heapq

from numpy import float32, float64, ndarray


def traverse_to_get_cycle(
    builder: Builder, startNuc: Nucleotide, endNuc: Nucleotide
) -> Tuple[Optional[List[Tuple[Nucleotide, Nucleotide]]], float64]:
    """
    The idea is to get the shortest path from startNuc to startNuc,
    but without passing the startNuc->endNuc link in that direction (but
    allowing the endNuc->startNuc direction), so that we find the shortest
    cycle containing that link.
    """
    chain_of_nucleotides = builder.chain_of_nucleotides

    prevMap = {}

    queue = [(0.0, startNuc, None)]

    score = None
    while len(queue):
        dist, nuc, prev = heapq.heappop(queue)

        if nuc in prevMap:
            continue

        if prev is not None:
            prevMap[nuc] = prev

        # print(dist, nuc, prev)

        if (
            nuc is startNuc
            and prev is not None
            and check_cycle_not_degenerate(prevMap, startNuc)
        ):
            prevMap[nuc] = prev
            score = dist
            break

        for nn in nuc.rigid.nucleotides:
            dist2 = dist + np.linalg.norm(
                get_c3prime(nn.model) - get_c3prime(nuc.model)
            )
            for j in (nn.pos + 1, nn.pos - 1):
                if j not in chain_of_nucleotides:
                    continue
                nnn = chain_of_nucleotides[j]
                if nnn.rigid is nn.rigid:
                    continue

                if nn is startNuc and nnn is endNuc:
                    continue

                dist3 = dist2 + NUCLEOTIDE_DISTANCE

                if nnn not in prevMap:
                    # print("nnnnn",nuc, nn, nnn)
                    heapq.heappush(queue, (dist3, nnn, (nuc, nn)))

    if score is None:
        return None, math.inf

    ret = []
    ret.append(prevMap[startNuc])
    current = prevMap[startNuc][0]
    while current is not startNuc:
        ret.append(prevMap[current])
        current = prevMap[current][0]

    ret.reverse()

    return ret, score


# Finds the next we want to process
def get_next_cycle(builder: Builder) -> Optional[List[Tuple[Nucleotide, Nucleotide]]]:
    chain_of_nucleotides = builder.chain_of_nucleotides
    # print(chain_of_nucleotides)
    bestScore = math.inf
    bestCycle = None
    for i in builder.valid_nuc_ids:
        if is_rigid_edge_nuc(builder, i):
            nuc = chain_of_nucleotides[i]

            for j in (i - 1, i + 1):
                if j not in chain_of_nucleotides:
                    continue
                nn = chain_of_nucleotides[j]
                if nn.rigid is nuc.rigid:
                    continue

                cycle, score = traverse_to_get_cycle(builder, nuc, nn)
                if score < bestScore:
                    bestScore = score
                    bestCycle = cycle
    return bestCycle


def is_terminal_rigid(builder: Builder, r: Rigid) -> bool:
    # Get all the terminal nucleotides
    # TODO: do this only once
    fp = builder.valid_nuc_ids - builder.prev_dict.keys()
    tp = builder.valid_nuc_ids - builder.next_dict.keys()
    term_nucs = fp | tp
    print("term_nucs", term_nucs)

    return r in (builder.chain_of_nucleotides[i].rigid for i in term_nucs)


def traverse_to_get_path(
    builder: Builder, startNuc: Nucleotide
) -> List[Tuple[Nucleotide, Nucleotide]]:
    """Traverse to get a path between two rigids that are connected to only one edge each
    (to get an "outer loop")
    """
    chain_of_nucleotides = builder.chain_of_nucleotides

    assert is_terminal_rigid(builder, startNuc.rigid)
    assert is_rigid_edge_nuc(builder, startNuc.pos)

    visitedRigids = set()
    out = []
    nuc = startNuc
    visitedRigids.add(nuc.rigid)
    while nuc is startNuc or not is_terminal_rigid(builder, nuc.rigid):
        # print(nuc)
        # print(out)

        result = None

        for nn in nuc.rigid.nucleotides:
            # dist2 = dist + np.linalg.norm(get_c3prime(nn.model) - get_c3prime(nuc.model))
            # print(nuc, nn)
            for j in (nn.pos + 1, nn.pos - 1):
                if j not in chain_of_nucleotides:
                    continue
                nnn = chain_of_nucleotides[j]
                # print(nuc, nn, nnn)
                if nnn.rigid is nn.rigid:
                    continue
                # dist3 = dist2 + NUCLEOTIDE_DISTANCE
                if nnn.rigid not in visitedRigids:
                    result = (nuc, nn, nnn)
                    break
            if result is not None:
                break

        assert result is not None

        nuc, nn, nnn = result
        out.append((nuc, nn))
        nuc = nnn
        visitedRigids.add(nuc)
    out.append((nuc, nuc))
    return out


def get_next_path(builder: Builder) -> Optional[List[Tuple[Nucleotide, Nucleotide]]]:
    """Assumes that all cycles have been found.
    Finds the next "outer loop"
    """
    for i in builder.valid_nuc_ids:
        if is_rigid_edge_nuc(builder, i):
            nuc = builder.chain_of_nucleotides[i]
            return traverse_to_get_path(builder, nuc)

    return None


from . import inscribed_polygon


def get_c3prime(residue: Residue) -> Atom:
    if "C3'" in residue:
        return residue["C3'"]
    if "C3*" in residue:
        return residue["C3*"]
    return None


def np_distance(a_coord: ndarray, b_coord: ndarray) -> Union[float32, float64]:
    return np.linalg.norm(a_coord - b_coord)


import math


def ortho_project(s: ndarray, d: ndarray) -> ndarray:
    return (s.dot(d) / (np.linalg.norm(d) ** 2)) * d


def create_rotation(s1: ndarray, s2: ndarray, d1: ndarray, d2: ndarray) -> ndarray:
    """Create a rotation matrix that transforms the vector s1 to align with the vector d1
    while having the vector s2 be aligned with the vector d2 as much as possible"""

    # Construct two orthonormal bases
    s1 = s1 / np.linalg.norm(s1)
    d1 = d1 / np.linalg.norm(d1)

    # print("s1 s2", s1, s2)
    # print("d1 d2", d1, d2)
    s2 = s2 - ortho_project(s2, s1)
    d2 = d2 - ortho_project(d2, d1)

    # print("s1 s2", s1, s2)
    # print("d1 d2", d1, d2)
    s2 = s2 / np.linalg.norm(s2)
    d2 = d2 / np.linalg.norm(d2)

    s3 = np.cross(s1, s2)
    d3 = np.cross(d1, d2)

    S = np.array([s1, s2, s3])
    D = np.array([d1, d2, d3])

    # Solve SX = D
    T = np.linalg.solve(S, D)
    T = T.transpose()

    # print(s1, d1, T.dot(s1))
    return T


def build_cycles(builder: Builder) -> None:
    while True:
        cycle = get_next_cycle(builder)

        if cycle is not None:
            # print("Cycle", cycle)
            isPath = False
        else:
            # for i,nuc in sorted(chain_of_nucleotides.items()):
            # print(i+1, nuc.rigid.id)
            cycle = get_next_path(builder)
            # print("Chain", cycle)
            isPath = True

        if cycle is None:
            break

        rigid_distances = []
        for (a, b) in cycle:
            a_atom = get_c3prime(a.model)
            b_atom = get_c3prime(b.model)
            rigid_distances.append(np_distance(a_atom.coord, b_atom.coord))

        distances = []
        for d in rigid_distances:
            distances.append(d)
            distances.append(NUCLEOTIDE_DISTANCE)

        if isPath:
            distances.pop()

        # print(distances)

        # To simplify things, if what I have is a path, I'll place it on a half-circle
        # It's guaranteed that a half circle is not going to be degenerate
        # (But the path needs to have at least 3 nodes, otherwise it's not really
        # possible to put it on a half-circle... (I think?) )
        is_degenerate = not isPath and (sum(distances) <= 2 * max(distances))
        if is_degenerate:
            # the polygon is going to be degenerate or impossible
            # might need to stretch the distance between nucleotides
            assert max(distances) != NUCLEOTIDE_DISTANCE
            adj_nuc_dist = NUCLEOTIDE_DISTANCE + (
                2 * max(distances) - sum(distances)
            ) / len(cycle)
            print("stretch distance to", adj_nuc_dist)
            pos = 0.0
            points = []
            for d in rigid_distances:
                points.append((pos, 0.0))
                if d == max(distances):
                    pos -= d
                else:
                    pos += d

                points.append((pos, 0.0))

                pos += adj_nuc_dist
                print("bridge", points)
        else:
            if not isPath:
                points = inscribed_polygon.construct_polygon(distances)
            else:
                points = inscribed_polygon.construct_polygon(distances, math.tau * 0.33)
                points.append(points[-1])

        points = [np.array([x, y, 0]) for (x, y) in points]
        # print(points)

        # place the pieces onto the polygon
        assert len(points) % 2 == 0
        for i in range(len(points) // 2):
            dpa = points[i * 2]
            dpb = points[i * 2 + 1]
            ra = cycle[i][0]
            rb = cycle[i][1]
            assert ra.rigid is rb.rigid

            spa = get_c3prime(ra.model).coord
            spb = get_c3prime(rb.model).coord

            if (spb == spa).all():
                sv1 = np.array([1.0, 0, 0])
                if "P" in ra.model:
                    if (i == 0 and cycle[i + 1][0].pos > cycle[i][1].pos) or cycle[i][
                        0
                    ].pos > cycle[i - 1][1].pos:
                        sv1 = spa - ra.model["P"].coord
                    else:
                        sv1 = ra.model["P"].coord - spa
            else:
                sv1 = spb - spa

            if is_degenerate:
                if rigid_distances[i] != max(rigid_distances):
                    dv1 = np.array([1, 0, 0])
                    dv2 = np.array([0, 1, 0])
                else:
                    dv1 = np.array([-1, 0, 0])
                    dv2 = np.array([0, -1, 0])

            else:
                if (spb == spa).all():
                    dv1 = np.array([-dpa[1], dpb[0], 0])
                else:
                    dv1 = dpb - dpa

                dv2 = (dpa + dpb) / 2
                # For the case of bridges, inverse the side that goes in the other direction
                # To check if the side goes in the other direction, compute CCW
                if (dpa[0] * dpb[1] - dpa[1] * dpb[0]) < 0:
                    dv2 *= -1

            # Find the centroid of the rigid to orient it away from the loop
            acc = sum(
                sum(atom.coord for atom in nuc.model) for nuc in ra.rigid.nucleotides
            )
            num = sum(sum(1 for atom in nuc.model) for nuc in ra.rigid.nucleotides)
            centroid = acc / num

            # print("centroid", centroid, spa)

            # if it is colinear, switch it to something else
            sv2 = centroid - spa
            if (
                np.linalg.norm(sv2) <= 1e-10
                or abs(
                    abs(sv1.dot(sv2)) / np.linalg.norm(sv1) / np.linalg.norm(sv2) - 1
                )
                <= 1e-10
            ):
                sv2 = np.array([1, 0, 0])
            if (
                abs(abs(sv1.dot(sv2)) / np.linalg.norm(sv1) / np.linalg.norm(sv2) - 1)
                <= 1e-10
            ):
                sv2 = np.array([0, 1, 0])

            rotation_matrix = create_rotation(sv1, sv2, dv1, dv2)
            # print(rotation_matrix)

            for nuc in ra.rigid.nucleotides:
                for atom in nuc.model:
                    coord = atom.coord
                    atom.coord = rotation_matrix.dot((coord - spa)) + dpa

        # combine all rigid bodies
        first = cycle[0][0]
        for (a, b) in cycle[1:]:
            assert a.rigid is not first
            assert a.rigid is b.rigid
            first.rigid.merge(a.rigid)


def visualize_stuff(builder: Builder):
    # visualize stuff
    # Used for debug purposes a long time ago
    # TODO: Does it even work?
    import vpython as vp

    def toVpVec(a):
        return vp.vector(a[0], a[1], a[2])

    for k, residue in builder.chain_of_nucleotides.items():
        for atom in residue.model:
            vp.sphere(pos=toVpVec(atom.coord), radius=0.3, color=vp.vector(0, 0.5, 0))
            if canonical_atom_name(atom.name) == "P":
                vp.sphere(
                    pos=toVpVec(atom.coord), radius=0.5, color=vp.vector(1, 0.25, 0)
                )


# Now all we need is to write a pdb output!
def generate_biopython_structure(builder: Builder) -> Structure:
    chain_of_nucleotides = builder.chain_of_nucleotides

    output_structure = PDB.Structure.Structure("output")
    output_model = PDB.Model.Model(0)
    output_chain = PDB.Chain.Chain("A")

    atom_serial_inc = 1
    for (i, rw) in sorted(chain_of_nucleotides.items()):
        r = rw.model
        output_residue = PDB.Residue.Residue((" ", i, " "), r.resname, "    ")
        for atom in r:
            # 'name', 'coord', 'bfactor', 'occupancy', 'altloc', 'fullname', and 'serial_number'
            output_atom = PDB.Atom.Atom(
                atom.name,
                atom.coord,
                atom.bfactor,
                atom.occupancy,
                atom.altloc,
                atom.fullname,
                atom_serial_inc,
                atom.element,
            )
            output_residue.add(output_atom)
        output_chain.add(output_residue)

    output_model.add(output_chain)
    output_structure.add(output_model)
    return output_structure


def get_output_filename(in_filename: Optional[str]) -> str:
    if in_filename is None:
        return "out.pdb"
    else:
        return in_filename + ".pdb"


def write_output_pdb_file(output_structure: Structure, out_filename: str) -> None:
    io = PDB.PDBIO()
    io.set_structure(output_structure)

    io.save(out_filename)
    print("Wrote to", repr(out_filename))


def process_rass_filename(rass_filename):
    if rass_filename is not None:
        f = open(rass_filename)
    else:
        f = sys.stdin

    rass_data = parse_rass_file(f)
    f.close()
    pairing = parse_parens(rass_data.dot_bracket)

    builder = Builder(rass_data)

    generate_motif_nodes(builder)
    generate_nodes_from_rnamoip_components(builder)
    generate_auto_nodes(builder, pairing)

    user_motifs_dir: str
    if rass_filename is not None:
        user_motifs_dir = os.path.dirname(rass_filename)
    else:
        user_motifs_dir = os.getcwd()

    assign_models(builder, user_motifs_dir)
    assign_priorities(builder)

    assemble(builder)
    generate_unplaced_nucleotides(
        builder
    )  # Probably unneeded because all lonely nucleotides should have been generated

    # for _, b in sorted(builder.chain_of_nucleotides.items()):
    #     print(b, b.rigid)

    build_cycles(builder)

    # rigids_seen = set()
    # for (_, nuc) in sorted(builder.chain_of_nucleotides.items()):
    #     rigids_seen.add(nuc.rigid)
    #     print(nuc, nuc.priority)

    # print(len(rigids_seen))

    if False:
        visualize_stuff()
    for a in builder.models_used:
        print(a)

    output_structure = generate_biopython_structure(builder)
    out_filename = get_output_filename(rass_filename)
    write_output_pdb_file(output_structure, out_filename)


def main(argv: List[str]) -> None:
    rass_filenames: List[str] = get_input_filenames(argv)

    if len(rass_filenames) == 0:
        process_rass_filename(None)
    else:
        for i, rass_filename in enumerate(rass_filenames):
            if len(rass_filenames) > 1:
                if i > 0:
                    print()
                print("Processing", rass_filename)
            process_rass_filename(rass_filename)


def main_wrapper():
    sys.setrecursionlimit(10000)  # Was testing something; probably remove
    main(sys.argv)


if __name__ == "__main__":
    main_wrapper()
