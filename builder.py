import sys
import numpy as np
from dataclasses import dataclass
from typing import Optional, Union, List, Tuple, Dict, TextIO, Literal
import io

# internally, I use 0-based indexing for sequences


def get_input_filename(argv: List[str]) -> Optional[str]:
    if len(sys.argv) >= 2:
        return sys.argv[1]
    else:
        return None


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
            first -= 1
            last -= 1
            comp -= 1
            rnamoip_components.append(
                RnamoipCompInfo(name=name, comp_id=comp, first_nuc=first, last_nuc=last)
            )
        elif line.startswith("motif:"):
            _, name, ranges = line.split(":")
            nucs = parse_and_expand_ranges(ranges)
            nucs = [n - 1 for n in nucs]
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
            a -= 1
            b -= 1
            for i in range(a, b):
                stacks.append((i, i + 1))

    f.close()
    return RassData(
        seq=seq,
        original_seq=original_seq,
        dot_bracket=dot_bracket,
        user_motifs=user_motifs,
        stacks=stacks,
        rnamoip_components=rnamoip_components,
    )


import Project.parser


def parse_parens(dot_bracket: str) -> List[Optional[int]]:
    pairing = [None] * len(dot_bracket)

    pos = 0

    for a, b in Project.parser.parseParens(dot_bracket):
        a -= 1
        b -= 1
        pairing[a] = b
        pairing[b] = a

    return pairing


NodeKind = Union[
    Tuple[Literal["ncm"], str, str],
    Tuple[Literal["module"], str],
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

        self.visited = False

        self.model: Optional[Model] = None
        self.model_filename: Optional[ModelSourceInfo] = None
        self.residue_mapping: Optional[Dict[int, Residue]] = None


class Builder:
    def __init__(self, rass_data: RassData) -> None:

        self.rass_data: RassData = rass_data
        self.seq: str = rass_data.seq

        # For each nucleotide, a list of (module, residue_in_module)
        self.module_assignment: List[List[Tuple[Node, int]]] = [
            [] for _ in range(len(rass_data.seq))
        ]
        self.nodes: List[Node] = []

        self.rigids: List["Rigid"] = []
        self.chain_of_nucleotides: Dict[int, "Nucleotide"] = {}
        self.rigidCounter: int = 0

        self.models_used: List[Tuple[List[int], ModelSourceInfo]] = []


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
        print(name, comp, first, last)
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


def generate_auto_nodes(builder: Builder, pairing: List[Optional[int]]) -> None:
    """Generate and add nodes related to ncm's, loose ends, lonely pairs and
    lonely nucleotides"""
    # TODO: split this function
    module_assignment = builder.module_assignment
    nodes = builder.nodes
    rass_data = builder.rass_data
    added_ncms = set()
    seq = rass_data.seq
    stacks = rass_data.stacks
    for i in range(0, len(seq)):
        if pairing[i] is not None:
            j = pairing[i]
            ii = i + 1
            jj = j - 1
            iii = i + 2
            jjj = j - 2
            module = None
            if ii < len(pairing) and pairing[ii] == jj:
                ts = tuple(sorted((i, ii, jj, j)))
                if ts in added_ncms or have_overlapping_module(builder, ts):
                    continue
                added_ncms.add(ts)
                module = Node(
                    ("ncm", "2_2", seq[i] + seq[ii] + seq[jj] + seq[j]), [i, ii, jj, j]
                )
                module_assignment[i].append((module, 0))
                module_assignment[ii].append((module, 1))
                module_assignment[jj].append((module, 2))
                module_assignment[j].append((module, 3))

            elif ii < len(pairing) and pairing[ii] == jjj:
                ts = tuple(sorted((i, ii, jjj, jj, j)))
                if ts in added_ncms or have_overlapping_module(builder, ts):
                    continue
                added_ncms.add(ts)
                module = Node(
                    ("ncm", "2_3", seq[i] + seq[ii] + seq[jjj] + seq[jj] + seq[j]),
                    [i, ii, jjj, jj, j],
                )
                module_assignment[i].append((module, 0))
                module_assignment[ii].append((module, 1))
                module_assignment[jjj].append((module, 2))
                module_assignment[jj].append((module, 3))
                module_assignment[j].append((module, 4))

            elif iii < len(pairing) and pairing[iii] == jj:
                ts = tuple(sorted((i, iii, iii, jj, j)))
                if ts in added_ncms or have_overlapping_module(builder, ts):
                    continue
                added_ncms.add(ts)
                module = Node(
                    ("ncm", "3_2", seq[i] + seq[ii] + seq[iii] + seq[jj] + seq[j]),
                    [i, ii, iii, jj, j],
                )
                module_assignment[i].append((module, 0))
                module_assignment[ii].append((module, 1))
                module_assignment[iii].append((module, 2))
                module_assignment[jj].append((module, 3))
                module_assignment[j].append((module, 4))
            elif iii < len(pairing) and pairing[iii] == jjj:
                ts = tuple(sorted((i, ii, iii, jjj, jj, j)))
                if ts in added_ncms or have_overlapping_module(builder, ts):
                    continue
                added_ncms.add(ts)
                module = Node(
                    (
                        "ncm",
                        "3_3",
                        seq[i] + seq[ii] + seq[iii] + seq[jjj] + seq[jj] + seq[j],
                    ),
                    [i, ii, iii, jjj, jj, j],
                )
                module_assignment[i].append((module, 0))
                module_assignment[ii].append((module, 1))
                module_assignment[iii].append((module, 2))
                module_assignment[jjj].append((module, 3))
                module_assignment[jj].append((module, 4))
                module_assignment[j].append((module, 5))

            if module:
                nodes.append(module)

    # Add isolated pairs
    for i in range(0, len(seq)):
        if pairing[i] is not None:
            j = pairing[i]
            ts = tuple(sorted((i, j)))
            if ts in added_ncms or have_overlapping_module(builder, ts):
                continue
            added_ncms.add(ts)
            module = Node(("pair", seq[i] + seq[j]), [i, j])
            module_assignment[i].append((module, 0))
            module_assignment[j].append((module, 1))
            nodes.append(module)

    def add_helix_stack_module(i, j):
        assert j == i + 1
        ts = tuple(sorted((i, j)))
        if ts in added_ncms or have_overlapping_module(builder, ts):
            # TODO: check if the overlap is a helix fragment or not
            return
        added_ncms.add(ts)
        module = Node(("helix_stack", seq[i] + seq[j]), [i, j])
        module_assignment[i].append((module, 0))
        module_assignment[j].append((module, 1))
        nodes.append(module)

    # Add stacking that was given by the instructions
    for i, j in stacks:
        add_helix_stack_module(i, j)

    # Turn dangling ends into stacks
    fp_start = 0
    for i in range(0, len(seq) - 1):
        if len(module_assignment[i]) == 0:
            fp_start = i + 1
        else:
            break
    print("fp_start", fp_start)
    for i in range(0, fp_start):
        add_helix_stack_module(i, i + 1)

    tp_start = len(seq) - 1
    for i in range(len(seq) - 1, 0, -1):
        if len(module_assignment[i]) == 0:
            tp_start = i - 1
        else:
            break
    for i in range(len(seq) - 1, tp_start, -1):
        add_helix_stack_module(i - 1, i)

    # Add single nucleotides
    for i in range(0, len(seq)):
        if len(module_assignment[i]) == 0:
            ts = (i,)
            if ts in added_ncms or have_overlapping_module(builder, ts):
                continue
            added_ncms.add(ts)
            module = Node(("nucleotide", seq[i]), [i])
            module_assignment[i].append((module, 0))
            nodes.append(module)


import vpython as vp
import Bio.PDB as PDB

import os

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


def load_struct_file(filename: str, struct_name: Optional[str] = None) -> Structure:
    if struct_name is None:
        struct_name = filename

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


# returns (model, filename)
# TODO: put the filename selection into a different file
# TODO: add a caching layer
# TODO: instead of passing the rass_filename, pass the directory to search for
# models in the "./" case
def load_model_kind(
    kind: NodeKind, rass_filename: Optional[str]
) -> Tuple[Model, ModelSourceInfo]:
    if kind[0] == "module":
        name = kind[1]
        directory = name
        if directory.startswith("./") and rass_filename is not None:
            directory = os.path.join(os.path.dirname(rass_filename), directory)
        if os.path.isdir(directory):
            options = [
                os.path.join(directory, x)
                for x in os.listdir(directory)
                if is_struct_filename(x) and not x.startswith(".")
            ]
            choice = choose_or_not(options)
        else:
            choice = directory
        model = load_struct_file(choice)[0]
        model = canonicalize_model(model)
        return model, choice
    if kind[0] == "ncm":
        directory = None
        for folder in sorted(os.listdir("mcsym-db/" + kind[1] + "/")):
            if folder.upper() == kind[2]:
                directory = "mcsym-db/" + kind[1] + "/" + folder
                break

        assert directory is not None

        # if the theoretical model exists:
        theoretical_filename = os.path.join(
            directory, "1-" + kind[1] + "-" + kind[2] + "_t.pdb.gz"
        )
        if os.path.isfile(theoretical_filename):
            filename = theoretical_filename
        else:
            options = [
                os.path.join(directory, x)
                for x in os.listdir(directory)
                if (x.endswith(".pdb.gz") and not x.startswith("."))
            ]
            filename = choose_or_not(options)
        model = load_struct_file(filename)[0]
        model = canonicalize_model(model)
        return model, filename
    if kind[0] == "pair":
        # load 2_2 from mcsym-db, and truncate
        model, filename = load_model_kind(
            ("ncm", "2_2", kind[1][0] + "GC" + kind[1][1]), None
        )
        i = 0
        for chain in model:
            residue_ids_to_delete = []
            for residue in chain:
                if i in (1, 2):
                    residue_ids_to_delete.append(residue.id)
                i += 1
            for id in residue_ids_to_delete:
                chain.detach_child(id)

        model = canonicalize_model(model)
        return model, (filename, "trimmed")
    if kind[0] == "helix_stack":
        complementary = {"A": "U", "U": "A", "G": "C", "C": "G"}
        a = kind[1][0]
        b = kind[1][1]
        c = complementary[b]
        d = complementary[a]
        # load 2_2 from mcsym-db, and truncate
        model, filename = load_model_kind(("ncm", "2_2", a + b + c + d), None)
        i = 0
        for chain in model:
            residue_ids_to_delete = []
            for residue in chain:
                if i >= 2:
                    residue_ids_to_delete.append(residue.id)
                i += 1
            for id in residue_ids_to_delete:
                chain.detach_child(id)

        model = canonicalize_model(model)
        return model, (filename, "trimmed")
    if kind[0] == "nucleotide":
        name = kind[1]
        directory = name + ".pdb"
        model = load_struct_file(directory)[0]
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

    # TODO: perhaps not enforce this, for example to allow adding ions in the fragments
    assert len(pos_list) == len(res_list)

    ret = {}
    for pos, res in zip(pos_list, res_list):
        ret[pos] = res

    return ret


def load_substitution_model() -> Model:
    return PDB.PDBParser().get_structure("substitution_model", "bases.pdb")[0]


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
    # I use remove_list because I got the impression that removing stuff while
    # iterating causes problems.
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


def assign_models(builder: Builder, seq: str, rass_filename: str) -> None:
    nodes = builder.nodes
    models_used = builder.models_used
    for node in nodes:
        print(node.kind, node.nucs)
        node.model, model_filename = load_model_kind(node.kind, rass_filename)
        node.model_filename = model_filename
        models_used.append((node.nucs, model_filename))
        node.residue_mapping = map_residues(node)
        print(node.residue_mapping)

        did_substitute = False
        # TODO: instead of passing seq, add the sequence information with the Node object
        for k, v in node.residue_mapping.items():
            if seq[k].isupper() and seq[k] != canonical_residue_name(v.resname):
                did_substitute = True
                substitute_base(v, seq[k])


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
        print("Some atoms were not shared", diff)
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
    def __init__(self, model: Residue, rigid: "Rigid", pos: int) -> None:
        self.model = model
        self.rigid = rigid
        self.pos = pos
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

    fixed_residues = []
    aligning_residues = []
    to_place_id = []
    for nuc_id in node.nucs:
        if nuc_id in builder.chain_of_nucleotides:
            fixed_residues.append(builder.chain_of_nucleotides[nuc_id].model)
            print(node.residue_mapping[nuc_id])
            aligning_residues.append(node.residue_mapping[nuc_id])
        else:
            to_place_id.append(nuc_id)

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

        assert node.residue_mapping is not None
        superimposer.apply(node.residue_mapping.values())

    for k in to_place_id:
        v = node.residue_mapping[k]
        wrapper = Nucleotide(v.copy(), currentRigid, k)
        builder.chain_of_nucleotides[k] = wrapper

    for nuc_id in node.nucs:
        for (other, _) in builder.module_assignment[nuc_id]:
            if not other.visited:
                traverse_and_stack(builder, other, currentRigid)


def assemble(builder: Builder) -> None:
    rigids = builder.rigids
    while True:
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
        traverse_and_stack(builder, central_node, currentRigid)


def generate_c3prime_nuc(builder: Builder, coord, letter, pos):
    output_residue = PDB.Residue.Residue((" ", 0, " "), letter, "    ")
    # 'name', 'coord', 'bfactor', 'occupancy', 'altloc', 'fullname', 'serial_number', element=None
    output_atom = PDB.Atom.Atom("C3'", coord, 0.0, 1.0, " ", "C3'", 1, "C")
    output_residue.add(output_atom)
    rigid = Rigid()
    builder.rigids.append(rigid)
    nucWrapper = Nucleotide(output_residue, rigid, pos)
    return nucWrapper


def generate_unplaced_nucleotides(builder: Builder, seq: str) -> None:
    chain_of_nucleotides = builder.chain_of_nucleotides
    # Find nucleotides that haven't been placed and turn them into rigids
    for i in range(0, len(seq)):
        # assert i in chain_of_nucleotides
        if i not in chain_of_nucleotides:
            print(i, "is not in chain of nucleotides for some reason")

            builder.rigidCounter += 1
            coord = np.array([0.0, 0, 0]) - builder.rigidCounter * 20.0
            chain_of_nucleotides[i] = generate_c3prime_nuc(builder, coord, seq[i], i)


def is_rigid_edge_nuc(builder: Builder, i: int) -> bool:
    chain_of_nucleotides = builder.chain_of_nucleotides
    if (i + 1) in chain_of_nucleotides:
        if chain_of_nucleotides[i + 1].rigid is not chain_of_nucleotides[i].rigid:
            return True
    if (i - 1) in chain_of_nucleotides:
        if chain_of_nucleotides[i - 1].rigid is not chain_of_nucleotides[i].rigid:
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
    seq = builder.seq
    chain_of_nucleotides = builder.chain_of_nucleotides
    print(chain_of_nucleotides)
    bestScore = math.inf
    bestCycle = None
    for i in range(len(seq)):
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
    # TODO: make this work even when there's more than 1 strand
    return r in (
        builder.chain_of_nucleotides[len(builder.seq) - 1].rigid,
        builder.chain_of_nucleotides[0].rigid,
    )


def traverse_to_get_path(
    builder: Builder, startNuc: Nucleotide
) -> List[Tuple[Nucleotide, Nucleotide]]:
    """Traverse to get a path between two rigids that are connected to only one edge each
    (to get an "outer loop")
    """
    seq = builder.seq
    chain_of_nucleotides = builder.chain_of_nucleotides

    assert is_terminal_rigid(builder, startNuc.rigid)
    assert is_rigid_edge_nuc(builder, startNuc.pos)

    visitedRigids = set()
    out = []
    nuc = startNuc
    visitedRigids.add(nuc.rigid)
    while nuc is startNuc or not is_terminal_rigid(builder, nuc.rigid):
        print(nuc)
        print(out)

        result = None

        for nn in nuc.rigid.nucleotides:
            # dist2 = dist + np.linalg.norm(get_c3prime(nn.model) - get_c3prime(nuc.model))
            print(nuc, nn)
            for j in (nn.pos + 1, nn.pos - 1):
                if j not in chain_of_nucleotides:
                    continue
                nnn = chain_of_nucleotides[j]
                print(nuc, nn, nnn)
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
    for i in range(len(builder.seq)):
        if is_rigid_edge_nuc(builder, i):
            nuc = builder.chain_of_nucleotides[i]
            return traverse_to_get_path(builder, nuc)

    return None


import inscribed_polygon


def get_c3prime(residue: Residue) -> Atom:
    if "C3'" in residue:
        return residue["C3'"]
    if "C3*" in residue:
        return residue["C3*"]
    return None


def np_distance(a_coord: ndarray, b_coord: ndarray) -> Union[float32, float64]:
    return np.linalg.norm(a_coord - b_coord)


import math
from io import TextIOWrapper


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
            print("Cycle", cycle)
            isPath = False
        else:
            # for i,nuc in sorted(chain_of_nucleotides.items()):
            # print(i+1, nuc.rigid.id)
            cycle = get_next_path(builder)
            print("Chain", cycle)
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

        print(distances)

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
        output_residue = PDB.Residue.Residue((" ", i + 1, " "), r.resname, "    ")
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


def main(argv: List[str]) -> None:
    rass_filename = get_input_filename(argv)
    if rass_filename is not None:
        f = open(rass_filename)
    else:
        f = sys.stdin
    rass_data = parse_rass_file(f)
    pairing = parse_parens(rass_data.dot_bracket)

    builder = Builder(rass_data)

    generate_motif_nodes(builder)
    generate_nodes_from_rnamoip_components(builder)
    generate_auto_nodes(builder, pairing)

    assign_models(builder, builder.rass_data.seq, rass_filename)
    assemble(builder)
    generate_unplaced_nucleotides(
        builder, rass_data.seq
    )  # Probably unneeded because all lonely nucleotides should have been generated
    for _, b in sorted(builder.chain_of_nucleotides.items()):
        print(b, b.rigid)

    build_cycles(builder)

    rigids_seen = set()
    for nuc in builder.chain_of_nucleotides.values():
        rigids_seen.add(nuc.rigid)
    print("rigids seen", [r.id for r in rigids_seen])

    if False:
        visualize_stuff()
    for a in builder.models_used:
        print(a)

    output_structure = generate_biopython_structure(builder)
    out_filename = get_output_filename(rass_filename)
    write_output_pdb_file(output_structure, out_filename)


if __name__ == "__main__":
    main(sys.argv)
