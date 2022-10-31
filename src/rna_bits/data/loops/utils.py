import sys, os
from Bio import PDB


def _get_conv(d):
    with open(d) as f:
        return dict(
            l.strip().split()
            for l in f
            if not l.startswith("#") and len(l.strip()) != 0
        )


res_dict = _get_conv("data/residues.list")
atom_dict = _get_conv("data/atoms.list")


def canonicalize_res_name(name):
    if name in res_dict and res_dict[name] != "-":
        return res_dict[name]
    else:
        return name


def canonicalize_atom_name(name):
    if name in atom_dict and atom_dict[name] != "-":
        return atom_dict[name]
    else:
        return name


def remove_string_end(s, end):
    assert s.endswith(end)
    return s[: -len(end)]


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
                out_residue = PDB.Residue.Residue(
                    residue.id, canonicalize_res_name(residue.resname), residue.segid
                )
                out_chain.add(out_residue)

                for atom in residue:
                    if is_nuc(residue):
                        can_name = canonicalize_atom_name(atom.name)
                    else:
                        can_name = atom.name
                    full_can_name = can_name
                    out_atom = PDB.Atom.Atom(
                        name=can_name,
                        coord=atom.coord,
                        bfactor=atom.bfactor,
                        occupancy=atom.occupancy,
                        altloc=atom.altloc,
                        fullname=full_can_name,
                        serial_number=atom.serial_number,
                        element=atom.element,
                    )
                    out_residue.add(out_atom)

    return out_struct


def is_nuc(residue):
    if res_dict.get(residue.resname, "-") in "AUCG":
        return True
    else:
        return False


def mc_name_to_tuple(s):
    chain_id = s[:1]
    s = s[1:]
    if "." in s:
        a, b = s.split(".")
        res_id = int(a)
        insertion_code = b
    else:
        res_id = int(s)
        insertion_code = " "

    return (chain_id, res_id, insertion_code)


def tuple_to_mc_name(t):
    (chain_id, res_id, insertion_code) = t
    if insertion_code in (" ", ""):
        return chain_id + str(res_id)
    else:
        return chain_id + str(res_id) + "." + insertion_code


def get_mc_style_name(residue):
    (
        _struct_id,
        _model_id,
        chain_id,
        (_hetero, res_id, insertion_code),
    ) = residue.get_full_id()
    if insertion_code == " ":
        return chain_id + str(res_id)
    else:
        return chain_id + str(res_id) + "." + insertion_code


def mc_name_to_tuple(s):
    chain_id = s[:1]
    s = s[1:]
    if "." in s:
        a, b = s.split(".")
        res_id = int(a)
        insertion_code = b
    else:
        res_id = int(s)
        insertion_code = " "

    return (chain_id, res_id, insertion_code)


def tuple_to_mc_name(t):
    (chain_id, res_id, insertion_code) = t
    if insertion_code in (" ", ""):
        return chain_id + str(res_id)
    else:
        return chain_id + str(res_id) + "." + insertion_code


def query_mc_name(model, name):
    (chain_id, res_id, insertion_code) = mc_name_to_tuple(name)
    return model[chain_id][(" ", res_id, insertion_code)]
