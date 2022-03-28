import vpython as vp
import Bio.PDB as PDB

hairpin = PDB.PDBParser().get_structure("bruh", "RNAMoIP/No_Redondance_VIEW3D/437D.A.1/1L2X.A.1.pdb")[0];
helix_fragment = PDB.PDBParser().get_structure("UGCA", "mcsym-db/2_2/UGCA/ds-2_2-UGCA_x.pdb")[22];
helix_fragment_2 = PDB.PDBParser().get_structure("UUAA", "mcsym-db/2_2/UUAA/ds-2_2-UUAA_x.pdb")[22];

# Reminder: Structure>Model>Chain>Residue>Atom

def residueToCanonicalDict(res):
    return {canonical_atom_name(atom.name):atom for atom in res}

# Given two nucleotides of the same type,
# returns two lists of corresponding atoms
# (references to the existing atoms, not copies)
# Tries to resolve atom naming convention differences
def correspondNucleotideAtoms(a, b):
    a_canon = {canonical_atom_name(atom.name):atom for atom in a}
    b_canon = {canonical_atom_name(atom.name):atom for atom in b}
    assert(len(a_canon) == len(a))
    assert(len(b_canon) == len(b))
    all_atoms = set(a_canon.keys()) | set(b_canon.keys())
    common_atoms = set(a_canon.keys()) & set(b_canon.keys())
    diff = all_atoms-common_atoms
    if (diff):
        print("Some atoms were not shared", diff)
    a_ret = []
    b_ret = []
    for atom_name in common_atoms:
        a_ret.append(a_canon[atom_name])
        b_ret.append(b_canon[atom_name])
    return (a_ret, b_ret)

def canonical_atom_name(name):
    if name == "O1P":
        return "OP1"
    if name == "O2P":
        return "OP2"
    return name.replace("*", "'")


fixed_1,moving_1 = correspondNucleotideAtoms( hairpin["A"][7],  helix_fragment["A"][2])
fixed_2,moving_2 = correspondNucleotideAtoms( hairpin["A"][14],  helix_fragment["A"][3])

superimposer = PDB.Superimposer()
superimposer.set_atoms(fixed_1+fixed_2, moving_1+moving_2)

superimposer.apply(helix_fragment.get_atoms())

fixed_1,moving_1 = correspondNucleotideAtoms(helix_fragment["A"][1], helix_fragment_2["A"][2])
fixed_2,moving_2 = correspondNucleotideAtoms(helix_fragment["A"][4], helix_fragment_2["A"][3])

superimposer = PDB.Superimposer()
superimposer.set_atoms(fixed_1+fixed_2, moving_1+moving_2)

superimposer.apply(helix_fragment_2.get_atoms())

BB_ATOM_NAMES = ["P","O5'","C5'","C4'","C3'","O3'"];


def toVpVec(a):
    return vp.vector(a[0],a[1],a[2])

def drawStructure(struct, curve_color):
    bb_chains =[];
    for chain in struct.get_chains():
        bb_chains.append([])
        for residue_ in chain.get_residues():
            residue = residueToCanonicalDict(residue_)
            for name in BB_ATOM_NAMES:
                bb_chains[-1].append(residue[name])
            for atom in residue_:
                vp.sphere(pos=toVpVec(atom.coord), radius=0.2, color=curve_color)

    for chain in bb_chains:
        for i in range(len(chain) - 1):
            vp.curve(toVpVec(chain[i].coord), toVpVec(chain[i+1].coord), color=curve_color)
            atom = chain[i]
            if atom.name == "P":
                vp.sphere(pos=toVpVec(atom.coord), radius=0.3, color=vp.vector(1,0.25,0))


drawStructure(hairpin, vp.vector(1,1,1))
drawStructure(helix_fragment, vp.vector(0,1,2))
drawStructure(helix_fragment_2, vp.vector(1,0,2))


