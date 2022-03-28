import logging


class Atom(object):
    def __init__(self, atom_type, atom_label, x, y, z):
        self.atom_type = atom_type
        self.atom_label = atom_label
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"ATOM: {self.atom_type} {self.atom_label} {self.x} {self.y} {self.z}"


class Nucleotide(object):
    def __init__(self, pos, nt, real_nt,
                 chemically_modified, pdb_pos, pdb_pos_ins):
        self.pos = pos
        self.nt = nt
        self.real_nt = real_nt
        self.chemically_modified = chemically_modified
        self.pdb_pos = pdb_pos
        self.pdb_pos_ins = pdb_pos_ins
        self.auth_seq_id = None
        self.atoms = []

    def add_atom(self, atom):
        if not isinstance(atom, Atom):
            logging.debug(
                f"Trying to insert non Atom object in\
                          nucleotide at pos {self.pos}.\nThe atom to be inserted is:\n{self.atom}")
            raise Exception("This is not an Atom object")
        self.atoms.append(atom)

    def __repr__(self):
        return f"NUCLEOTIDE: {self.nt} {self.real_nt} {self.chemically_modified} {self.pdb_pos} {self.pdb_pos_ins} {self.atoms}"


class Strand(dict):

    def __init__(self, *args, **kwargs):  # , name, entity_id, description):
        super().__init__()
        if kwargs:
            self.name = kwargs.get('name', '')
            self.description = kwargs.get('description', '')
            self.entity_id = kwargs.get('entity_id', '')

    def __setitem__(self, key, val):
        if not isinstance(val, Nucleotide):
            logging.debug(
                f"Trying to add a non Nucleotide object to Strand.\n{val}")
            raise Exception(
                f"Trying to add a non Nucleotide object to RNA_Strand.\n{val}")
        if not isinstance(key, int):
            logging.debug(
                f"Trying to add a Nucleotide with non-integer key to strand\n{key}")
            raise Exception(
                f"Trying to add a Nucleotide with non-integer key to strand\n{key}")
        super().__setitem__(key, val)


class RNA_Molecule(dict):

    def __init__(self, pdb_id, title):
        self.pdb_id = pdb_id
        self.title = title
        self.fr3d_graph = None
        super().__init__()

    def __setitem__(self, key, val):
        if not isinstance(val, Strand):
            logging.debug(
                f"Trying to add a non Strand object to RNA_Molecule.\n{val}")
            raise Exception(
                f"Trying to add a non Strand object to RNA_Molecule.\n{val}")
        super().__setitem__(key, val)

    def __repr__(self):
        return f"RNA MOLECULE: {self.pdb_id} {self.title}"

"""
The tuple is of the type (str, int) where the string is the chain, and the integer the position

So if you have your graph, the information is always in the attribute of the nodes / edges . 

For the nodes data you can do something like that :
#if g is your fr3d_graph / x3dna_graph
for (strand, position), data in g.nodes(data=True):
  nucleotide = data['nucleotide']
  print(f"Strand {strand} position {position} has a nucleotide {nucleotide.real_nt} ({nucleotide.nt}) at position {nucleotide.pos} ({nucleotide.pdb_pos + nucleotide.pdb_pos_ins})")
now the differences. 

nucleotide.real_nt vs nucleotide.nt the latter is always in ACGU, the other can be a code indicating chemical modifications. 

nucleotide.pos + nucleotide.pdb_pos + nucleotide.pdb_pos_ins the first one is a rational ordering, the other is the old numbers given by the others in the PDB (composed of an integer pdb_pos and potentially a string of characters pdb_pos_ins)

"""

