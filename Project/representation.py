from typing import List, Tuple, Union

import Bio.PDB as PDB

class Nucleotide:
    def __init__(self, index: int, pos: Tuple[int,int], letter: str) -> None:
        self.index: int = index
        self.pos: Tuple[int, int] = pos
        self.letter: str = letter

        self.next: Nucleotide = None
        self.prev: Nucleotide = None

        self.pair: Pair = None
        self.paired: Nucleotide = None

        self.motifPositions: List[Tuple[Motif, int]] = None

class Pair:
    """Represents 2 nucleotides paired using cis W/W
    (intended to be used from the secondary structure, mostly for canonical or wobble)
    """
    def __init__(self, a, b):
        self.a: Nucleotide = a
        self.b: Nucleotide = b
        self.nextStack: Stack = None
        self.nextStacked: Pair = None
        self.prevStack: Stack = None
        self.prevStacked: Pair = None

class Stack:
    """ Represents a pair of stacked base pairs
    Intended to create a helix.
    """
    def __init__(self, a, b):
        self.a: Pair = a
        self.b: Pair = b
        self.next: Stack = None
        self.prev: Stack = None

class FreeLink:
    """ Represents that two nucleotides are connected by the backbone and have
    no other constraint (Stack or Motif) imposed on them
    """
    def __init__(self, a, b):
        self.a: Nucleotide = a
        self.b: Nucleotide = b

class Motif:
    """ Represents that we want to insert a motif at a certain position in the structure
    Contains the 3d model(s) we will try to insert
    """
    def __init__(self, nucleotides: List[Nucleotide], models) -> None:
        self.nucleotides: List[Nucleotide] = nucleotides
        # TODO: instead of a list, perhaps make this an abstract object to
        # allow other ways of sampling models
        self.models: List[PDB.Model] = models
