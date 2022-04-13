import sys;
import numpy as np

f = sys.stdin

# internally, I use 0-based indexing for now

seq = f.readline().strip()
dot_bracket = f.readline().strip()
assert(len(seq) == len(dot_bracket))
module_components = []
while True:
    line = f.readline()
    if not line:
        break;
    line = line.strip()
    info = line.split("-")
    if (info[0] != 'C'):
        continue;
    _,name,first,last,comp = info
    first,last,comp = map(int, (first,last,comp))
    first-=1
    last-=1
    comp-=1
    module_components.append((name, comp, first, last))
f.close()

pairing = [None]*len(seq)

pos = 0
def find_pairing():
    global pos
    while(pos<len(dot_bracket)):
        if dot_bracket[pos]=='.':
            pos+=1
            continue;
        elif dot_bracket[pos]==')':
            return
        elif dot_bracket[pos]=='(':
            open_pos = pos
            pos+=1
            find_pairing()
            assert(pos<len(dot_bracket) and dot_bracket[pos]==')')
            pairing[open_pos]=pos
            pairing[pos]=open_pos
            pos+=1
        else:
            assert False, "Unexpected character in dot-bracket string: "+repr(dot_bracket[pos])


import Project.parser
for a,b in Project.parser.parseParens(dot_bracket):
    a-=1
    b-=1
    pairing[a]=b
    pairing[b]=a


#find_pairing()

# We build a tree where the nodes are NCMs or modules, and the edges are overlaps
# Currently NCMs we look at are 2_2, 2_3 (and 3_2) and 3_3

# TODO: deal with the same module coming up in two places

class Node:
    id_increment = 10000
    def __init__(self, kind):
        self.kind = kind
        self.components = []
        self.edges = []
        self.nucls = set()
        self.id = Node.id_increment
        Node.id_increment+=1
        self.visited = False
        self.model = None
        self.user_added = (kind[0] == "module")

module_assignment = [[] for _ in range(len(seq))]

seen_comp = set()
modules_by_name = {}

nodes = []

module_components.sort() # just to make sure in order
for name,comp,first,last in module_components:
    assert (not (name,comp) in seen_comp), "Currently we only deal with one instance of a module"
    if name in modules_by_name:
        module = modules_by_name[name]
    else:
        module = Node(("module", name))
        nodes.append(module)
        modules_by_name[name] = module
    assert(len(module.components) == comp), "Index of module component does not match order"
    module.components.append((first,last))
    for i in range(first, last+1):
        module_assignment[i].append((module, comp, i-first))
        module.nucls.add(i)

added_ncms = set()

def have_overlapping_module(nucls):
    b = [
            {node.id for (node,_,_) in module_assignment[id]} for id in nucls
    ]
    a = set.intersection(*b)
    return (len(a) >= 1)

for i in range(0, len(seq)):
    if pairing[i] is not None:
        j = pairing[i]
        ii = i+1
        jj = j-1
        iii = i+2
        jjj = j-2
        module = None
        if ii < len(pairing) and pairing[ii] == jj:
            ts = tuple(sorted((i,ii,jj,j)))
            if (ts in added_ncms or have_overlapping_module(ts)):
                continue
            added_ncms.add(ts)
            module = Node(("ncm", "2_2", seq[i]+seq[ii]+seq[jj]+seq[j]))
            module_assignment[i].append((module,0,0))
            module_assignment[ii].append((module,0,1))
            module_assignment[jj].append((module,1,0))
            module_assignment[j].append((module,1,1))
            module.nucls.update({i,ii,jj,j})
            module.components = [(i,ii),(jj,j)]

        elif ii < len(pairing) and pairing[ii] == jjj:
            ts = tuple(sorted((i,ii,jjj,jj,j)))
            if (ts in added_ncms or have_overlapping_module(ts)):
                continue
            added_ncms.add(ts)
            module = Node(("ncm", "2_3", seq[i]+seq[ii]+seq[jjj]+seq[jj]+seq[j]))
            module_assignment[i].append((module,0,0))
            module_assignment[ii].append((module,0,1))
            module_assignment[jjj].append((module,1,0))
            module_assignment[jj].append((module,1,1))
            module_assignment[j].append((module,1,2))
            module.nucls.update({i,ii,jjj,jj,j})
            module.components = [(i,ii),(jjj,j)]

        elif iii < len(pairing) and  pairing[iii] == jj:
            ts = tuple(sorted((i,iii,iii,jj,j)))
            if (ts in added_ncms or have_overlapping_module(ts)):
                continue
            added_ncms.add(ts)
            module = Node(("ncm", "3_2", seq[i]+seq[ii]+seq[iii]+seq[jj]+seq[j]))
            module_assignment[i].append((module,0,0))
            module_assignment[ii].append((module,0,1))
            module_assignment[iii].append((module,0,2))
            module_assignment[jj].append((module,1,0))
            module_assignment[j].append((module,1,1))
            module.nucls.update({i,ii,iii,jj,j})
            module.components = [(i,iii),(jj,j)]

        elif iii < len(pairing) and  pairing[iii] == jjj:
            ts = tuple(sorted((i,ii,iii,jjj,jj,j)))
            if (ts in added_ncms or have_overlapping_module(ts)):
                continue
            added_ncms.add(ts)
            module = Node(("ncm", "3_3", seq[i]+seq[ii]+seq[iii]+seq[jjj]+seq[jj]+seq[j]))
            module_assignment[i].append((module,0,0))
            module_assignment[ii].append((module,0,1))
            module_assignment[iii].append((module,0,2))
            module_assignment[jjj].append((module,1,0))
            module_assignment[jj].append((module,1,1))
            module_assignment[j].append((module,1,2))
            module.nucls.update({i,ii,iii,jjj,jj,j})
            module.components = [(i,iii),(jjj,j)]

        if module:
            nodes.append(module)

# Add isolated pairs
for i in range(0, len(seq)):
    if pairing[i] is not None:
        j = pairing[i]
        ts = tuple(sorted((i,j)))
        if (ts in added_ncms or have_overlapping_module(ts)):
            continue
        added_ncms.add(ts)
        module = Node(("pair", seq[i]+seq[j]))
        module_assignment[i].append((module,0,0))
        module_assignment[j].append((module,1,0))
        module.nucls.update({i,j})
        module.components = [(i,i), (j,j)]
        nodes.append(module)


# Add single nucleotides
for i in range(0, len(seq)):
    if len(module_assignment[i]) == 0:
        ts = (i,)
        if (ts in added_ncms or have_overlapping_module(ts)):
            continue
        added_ncms.add(ts)
        module = Node(("nucleotide", seq[i]))
        module_assignment[i].append((module,0,0))
        module.nucls.update({i})
        module.components = [(i,i)]
        nodes.append(module)



class Edge:
    def __init__(self,a,b):
        self.a=a
        self.b=b
        self.nucls=set()

edges_by_nodes = {}

for i in range(0,len(seq)):
    ass = module_assignment[i]
    for mod, comp, pos in ass:
        print(mod.kind, mod.components)
    # I comment this line out to start implementing arcs:
    # assert len(ass)>=1, "nucleotide "+repr(i)+" isn't assigned to a node"
    assert len(ass)<=2, "nucleotide "+repr(i)+" is assigned to more than 2 components"
    if len(ass) == 2:
        a = ass[0][0]
        b = ass[1][0]
        ts = tuple(sorted((a.id,b.id)))
        if ts in edges_by_nodes:
            edge = edges_by_nodes[ts]
        else:
            edge = Edge(a,b)
            a.edges.append(edge)
            b.edges.append(edge)
            edges_by_nodes[ts] = edge
        edge.nucls.add(i)


for edge in edges_by_nodes.values():
    #assert len(edge.nucls) == 2, "edge has "+repr(len(edge.nucls))+" nucleotides"
    assert len(edge.nucls) in (1,2), "edge has "+repr(len(edge.nucls))+" nucleotides"


import vpython as vp
import Bio.PDB as PDB

import os
import random

import gzip


# TODO: take the correspondences from the RNA-Puzzles assessment source
def canonical_atom_name(name):
    if name == "O1P":
        return "OP1"
    if name == "O2P":
        return "OP2"
    return name.replace("*", "'")

def canonical_residue_name(name):
    # TODO
    # Note: in this function we won't rename modified nucleotide names
    return name

def canonical_unmodified_residue_name(name):
    # TODO
    # Note: in this function we will rename modified nucleotide names
    return name

# Creates a new model that has canonical atom representation
def canonicalize_model(model):
    output_model = PDB.Model.Model(0)
    for chain in model:
        output_chain = PDB.Chain.Chain(chain.id)
        output_model.add(output_chain)
        for residue in chain:
            output_residue = PDB.Residue.Residue(residue.id, canonical_residue_name(residue.resname), residue.segid)
            output_chain.add(output_residue)
            for atom in residue:
                can_name = canonical_atom_name(atom.name)
                full_can_name = can_name
                output_atom = PDB.Atom.Atom(
                        name = can_name,
                        coord = atom.coord,
                        bfactor = atom.bfactor,
                        occupancy = atom.occupancy,
                        altloc = atom.altloc,
                        fullname = full_can_name,
                        serial_number = atom.serial_number,
                        element = atom.element,
                        )
                output_residue.add(output_atom)

    return  output_model

def is_struct_filename(filename):
    return (
            filename.endswith(".pdb") or
            filename.endswith(".pdb.gz") or
            filename.endswith(".cif") or
            filename.endswith(".cif.gz")
            )

def load_struct_file(filename, struct_name=None):
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

    assert(False), "Not a proper filename"

# returns (model, filename)
# TODO: put the filename selection into a different file
# TODO: add a caching layer
def load_model_kind(kind):
    if kind[0] == "module":
        name = kind[1]
        directory = name
        if (os.path.isdir(directory)) :
            options = [os.path.join(directory, x) for x in os.listdir(directory) if is_struct_filename(x) and
                    not x.startswith(".")]
            choice = random.choice(options)
        else:
            choice = directory
        model = load_struct_file(choice)[0]
        model = canonicalize_model(model)
        return model,choice
    if kind[0] == "ncm":
        directory = "mcsym-db/"+kind[1]+"/"+kind[2]
        # if the theoretical model exists:
        theoretical_filename = os.path.join(directory, "1-"+kind[1]+"-"+kind[2] +"_t.pdb.gz")
        if os.path.isfile(theoretical_filename):
            filename = theoretical_filename
        else:
            options = [os.path.join(directory, x) for x in os.listdir(directory)
                       if ( x.endswith('.pdb.gz') and not x.startswith("."))]
            filename = random.choice(options)
        model = load_struct_file(filename)[0]
        model = canonicalize_model(model)
        return model,filename
    if kind[0] == "pair":
        # load 2_2 from mcsym-db, and truncate
        model, filename = load_model_kind(("ncm", "2_2", kind[1][0]+"GC"+kind[1][1]))
        i = 0
        for chain in model:
            residue_ids_to_delete = []
            for residue in chain:
                if i in (1,2):
                    residue_ids_to_delete.append(residue.id)
                i+=1
            for id in residue_ids_to_delete:
                chain.detach_child(id)

        model = canonicalize_model(model)
        return model, (filename, "trimmed")
    if kind[0] == "nucleotide":
        name = kind[1]
        directory = name+".pdb"
        model = load_struct_file(directory)[0]
        model = canonicalize_model(model)
        return model,directory



MAX_COMP_GAP = 4;

def get_model_components(model, kind):
    if (kind[0] == "module"):
        components = []
        for chain in model:
            prev = None
            components.append([])
            for residue in chain:
                if prev is not None:
                    diff = residue.id[1]-prev.id[1]
                    if diff > MAX_COMP_GAP+1:
                        components.append([])
                    elif diff > 1:
                        for _ in range(diff-1):
                            components[-1].append(None)

                components[-1].append(residue)
                prev = residue
        return components

    # TODO refactor this stuff
    elif (kind[0] == "ncm"):
        components = []
        lengths = list(map(int,kind[1].split("_")))
        components = [[] for _ in lengths]
        li = 0
        lj = 0
        for chain in model:
            for residue in chain:
                components[li].append(residue)
                lj+=1
                if(lj>=lengths[li]):
                    lj=0
                    li+=1
        return components

    elif (kind[0] == "pair"):
        components = []
        lengths = [1,1]
        components = [[] for _ in lengths]
        li = 0
        lj = 0
        for chain in model:
            for residue in chain:
                components[li].append(residue)
                lj+=1
                if(lj>=lengths[li]):
                    lj=0
                    li+=1
        return components

    elif (kind[0] == "nucleotide"):
        components = []
        lengths = [1]
        components = [[] for _ in lengths]
        li = 0
        lj = 0
        for chain in model:
            for residue in chain:
                components[li].append(residue)
                lj+=1
                if(lj>=lengths[li]):
                    lj=0
                    li+=1
        return components





# returns a map nucleotide_id -> residue_in_the_model
def map_node_components(node):
    model = node.model

    components = get_model_components(model, node.kind)

    print(node.kind)
    print(node.components, components)
    assert(len(node.components) == len(components)), "number of fragment components doesn't match"
    ret = {}
    for nc, c in zip(node.components, components):
        assert(nc[1]-nc[0]+1 == len(c)), "fragment doesn't match " + str(node.kind)
        for i,r in zip(range(nc[0], nc[1]+1), c):
            if r is not None:
                # r can be None in the case of RNAMoIP-style input where skips of up to 4 nts are allowed
                # TODO: make this coherent with module_assignment
                ret[i] = r
    return ret

models_used = []

# Model contains only the 4 bases, used for base substitution
substitution_model = PDB.PDBParser().get_structure("substitution_model", "bases.pdb")[0]

# Returns (a reference to) the reference base of the given letter, used for base substitution
def get_reference_base(letter):
    return substitution_model[letter][1]


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


def substitute_base(residue, newResname):
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
    superimposer.set_atoms(fixed = fixed_list, moving = moving_list)

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



for node in nodes:
    print(node.kind, node.nucls)
    node.model, model_filename = load_model_kind(node.kind)
    node.model_filename = model_filename
    models_used.append((node.components, model_filename))
    node.component_mapping = map_node_components(node)
    print(node.component_mapping)

    did_substitute = False
    for k,v in node.component_mapping.items():
        if seq[k].isupper() and seq[k] != canonical_residue_name(v.resname):
            did_substitute = True
            substitute_base(v, seq[k])
    #exit()


# Before turning data into the BIO classes, I will represent the chain as a dict of nucleotides first
chain_of_nucleotides = {}



def residueToCanonicalDict(res):
    return {canonical_atom_name(atom.name):atom for atom in res}

BB_ATOM_NAMES = ["P","O5'","C5'","C4'","C3'","O3'"];


# Given two nucleotides of the same type,
# returns two lists of corresponding atoms
# (references to the existing atoms, not copies)
def correspondNucleotideAtoms(a, b, only_bb=False):
    a_canon = {atom.name:atom for atom in a}
    b_canon = {atom.name:atom for atom in b}
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
        if (not only_bb) or atom_name in BB_ATOM_NAMES:
            a_ret.append(a_canon[atom_name])
            b_ret.append(b_canon[atom_name])
    return (a_ret, b_ret)

def superimpose_nodes(fixed, moving, nucls):
    fixed_list = []
    moving_list = []
    for n in nucls:
        fixed_res = fixed.component_mapping[n]
        moving_res = moving.component_mapping[n]
        fixed_atoms,moving_atoms = correspondNucleotideAtoms(fixed_res, moving_res, False)
        fixed_list+=fixed_atoms
        moving_list+=moving_atoms

    superimposer = PDB.Superimposer()
    superimposer.set_atoms(fixed_list, moving_list)

    superimposer.apply(moving.model.get_atoms())

class Nucleotide:
    def __init__(self, model, rigid, pos):
        self.model = model
        self.rigid = rigid
        self.pos = pos
        rigid.nucleotides.add(self)
    def __repr__(self):
        return "Nuc{}".format(self.pos + 1)
    def __lt__(self, other):
        return self.pos < other.pos

class Rigid:
    id_increment = 1000
    def __init__(self):
        # self.translation = nunmpy.array([0.,0.,0.])
        # self.rotation = nunmpy.array([0.,0.,0.])
        self.nucleotides = set()
        self.parent = None
        self.id = Rigid.id_increment
        Rigid.id_increment+=1
    
    # Merge other into self
    def merge(self, other):
        if (self is other):
            return
        for n in other.nucleotides:
            n.rigid = self
        self.nucleotides.update(other.nucleotides)
        other.nucleotides = None
        other.parent = self

def traverse_and_stack(node, currentRigid):
    if node.visited:
        return
    print("Putting down", node.model_filename)
    node.visited=True
    for (k,v) in node.component_mapping.items():
        if k not in chain_of_nucleotides:
            print(k,v)
            wrapper = Nucleotide(v.copy(), currentRigid, k)
            chain_of_nucleotides[k] = wrapper

    for edge in node.edges:
        other = edge.a if edge.b is node else edge.b
        if other.visited:
            continue
        superimpose_nodes(node,other,edge.nucls)
        traverse_and_stack(other, currentRigid)

rigidCounter = 0

rigids = []

while(True):
    central_node = None
    for node in nodes:
        if not node.visited and (central_node is None or (
                len(node.components)>len(central_node.components))):
            central_node = node
    if central_node is None:
        break

    # center the first model, just because
    acc = sum(atom.coord for atom in central_node.model.get_atoms());
    tot = sum(1 for atom in central_node.model.get_atoms());
    centroid = acc/tot;
    rigidCounter+=1

    currentRigid = Rigid()
    rigids.append(currentRigid)

    for atom in central_node.model.get_atoms():
        atom.coord-=centroid + rigidCounter*20.
    traverse_and_stack(central_node, currentRigid)

def createC3PrimeNuc(coord, letter, pos):
    output_residue = PDB.Residue.Residue((" ",0," "), letter, "    ")
    # 'name', 'coord', 'bfactor', 'occupancy', 'altloc', 'fullname', 'serial_number', element=None
    output_atom = PDB.Atom.Atom(
        "C3'", coord, 0., 1.0, " ", "C3'", 1, "C")
    output_residue.add(output_atom)
    rigid = Rigid()
    rigids.append(rigid)
    nucWrapper = Nucleotide(output_residue, rigid, pos)
    return nucWrapper

# Load models for A,U,G,C nucleotides

# Find nucleotides that haven't been placed and turn them into rigids
for i in range(0, len(seq)):
    #assert i in chain_of_nucleotides
    if i not in chain_of_nucleotides:
        print(i, "is not in chain of nucleotides for some reason")

        rigidCounter+=1
        coord = np.array([0.,0,0]) - rigidCounter*20.
        chain_of_nucleotides[i] = createC3PrimeNuc(coord, seq[i], i)

# Find the shortest cycle:
# shortest_cycle = None
# shortest_cycle_length = math.inf


def isRigidEdgeNuc(i):
    if (i+1) in chain_of_nucleotides:
        if (chain_of_nucleotides[i+1].rigid is not chain_of_nucleotides[i].rigid):
            return True
    if (i-1) in chain_of_nucleotides:
        if (chain_of_nucleotides[i-1].rigid is not chain_of_nucleotides[i].rigid):
            return True
    return False


def cycleNotDegenerate(stack):
    return len(stack)>2 or stack[0][0] != stack[0][1] or stack[1][0] != stack[1][1]

def checkCycleNotDegenerate(prevMap, startNuc):
    assert startNuc in prevMap
    c, d = prevMap[startNuc]
    a, b = prevMap[c]
    return not ((c is d) and (a is b) and a is startNuc)

# Distance between C3' atoms.
# Took from NAST
NUCLEOTIDE_DISTANCE = 5.78

# Some arbitrary distance
#NUCLEOTIDE_DISTANCE = 7.5

import heapq

        # if want_chain and isTerminalRigid(nuc.rigid) and nuc.rigid is not startNuc.rigid:
        #     prevMap[nuc] = prev
        #     endNuc = nuc
        #     score = dist
        #     break;

def traverseToGetCycle(startNuc, endNuc):
    # The idea is to get the shortest path from startNuc to startNuc,
    # but without passing the startNuc->endNuc link in that direction (but
    # allowing the endNuc->startNuc direction), so that we find the shortest
    # cycle containing that link.

    prevMap = {}

    queue = [(0., startNuc, None)]

    score = None
    while len(queue):
        dist, nuc, prev = heapq.heappop(queue)

        if nuc in prevMap:
            continue

        if prev is not None:
            prevMap[nuc] = prev

        #print(dist, nuc, prev)

        if nuc is startNuc and prev is not None and checkCycleNotDegenerate(prevMap, startNuc):
            prevMap[nuc] = prev
            score = dist
            break;



        for nn in nuc.rigid.nucleotides:
            dist2 = dist + np.linalg.norm(getC3Prime(nn.model) - getC3Prime(nuc.model))
            for j in (nn.pos+1, nn.pos-1):
                if j not in chain_of_nucleotides:
                    continue
                nnn = chain_of_nucleotides[j]
                if nnn.rigid is nn.rigid:
                    continue

                if nn is startNuc and nnn is endNuc:
                    continue


                dist3 = dist2 + NUCLEOTIDE_DISTANCE

                if nnn not in prevMap:
                    #print("nnnnn",nuc, nn, nnn)
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
# TODO: make this return cycles in the order of their lengths
def getNextCycle():
    print (chain_of_nucleotides)
    bestScore = math.inf
    bestCycle = None
    for i in range(len(seq)):
        if isRigidEdgeNuc(i):
            nuc = chain_of_nucleotides[i]

            for j in (i-1, i+1):
                if j not in chain_of_nucleotides:
                    continue
                nn = chain_of_nucleotides[j]
                if nn.rigid is nuc.rigid:
                    continue

                cycle, score = traverseToGetCycle(nuc, nn)
                if score<bestScore:
                    bestScore = score
                    bestCycle = cycle
    return bestCycle


def isTerminalRigid(r):
    # TODO: make this work even when there's more than 1 strand
    return r in (
            chain_of_nucleotides[len(seq)-1].rigid,
            chain_of_nucleotides[0].rigid
            )

def traverseToGetPath(startNuc):
    """ Traverse to get a path between two rigids that are connected to only one edge each
    (to get an "outer loop")
    """
    assert(isTerminalRigid(startNuc.rigid))
    assert(isRigidEdgeNuc(startNuc.pos))

    visitedRigids = set()
    out = []
    nuc = startNuc
    visitedRigids.add(nuc.rigid)
    while(nuc is startNuc or not isTerminalRigid(nuc.rigid)):
        print(nuc)
        print(out)

        result = None

        for nn in nuc.rigid.nucleotides:
            #dist2 = dist + np.linalg.norm(getC3Prime(nn.model) - getC3Prime(nuc.model))
            print(nuc, nn)
            for j in (nn.pos+1, nn.pos-1):
                if j not in chain_of_nucleotides:
                    continue
                nnn = chain_of_nucleotides[j]
                print(nuc, nn, nnn)
                if nnn.rigid is nn.rigid:
                    continue
                #dist3 = dist2 + NUCLEOTIDE_DISTANCE
                if nnn.rigid not in visitedRigids:
                    result = (nuc, nn, nnn)
                    break;
            if result is not None:
                break;

        nuc, nn, nnn = result
        out.append((nuc, nn))
        nuc = nnn
        visitedRigids.add(nuc)
    out.append((nuc, nuc))
    return out

def getNextPath():
    """ Assumes that all cycles have been found.
    Finds the next "outer loop"
    """
    for i in range(len(seq)):
        if isRigidEdgeNuc(i):
            nuc = chain_of_nucleotides[i]
            return traverseToGetPath(nuc)

    return None


def hasPrev(nuc_id):
    return (nuc_id >= 1) 

def hasNext(nuc_id):
    return (nuc_id < len(seq)-1)

def isNucValid(nuc_id):
    return nuc_id>=0 and nuc_id<len(seq)

import inscribed_polygon

def getC3Prime(residue):
    if "C3'" in residue:
        return residue["C3'"]
    if "C3*" in residue:
        return residue["C3*"]
    return None

def np_distance(a_coord,b_coord):
    return np.linalg.norm(a_coord - b_coord)

import math
import scipy

def ortho_project(s, d):
    return (s.dot(d)/(np.linalg.norm(d)**2)) * d

def create_rotation(s1, s2, d1, d2):
    """ Create a rotation matrix that transforms the vector s1 to align with the vector d1
    while having the vector s2 be aligned with the vector d2 as much as possible """

    # Construct two orthonormal bases
    s1 = s1/np.linalg.norm(s1)
    d1 = d1/np.linalg.norm(d1)

    # print("s1 s2", s1, s2)
    # print("d1 d2", d1, d2)
    s2 = s2 - ortho_project(s2, s1)
    d2 = d2 - ortho_project(d2, d1)

    # print("s1 s2", s1, s2)
    # print("d1 d2", d1, d2)
    s2 = s2/np.linalg.norm(s2)
    d2 = d2/np.linalg.norm(d2)

    s3 = np.cross(s1, s2)
    d3 = np.cross(d1, d2)

    S = np.array([s1,s2,s3])
    D = np.array([d1,d2,d3])

    # Solve SX = D
    T = np.linalg.solve(S, D)
    T = T.transpose()

    # print(s1, d1, T.dot(s1))
    return T


while True:
    cycle = getNextCycle()

    if cycle is not None:
        print("Cycle", cycle)
        isPath = False
    else:
        # for i,nuc in sorted(chain_of_nucleotides.items()):
            # print(i+1, nuc.rigid.id)
        cycle = getNextPath()
        print("Chain", cycle)
        isPath = True

    if cycle is None:
        break;

    rigid_distances = []
    for (a,b) in cycle:
        a_atom = getC3Prime(a.model)
        b_atom = getC3Prime(b.model)
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
    is_degenerate = not isPath and (sum(distances) <= 2*max(distances))
    if is_degenerate:
        # the polygon is going to be degenerate or impossible
        # might need to stretch the distance between nucleotides
        assert( max(distances) != NUCLEOTIDE_DISTANCE )
        adj_nuc_dist = NUCLEOTIDE_DISTANCE + (2*max(distances) - sum(distances))/len(cycle)
        print("stretch distance to", adj_nuc_dist)
        pos = 0.
        points = []
        for d in rigid_distances:
            points.append((pos, 0.))
            if d == max(distances):
                pos -= d
            else:
                pos += d

            points.append((pos, 0.))

            pos += adj_nuc_dist
            print("bridge", points)
    else:
        if not isPath:
            points = inscribed_polygon.construct_polygon(distances)
        else:
            points = inscribed_polygon.construct_polygon(distances, math.tau*0.33)
            points.append(points[-1])

    points = [np.array([x,y,0]) for (x,y) in points]
    # print(points)

    #place the pieces onto the polygon
    assert(len(points)%2 == 0)
    for i in range(len(points)//2):
        dpa = points[i*2]
        dpb = points[i*2+1]
        ra = cycle[i][0]
        rb = cycle[i][1]
        assert(ra.rigid is rb.rigid)

        spa = getC3Prime(ra.model).coord
        spb = getC3Prime(rb.model).coord

        if (spb == spa).all():
            sv1 = np.array([1.,0,0])
            if "P" in ra.model:
                if ( (i == 0 and cycle[i+1][0].pos > cycle[i][1].pos)
                        or cycle[i][0].pos > cycle[i-1][1].pos ):
                    sv1 = spa - ra.model["P"].coord
                else:
                    sv1 = ra.model["P"].coord - spa
        else:
            sv1 = spb-spa

        if is_degenerate:
            if rigid_distances[i] != max(rigid_distances):
                dv1 = np.array([1, 0, 0])
                dv2 = np.array([0, 1, 0])
            else:
                dv1 = np.array([-1, 0, 0])
                dv2 = np.array([0, -1, 0])

        else :
            if (spb == spa).all():
                dv1 = np.array([-dpa[1], dpb[0], 0])
            else:
                dv1 = dpb-dpa

            dv2 = (dpa+dpb)/2
            # For the case of bridges, inverse the side that goes in the other direction
            # To check if the side goes in the other direction, compute CCW
            if (dpa[0]*dpb[1] - dpa[1]*dpb[0])<0:
                dv2*=-1
        
        # Find the centroid of the rigid to orient it away from the loop
        acc = sum(sum(atom.coord for atom in nuc.model) for nuc in ra.rigid.nucleotides);
        num = sum(sum(1 for atom in nuc.model) for nuc in ra.rigid.nucleotides);
        centroid = acc/num;

        # print("centroid", centroid, spa)

        # if it is colinear, switch it to something else
        sv2 = centroid - spa
        if np.linalg.norm(sv2) <= 1e-10 or abs(abs(sv1.dot(sv2))/np.linalg.norm(sv1)/np.linalg.norm(sv2) - 1) <=1e-10:
            sv2 = np.array([1,0,0]) 
        if abs(abs(sv1.dot(sv2))/np.linalg.norm(sv1)/np.linalg.norm(sv2) - 1) <=1e-10:
            sv2 = np.array([0,1,0]) 


        rotation_matrix = create_rotation(sv1, sv2, dv1, dv2)
        # print(rotation_matrix)
        
        for nuc in ra.rigid.nucleotides:
            for atom in nuc.model:
                coord = atom.coord
                atom.coord = rotation_matrix.dot( (coord - spa) ) + dpa

    #combine all rigid bodies
    first = cycle[0][0]
    for (a,b) in cycle[1:]:
        assert(a.rigid is not first)
        assert(a.rigid is b.rigid)
        first.rigid.merge(a.rigid)




def stretch_chain(chain):
    if not (hasPrev(chain[0]) and hasNext(chain[-1])):
        return
    start_id = chain[0]-1
    end_id = chain[-1]+1
    start_res = residueToCanonicalDict(chain_of_nucleotides[start_id].model)
    end_res = residueToCanonicalDict(chain_of_nucleotides[end_id].model)
    start_coord = start_res["C3'"].coord
    end_coord = end_res["C3'"].coord
    diff_vec = end_coord-start_coord
    print("Length", np.linalg.norm(((1)/(1+len(chain)))*diff_vec))
    for i in range(0, len(chain)):
        j = chain[i]
        coord = start_coord + ((1+i)/(1+len(chain)))*diff_vec
        chain_of_nucleotides[j] = createC3PrimeNuc(coord, seq[j])
    #exit()

rigids_seen = set()
for nuc in chain_of_nucleotides.values():
    rigids_seen.add(nuc.rigid)

print("rigids seen", [r.id for r in rigids_seen])

if (False):
    # visualize stuff
    import vpython as vp

    def toVpVec(a):
        return vp.vector(a[0],a[1],a[2])

    for k,residue in chain_of_nucleotides.items():
        for atom in residue:
            vp.sphere(pos=toVpVec(atom.coord), radius=0.3, color=vp.vector(0,0.5,0))
            if canonical_atom_name(atom.name) == "P":
                vp.sphere(pos=toVpVec(atom.coord), radius=0.5, color=vp.vector(1,0.25,0))

# Now all we need is to write a pdb output!
print("Great! now to write it down!")
output_structure = PDB.Structure.Structure("output")
output_model = PDB.Model.Model(0)
output_chain = PDB.Chain.Chain("A")

atom_serial_inc = 1
for (i,rw) in sorted(chain_of_nucleotides.items()):
    r = rw.model
    output_residue = PDB.Residue.Residue((" ",i+1," "), r.resname, "    ")
    for atom in r:
        # 'name', 'coord', 'bfactor', 'occupancy', 'altloc', 'fullname', and 'serial_number'
        output_atom = PDB.Atom.Atom(
            atom.name, atom.coord, atom.bfactor, atom.occupancy, atom.altloc, atom.fullname, atom_serial_inc, atom.element)
        output_residue.add(output_atom)
    output_chain.add(output_residue)

output_model.add(output_chain)
output_structure.add(output_model)

io = PDB.PDBIO()
io.set_structure(output_structure)


if (len (sys.argv) >= 2):
    out = sys.argv[1]
else:
    out = "out.pdb"
io.save(out)
for a in (models_used):
    print(a)
