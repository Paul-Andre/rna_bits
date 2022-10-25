import utils.pdb
from Bio.PDB.PDBList import PDBList
from utils.data_path import get_path
import subprocess
import os

from utils.misc import remove_string_end

# TODO put it in its utility file
MCA_DIR = "/home/paul/MC-Annotate"

pdb_code = "1jj2"

in_path = get_path("interim/stacks/pdb")
out_path = get_path("interim/stacks/mcout", create=True)

filenames = [a for a in os.listdir(in_path) if a.endswith(".pdb")]
filenames.sort()
assert filenames

for i,fn in enumerate(filenames):
    print(fn, str(i+1)+"/"+str(len(filenames)))
    pdb_code = remove_string_end(fn, ".pdb")

    pdb_fn = os.path.join(in_path, pdb_code +".pdb")
    mcout_fn = os.path.join(out_path, pdb_code+".mcout")

    fp = subprocess.run([MCA_DIR, "-f", "0", pdb_fn], capture_output=True, encoding="ascii")
    mcout = fp.stdout
    assert("Invalid" not in mcout)

    with open(mcout_fn, "w") as f:
        f.write(mcout)
