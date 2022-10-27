import os

from Bio.PDB.PDBList import PDBList

from rna_bits.utils.data_path import get_path

# The codes that seem to have been used to generate stacks in Rosetta
pdb_codes = [
    "1jj2",
    "2pxp",
    "2pxl",
    "2pxe",
    "2pxt",
    "2pxf",
    "2pxd",
    "2pxq",
    "2pxb",
    "2pxu",
    "2pxk",
    "2pxv",
]

out_path = get_path("interim/stacks/pdb", create=True)

for pdb_code in pdb_codes:
    # I need to download the pdb version in order to pass it through MC-Annotate
    fn = PDBList(pdb=out_path).retrieve_pdb_file(
        pdb_code, pdir=out_path, file_format="pdb"
    )

    os.rename(fn, os.path.join(os.path.dirname(fn), pdb_code + ".pdb"))
