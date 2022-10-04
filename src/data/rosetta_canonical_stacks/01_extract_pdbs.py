import os
import glob
import shutil
from os.path import join as pjoin
import sys
import subprocess


if "ROSETTA3" in os.environ:
    ROSETTA = os.environ["ROSETTA3"]
elif "ROSETTA" in os.environ:
    ROSETTA = os.environ["ROSETTA"]
else:
    print("$ROSETTA3 or $ROSETTA not in environment. Install Rosetta and set $ROSETTA to Rosetta's main/source/ directory.")
    exit(1)
extract_pdb_executables = glob.glob(pjoin(ROSETTA, "bin/extract_pdb*"))
assert(len(extract_pdb_executables) == 1)
EXTRACT_PDB = extract_pdb_executables[0]

SCRIPT_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.abspath(pjoin(SCRIPT_DIR, "../../../data/"))

IN_DIR = pjoin(DATA_DIR, "original/rosetta_rna_db/base_pair_steps/canonical/standard")

INTERIM_DIR = pjoin(DATA_DIR, "interim/rosetta_canonical_stacks")
try:
    os.mkdir(INTERIM_DIR)
except FileExistsError:
    pass

OUT_DIR = pjoin(INTERIM_DIR, "01_pdbs")
try:
    os.mkdir(OUT_DIR)
except FileExistsError:
    shutil.rmtree(OUT_DIR)
    os.mkdir(OUT_DIR)


PAIRS = ("au", "ua", "gc", "cg", "gu", "ug")

def make_rosetta_name(a, b):
    return a[0]+b[0]+"_"+b[1]+a[1]

def make_mcsym_name(a, b):
    return (a[0]+b[0]+b[1]+a[1]).upper()

for a in PAIRS:
    for b in PAIRS:
        n = make_mcsym_name(a,b)
        os.mkdir(pjoin(OUT_DIR, n))
        shutil.copy(pjoin(IN_DIR, make_rosetta_name(a,b) + ".out"), pjoin(OUT_DIR, n, "silent.out"))
        os.chdir(pjoin(OUT_DIR, n))
        # Some of the files in the silent file need to get extracted to a kieft_GU folder.
        # Why? I have no idea. There wasn't any occurence of "kieft" in the Rosetta source code
        os.mkdir("kieft_GU")
        subprocess.run([EXTRACT_PDB, "-in:file:silent", "silent.out"])
        for f in glob.glob("kieft_GU/*.pdb"):
            shutil.move(f, ".")
        os.rmdir("kieft_GU")
        os.remove("silent.out")
