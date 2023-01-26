Provides the `rna_bits` python library.

Provides the `rna_builder`, `rna_insert_loops` and `bp2_bridge` command line utilities.

<!--
The builder, the loop_inserter and the BayesPairing2 bridge. As well as scripts that generate the data that is used by these 3 utilities



Running from requires downloading McAnnotate ...ads f...

Otherwise, ... running from precomputed ...
-->

# Requirements
Python >= 3.8

Python dependencies will be installed automatically.

Running certain data generation pipelines requires MC-Annotate, which needs to be downloaded separately.
For convenience, I have included the output of MC-Annotate for those pipelines. (TODO: MC-Annotate installation pointer.)

As of now, I have only tested on Linux.

# Installation
Make sure you're inside the directory and run:
```bash
pip3 install -e .
```
It will automatically install the python requirements, and make the python library and command line utilities available.

NOTE: Installing the python library is required in order to run the data generation pipelines.

# rna_builder
Builds a 3d RNA structure from a given .rass file and .pdb files of fragments.

## Example
```
cd examples/rna_builder
rna_builder *.rass
```

# rna_insert_loops
Given an RNA secondary structure and sequence, creates a .rass file with loops 

For it to run, loops need to be first generated.

# bp2_bridge
Used to turn an output from BayesPairing2 (https://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing2) into an input that can be passed to rna_builder.

Requires access to the 

## Example
First, install BayesPairing2 (follow instructions in above link.)

Next, run BayesPairing2 using the `BayesPairing` command (if you installed BayesPairing2 in a new conda environment, make sure that conda environment is active). Example command:
```
mkdir bp2b_example  # create an empty directory
cd bp2b_example
BayesPairing -seq "UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU" -samplesize 1000 -d RELIABLE
```

Next, use bp2_bridge to generate the rass files.

bp2_bridge has 2 required arguments: the output.json from BayesPairing2, and the models json file that BayesPairing2 used. (The `--chefs_choice` flag is optional, to restrict to using only the "chef's choice" motifs as output by BayesPairing2.)
```
bp2_bridge --bp2_result output.json --database ~/where_you_installed_bp2/rnabayespairing2/bayespairing/models/RELIABLE.json --chefs_choice
```
This might take some time because it will download PDB files to generate the 3d motifs.

Finally, assemble the 3d structures:
```
rna_builder bp2b_out/*.rass
```
