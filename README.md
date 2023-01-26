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

# rna_insert_loops
Given an RNA secondary structure and sequence, creates a .rass file with loops 

For it to run, loops need to be first generated.

# bp2_bridge
Used to turn an output from BayesPairing2 (https://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing2) into an input that can be passed to rna_builder. 
