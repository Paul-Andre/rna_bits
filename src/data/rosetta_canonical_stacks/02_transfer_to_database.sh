set -e
DATA=$(dirname $BASH_SOURCE)/../../../data

mkdir -p $DATA/database/rosetta_canonical
rm -rf $DATA/database/rosetta_canonical/2_2
cp -r --no-target-directory $DATA/interim/rosetta_canonical_stacks/01_pdbs $DATA/database/rosetta_canonical/2_2
