set -e
DATA=$(dirname $BASH_SOURCE)/../../../data

ROSETTA_DB=${ROSETTA3_DB:-$ROSETTA_DB}
if [[ -z $ROSETTA_DB ]]; then
  echo "\$ROSETTA3_DB or \$ROSETTA_DB not found. Download Rosetta and either set the environment variable and rerun this script, or manually copy the contents of Rosetta's main/database/sampling/rna/ into data/original/rosetta-db/"
  exit 1
fi

cp -r --no-target-directory $ROSETTA_DB/sampling/rna/ $DATA/original/rosetta_rna_db
echo $ROSETTA_DB/sampling/rna/ > $DATA/original/rosetta_rna_db/source.txt
