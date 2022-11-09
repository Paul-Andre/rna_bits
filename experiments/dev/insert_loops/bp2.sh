set -e  # exit on error

cd ~/rnabayespairing2/bayespairing/src
python3 parse_sequences.py \
  -seq "CGUCUUUAUAGCCCAAGGGUAGCCGUAACAAACGCCAAAGCUCCGUAGUAACUGAAAAGAAGAUAACUCAUGAGUAGUCACCCACCCACUAGCCACCAGG" \
  -ss  "..(((((.(((......((.((((((.....))).....))))).......))).)))))......((..((..((((.........))))..))..))."
bp2
