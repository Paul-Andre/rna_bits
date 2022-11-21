set -e  # exit on error

(
cd ~/rnabayespairing2/bayespairing/src
python3 parse_sequences.py \
  -seq "GGGGGGCCCCCC" \
  -ss  "((((....))))"
)

bp2_bridge -r ~/rnabayespairing2/bayespairing/output/output.json -d ~/rnabayespairing2/bayespairing/models/ALL.json --ss "((((....))))"

