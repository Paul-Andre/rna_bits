BayesPairing -seq GGCGAUACCAGCGAAACACGCCCUUGGCAGCGUC -ss "((((...((((((.....)))...)))...))))" -d ALL

bp2_bridge \
  --ss "((((...((((((.....)))...)))...))))" \
  -d "/home/paul/Masters/installation_try/rnabayespairing2/bayespairing/models/ALL.json" \
  -m "/home/paul/Masters/RNA/data/benchmark/bp2b_motifs/ALL/" \
  --output_absolute_path \
  -r output.json

