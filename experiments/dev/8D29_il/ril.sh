set -e
rna_insert_loops --seq GGCGAUACCAGCGAAACACGCCCUUGGCAGCGUC --ss "((((...((((((.....)))...)))...))))" --top 10 -o loops_sample_10 --sample 50 --exclude "NR_4.0_56838.1"
rna_insert_loops --seq GGCGAUACCAGCGAAACACGCCCUUGGCAGCGUC --ss "((((...((((((.....)))...)))...))))" --top 20 -o loops_sample_20 --sample 50 --exclude "NR_4.0_56838.1"
rna_insert_loops --seq GGCGAUACCAGCGAAACACGCCCUUGGCAGCGUC --ss "((((...((((((.....)))...)))...))))" --top 38 -o loops_sample_38 --sample 50 --exclude "NR_4.0_56838.1"
rna_insert_loops --seq GGCGAUACCAGCGAAACACGCCCUUGGCAGCGUC --ss "((((...((((((.....)))...)))...))))" --top 50 -o loops_sample_50 --sample 50 --exclude "NR_4.0_56838.1"
rna_builder loops_sample*/*.rass
