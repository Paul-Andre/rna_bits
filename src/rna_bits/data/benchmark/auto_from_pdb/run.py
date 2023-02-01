from rna_bits.utils.data_path import get_path
from rna_bits.data.benchmark.auto_from_pdb.select_simple_strutures import  generate_structure_files
from rna_bits.data.benchmark.auto_from_pdb.insert_loops import  run_rna_insert_loop
from rna_bits.data.benchmark.auto_from_pdb.insert_loops_force_native import  run_insert_loops_force_native
from rna_bits.data.benchmark.auto_from_pdb.build import run_builder
from rna_bits.data.benchmark.auto_from_pdb.evaluate import run_evaluate

# OUT_DIR_ALL = get_path("benchmark/auto_from_pdb/all/provided_ss", create=True)
# generate_structure_files(OUT_DIR_ALL, allow_multi_chain=False, allow_junctions=True)

# OUT_DIR_PB2 = get_path("benchmark/auto_from_pdb/bp2_limited/provided_ss", create=True)
# generate_structure_files(OUT_DIR_PB2, allow_multi_chain=False, allow_junctions=False)
