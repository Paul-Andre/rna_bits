if g is your fr3d_graph / x3dna_graph
for (strand, position), data in g.nodes(data=True):
  nucleotide = data['nucleotide']
  print(f"Strand {strand} position {position} has a nucleotide {nucleotide.real_nt} ({nucleotide.nt}) at position {nucleotide.pos} ({nucleotide.pdb_pos + nucleotide.pdb_pos_ins})")
