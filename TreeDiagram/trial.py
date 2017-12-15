import tree

file = 'seqHydrolase.fasta'
aligned_file = tree.alignment(file)
tree.draw_tree(aligned_file)


