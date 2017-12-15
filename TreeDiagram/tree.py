from Bio import AlignIO, SeqIO, Seq, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import pylab
import os

def alignment(input_file):
    """
    File(.fasta) - > File(.fasta)
    ---------------------------------------------
    This function takes in a DNA or Protein
    and performs a multiple sequence alignment
    and returns a the aligned file
    """

    records = SeqIO.parse(input_file, 'fasta')
    records = list(records) # make a copy, otherwise our generator
                            # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen,'.')
            record.seq = Seq.Seq(sequence)
            
    assert all(len(record.seq) == maxlen for record in records)

    # write to temporary file and do all alignment
    output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    alignment = AlignIO.read(output_file, 'fasta')

    return alignment

def get_label(leaf):
    return leaf.name

def draw_tree(alignment):
    """
    File(.fasta) -> None
    ----------------------------------------------------------
    Takes in a fasta file and produces a tree diagram based on it
    """

    
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    newick = input("Enter the name of the file to be saved: ")

    Phylo.write(tree, newick + '.nwk', 'newick')

    t = next(Phylo.parse(newick +'.nwk', 'newick'))
    t.rooted = True
    t.ladderize()

    Phylo.draw_ascii(t)
    Phylo.draw(t, branch_labels=lambda c: c.branch_length, label_func=get_label)
