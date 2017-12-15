from Bio import SeqIO
from CircularGenome import CircularGenome
from Gene import Gene
from collections import namedtuple


def paran(my_list):
    ret = ""
    for i in range(len(my_list)):
        if my_list[i] != '[' or my_list[i] != ']':
            ret += my_list[i]
    return ret


record = SeqIO.read(open("AF323668.gb", "r"), "genbank")
genome = CircularGenome()


def read_file():
    locus_tag = product = protein_id = translation = gen = "N/A"
    start = stop = codon_start = table = strand = 0
    Info = namedtuple('Info', 'locus, gene, protein_id, product, length')
    gene = Gene()
    for feature in record.features:
        if "CDS" in feature.type:
            for value in feature.qualifiers:
                if 'protein_id' in feature.qualifiers.keys():
                    protein_id = paran(feature.qualifiers['protein_id'])
                if 'translation' in feature.qualifiers.keys():
                    translation = paran(feature.qualifiers['translation'])
                if 'product' in feature.qualifiers.keys():
                    product = paran(feature.qualifiers['product'])
                if 'locus_tag' in feature.qualifiers.keys():
                    locus_tag = paran(feature.qualifiers['locus_tag'])
                if 'transl_table' in feature.qualifiers.keys():
                    table = paran(feature.qualifiers['transl_table'])
                if 'codon_start' in feature.qualifiers.keys():
                    codon_start = paran(feature.qualifiers['codon_start'])
                if 'gene' in feature.qualifiers.keys():
                    gen = paran(feature.qualifiers['gene'])
            strand = feature.strand
            start = feature.location.start.position
            stop = feature.location.end.position
            gene = Gene(locus_tag, gen, start, stop, codon_start, table, product, protein_id, translation, strand)  # noqa
            data = Info(locus=locus_tag, gene=gen, protein_id=protein_id, product=product, length=stop-start)  # noqa
            genome.add_node(data, gene)  # noqa


def get_gene():
    gene_set = genome.giveNode(4, "integrase")
    if(isinstance(gene_set, str)):
        print(gene_set)
    else:
        i = 0
        while i < len(gene_set):
            print(gene_set[i].record())
            i += 1


read_file()
get_gene()
