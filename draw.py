from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def draw_gene(sequence):
    """
    (genbank file) - > image (pdf, png etc)
    
    """

    record = SeqIO.read(sequence, "genbank")
    diagram = GenomeDiagram.Diagram(record.id)
    feature_track = diagram.new_track(1, name="Annotated Features")
    feature_set = feature_track.new_set()

    for feature in record.features:
        if feature.type != "gene":
            #Exclude this feature
            continue
        if len(feature_set) % 2 == 0:
            color = colors.blue
        else:
            color = colors.lightblue
        feature_set.add_feature(feature, sigil="ARROW", label_size=14,color=color,
                                label=True)

    diagram.draw(format="cirular", circular=True,
                 pagesize=(50*cm,50*cm), fragments=1, start=0,
                 end=len(record), circle_core=0.5)
    diagram.write("plasmid.pdf", "PDF")

def draw_chromosome(sequence):
    entries = [("Legionella Pneumophilia")]
    max_len = 30432563
    telomere_length = 1000000

    chr_diagram = BasicChromosome.Organism()
    chr_diagram.page_size = (29.7*cm, 21*cm) #A4 landscape

    for name, length in entries:
        cur_chromosome = BasicChromosome.Chromosome(name)
        cur_chromosome.scale_num = max_len + 2 * telomere_length

        start = BasicChromosome.TelomereSegment()
        start.scale = telomere_length
        cur_chromosome.add(start)

        body = BasicChromosome.ChromosomeSegment()
        body.scale = length
        cur_chromosome.add(body)

        end = BasicChromosome.TelomereSegement(inverted=True)
        end.scale = telomere_length
        cur_chromosome.add(end)

        chr_diagram.add(cur_chromosome)
    
    chr_diagram.draw("Chromosome.pdf", "Legionella Pneumophilia")

draw_gene("legionella.gb")











        
    
