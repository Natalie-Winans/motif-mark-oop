#!/usr/bin/env python

import argparse
import os
import re
import cairo

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_args():
    parser = argparse.ArgumentParser(
        prog="MotifMark",
        description="Tool for visualization of motifs on sequences",
        add_help=True)
    parser.add_argument(
        '-f', '--fasta', help="input FASTA file containing sequences of interest")
    parser.add_argument(
        '-m', '--motifs', help="input file containing motifs of interest")
    return parser.parse_args()


def extract_fasta_name(fasta_file):
    """Extract input filename without extension for figure title and output files."""
    fasta_name = os.path.basename(fasta_file).split('.')[0]
    return fasta_name


def convert_motifs(motif_file):
    """Take text file of motifs, convert to regex statements based on IUPAC conventions, 
    return list of converted motifs ready to be searched with regex."""
    base_dict = {
    'A': '[Aa]',
    'C': '[Cc]',
    'G': '[Gg]',
    'T': '[TtUu]',
    'U': '[UuTt]',
    'W': '[AaTtUu]',
    'S': '[CcGg]',
    'M': '[AaCc]',
    'K': '[GgTtUu]',
    'R': '[AaGg]',
    'Y': '[CcTtUu]',
    'B': '[CcGgTtUu]',
    'D': '[AaGgTtUu]',
    'H': '[AaCcTtUu]',
    'V': '[AaCcGg]',
    'N': '[AaCcGgTtUu]',
    'Z': '[]'}
    motif_dict = {}
    with open(motif_file, 'r') as fh:
        for motif in fh:
            re_motif = ""
            motif = motif.upper().strip()
            for base in motif:
                re_motif += base_dict[base]
            motif_dict[motif] = re_motif
        return motif_dict


def palette(num_motifs):
    """Take number of motifs (up to 6) and return 
    colorblind-friendly color palette."""
    if num_motifs == 1:
        pal = [.13, .53, .2, .8]
    if num_motifs == 2:
        pal = [[.13, .53, .2, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 3:
        pal = [[.27, .47, .67, .8],
               [.13, .53, .2, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 4:
        pal = [[.4, .8, .93, .8],
               [.13, .53, .2, .8],
               [.8, .73, .27, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 5:
        pal = [[.27, .47, .67, .8],
               [.4, .8, .93, .8],
               [.13, .53, .2, .8],
               [.8, .73, .27, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 6:
        pal = [[.27, .47, .67, .8],
               [.4, .8, .93, .8],
               [.13, .53, .2, .8],
               [.8, .73, .27, .8],
               [.93, .4, .47, .8],
               [.67, .2, .47, .8]]
    return pal


def parse_fasta(fasta_file):
    """Take fasta file and yield header and seq string objects.
    (generator function adapted from biopython)"""
    header, seq = None, []
    for line in fasta_file:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield(header, ''.join(seq))
            header, seq = line, []
        else:
            seq.append(line)
    if header:
        yield(header, ''.join(seq))


def gene_dictionary(fasta_file):
    """Take fasta file and return dictionary with header lines 
    as keys and concatenated sequences as values."""
    with open(fasta_file, 'r') as fh:
        gene_dict = {}
        for header, seq in parse_fasta(fh):
            gene_dict[header] = seq

    return(gene_dict)

def motif_spans(motif_dict, gene_dict):
    """
    Take motif file and list of gene/pre-mRNA sequences,
    return dictionary of dictonaries of lists of tuples.

    Outer dictionary:
        keys = gene/pre-mRNA sequence indexes
        values = inner dictionaries:
            keys = motifs
            values = lists of tuples representing span of 
                     motif in the sequence
    """
    all_motif_spans = {}
    for header, seq in gene_dict.items():
        seq_motif_spans = {}
        for j, motif in enumerate(motif_dict):
            motif_match = re.finditer(motif_dict[motif], seq)
            seq_motif_spans[motif] = [m.span() for m in motif_match]
        #print(seq_motif_spans)

        all_motif_spans[header] = seq_motif_spans

    return all_motif_spans


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CLASSES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class FastaHeader:
    """
    A class to represent a FASTA file header line.

    Attributes:
        text (str): text of a header line beginning with '>'
        gene_number (int): the gene index in the FASTA file
    """
    def __init__(self, text, gene_number):
        self.text = text
        self.gene_number = gene_number

    def write(self, ctx):
        """Write FASTA header above each gene object."""
        y = GENE_GROUP_HEIGHT * self.gene_number + GENE_Y_OFFSET
        ctx.move_to(LEFT_MARGIN, y - 35)
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(16)
        ctx.show_text(self.text)


class Gene:
    """A class to represent one gene/entry in a FASTA file.
    
    Attributes:
        length (int): length in base pairs of the gene sequence
        gene_number (int): the gene index in the FASTA file
    """
    line_width = 5
    def __init__(self, length, gene_number):
        self.length  = length
        self.gene_number = gene_number

    def draw(self, ctx):
        """Draw line representing gene to scale."""
        y = GENE_GROUP_HEIGHT * self.gene_number + GENE_Y_OFFSET
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_line_width(self.line_width)
        ctx.move_to(LEFT_MARGIN, y)
        ctx.line_to(LEFT_MARGIN + self.length, y)
        ctx.stroke()

class Exon:
    """A class to represent an exon within a gene/entry.

    Attributes:
        start (int): value representing the starting position
        end (int): value representing the ending position
        gene_number (int): the gene index in the FASTA file
    """
    line_width = 35
    def __init__(self, start, end, gene_number):
        self.start = start
        self.end = end
        self.gene_number = gene_number

    def draw(self, ctx):
        """Draw line representing the exon to scale."""
        y = GENE_GROUP_HEIGHT * self.gene_number + GENE_Y_OFFSET
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_line_width(self.line_width)
        ctx.move_to(LEFT_MARGIN + self.start, y)
        ctx.line_to(LEFT_MARGIN + self.end, y)
        ctx.stroke()
        

class Motif:
    """A class to represent a motif to draw within a gene/sequence.
    
    Attributes:
        start (int):value representing the starting position
        end (int): value representing the ending position
        gene_number (int): the gene index in the FASTA file
        num_motifs (int): number of motifs in input file
    """
    line_width = 35
    def __init__(self, spans, gene_number, num_motifs):
        # self.start = start
        # self.end = end
        self.spans = spans
        self.gene_number = gene_number
        self.pal = palette(num_motifs)
    

    def draw(self, ctx):
        """Draw all instances of motif found in gene to scale.""" 
        y = GENE_GROUP_HEIGHT * self.gene_number + GENE_Y_OFFSET
        ctx.set_line_width(self.line_width)
        # ctx.set_source_rgba(*self.pal[j])
        # ctx.set_source_rgb(0, 0, 0)
        # ctx.move_to(LEFT_MARGIN + self.start, y)
        # ctx.line_to(LEFT_MARGIN + self.end, y)
        # ctx.stroke()

        for span in self.spans:
            ctx.move_to(LEFT_MARGIN + span[0], y)
            ctx.set_line_width(35)
            ctx.line_to(LEFT_MARGIN + span[1], y)
            ctx.stroke()


class GeneGroup:
    """A class to represent a gene group.
    
    Attributes:
        header: Header object
        gene: Gene object
        exon: Exon object
        motif_spans: dictionary containing ...
    """
    def __init__(self, header, gene, exon, found_motifs):
        self.header = header
        self.gene = gene
        self.exon = exon
        self.found_motifs = found_motifs

    def draw_gene_group(self):
        """Draw gene group including introns, exon, motifs, and header line."""
        self.header.write(ctx)
        self.gene.draw(ctx)
        self.exon.draw(ctx)

        for i, motif in enumerate(self.found_motifs):
            spans = found_motifs[motif]
            motif = Motif(spans, gene_num, num_motifs)

            motif.draw(ctx)
        



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN CODE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

args = get_args()
fasta_file = args.fasta
motif_file = args.motifs

# global constants
GENE_GROUP_HEIGHT = 150
GENE_Y_OFFSET = 150
LEFT_MARGIN = 50

# for figure title
fasta_name = extract_fasta_name(fasta_file)

# read in information from input files
gene_dict = gene_dictionary(fasta_file)
motif_dict = convert_motifs(motif_file)
num_motifs = len(motif_dict)

# motif_spans = motif_spans(motif_dict, gene_dict)


# set surface dimensions based on number of sequences and max length
### make this a function too?
num_genes = len(gene_dict)
max_length = len(max(gene_dict.values(), key = len))
width = max_length + 250
height = num_genes * GENE_GROUP_HEIGHT + GENE_Y_OFFSET
surface = cairo.SVGSurface('%s.svg' % fasta_name, width, height)
ctx = cairo.Context(surface)
ctx.rectangle(0, 0, width, height)
ctx.set_source_rgba(1, 1, 1, 1)
ctx.fill()

# write figure title from fasta name
ctx.move_to(LEFT_MARGIN, 60)
ctx.set_source_rgb(0, 0, 0)
ctx.set_font_size(24)
ctx.show_text("FASTA file: %s" % fasta_name)

# draw gene groups
gene_num = 0
for header, seq in gene_dict.items():
    header = FastaHeader(header, gene_num)

    gene = Gene(len(seq), gene_num)

    exon_start = re.search('[A-Z]+', seq).span()[0]
    exon_end = re.search('[A-Z]+', seq).span()[1]
    exon = Exon(exon_start, exon_end, gene_num)

    found_motifs = {}
    for motif in motif_dict:
        motif_match = re.finditer(motif_dict[motif], seq)
        found_motifs[motif] = [m.span() for m in motif_match]
    print(found_motifs)

    group = GeneGroup(header, gene, exon, found_motifs)
    group.draw_gene_group()

    # draw motifs
    
    gene_num += 1


surface.write_to_png('%s.png' % fasta_name)






