# Motif Mark (OOP)
### A Python script to plot protein-binding motifs on an image of an exon and its 5' and 3' introns

* Handles up to six motifs, may contain degenerate bases (IUPAC symbols).
* Outputs a single figure in svg and png formats.

Call `motif-mark-oop.py` from the command line using the following arguments:

```
-f, --fasta     FASTA file containing sequences of interest
-m, --motifs    Text file containing up to six motifs, one motif per line
-h, --help      Display help message describing arguments

Example output:

![](https://github.com/Natalie-Winans/motif-mark-oop/blob/main/Figure_1.svg)
