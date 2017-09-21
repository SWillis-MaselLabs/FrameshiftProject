Artificial Frameshifter ReadMe:

The purpose of this program is to take in data from a tab-delimited file, and return two different nucleotide sequences in tab-delimited format. The two sequences will be the original sequence artificially frameshifted into the +1 and +2 frame. The results will be printed to two tab-delimited files. For a somewhat cleaner workspace, the output files will be saved one directory up from where the python script is stored.

The output file will be printed in the same format with the same entry layout as the input file

The format of the input/output files are tab-delimited and require the following entries:

1 - UID
2 - Gene Family
3 - Species
4 - Genome Accession Number
5 - Protein Accession Number
6 - Gene Designation
7 - Gene Name
8 - Nucleotide Sequence
9 - Protein Sequence
10 - Protein Sequence (No Cysteines)

_____________________________________________________________________________

This program functions by doing the following (with examples) to nucleotide sequences. Spaces are used to aid the reader and are not present in the actual code

1) Takes in a nucleotide sequence

ATG CGA AGG TTA GCA TAG

2) Removes the start and stop codon from the beginning and end of the sequence and stores them for later

[stored: ATG, TAG]

CGA AGG TTA GCA

3) Removes first (first and second) nucleotide(s) from the sequences for an artificial plus one (plus two) frameshift (plus one frameshift shown):

[stored: ATG, TAG]

GAA GGT TAG CA

4) Removes the minimum number of nucleotides from the end of the sequence so that the length of the sequences remains a multiple of 3:

[stored: ATG, TAG]

GAA GGT TAG

5) Removes any start codons that may have arisen in the body of the sequence as a result of the frameshift

[stored: ATG, TAG]

GAA GGT

6) Reattaches the start and stop codons to the end of the sequence

ATG GAA GGT TAG
