OverlapByReadingFrame ReadMe

The purpose of this file is to take in the nucleotide sequences of genes and a subsequence (for the purposes of this research, the subsequence that is shared by the input gene and an overlapping partner), and to trim the subsequence so that it’s length is a multiple of 3 and it is in the same reading frame as the overall gene.

The main program is called “OverlapByReadingFrame.py” and will take in a tab-delimited file with the following entries:

1) Sequence UID
2) Overlapping Sequence (n+)
3) Gene Sequence (n+)

Once the overlapping sequence is found in the correct frame, with the correct length, it will be printed to an output file:

	“OverlappingSequences_InFrame.txt”

which will be located in the same directory as the program. It will be printed in tab-delimited format with the following entries:

1) UID
2) Overlapping Sequence (n+) (correct frame and length)