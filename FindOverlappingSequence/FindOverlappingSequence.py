from Bio import SeqIO
from Bio.Seq import Seq
from LongestCommonSubstring import longest_common_substring

# The purpose of this file is to take in two overlapping
# nucleotide sequences in FASTA format and to find the
# nucleotide sequence that is shared by both.

# The two FASTA files (found in the directory where this script
# is stored) are read in using the SeqIO function from the
# Bio module. The sequences are stored as sequence1 and
# sequence2 respectively.
for seq_record in SeqIO.parse('Sequence1.fasta','fasta'):
    sequence1 = seq_record.seq
    SeqRecord1 = seq_record.id

for seq_record in SeqIO.parse('Sequence2.fasta', 'fasta'):
    sequence2 = seq_record.seq
    SeqRecord2 = seq_record.id

# We import the submodule LongestCommonSubstring and use the
# function longest_common_substring, which takes in two strings
# and returns the longest shared substring. In this case, that's
# the nucleotide sequence shared by both genes.
Overlap = longest_common_substring(sequence2, sequence1)

# The results are then printed with the two original sequences
# and their corresponding data from the FASTA file, along with the
# overlapping sequence and its length.
print('Sequence 1: %s' %SeqRecord1)
print(sequence1)
print('Sequence 2: %s' %SeqRecord2)
print(sequence2)
print('Overlapping Sequence:')
print(Overlap)
print('Length of Overlap: %s' %len(Overlap))

