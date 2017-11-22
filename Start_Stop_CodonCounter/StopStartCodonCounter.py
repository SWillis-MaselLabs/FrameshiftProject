import csv
import numpy as np
from scipy import stats
import sys

# The purpose of this file is to go through a flat file containing
# non-overlapping viral genes and to count the occurences of start and
# stop codons in the gene's alternative reading frames.

# Specify the name of the flat file containing the genes
filename = 'VirusGenes_Nonoverlapping_Controls_Complete.txt'

# All counters are defined and set to 0

# The first set of counters will keep track of the total number of
# codons in the +1 and +2 frames
PlusOneCodons = 0
PlusTwoCodons = 0

# The second set will keep track of the number of stop codons
# that occur in the +1 and +2 frames, respectively
PlusOneStops = 0
PlusTwoStops = 0

# The final set will keep track of the number of start codons that
# occur in the +1 and +2 frames, respectively
PlusOneStarts = 0
PlusTwoStarts = 0



# The flat file is opened using the csv module. The flat file should be
# in tab-delimited format. If it's not, the delimiter can be changed
# to match the format.
with open('%s'%filename, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        # The first element in the flat file is the database UID, followed by the
        # gene's nucleotide sequence.
        UID = row[0]
        Sequence_0 = row[1]
        # The nucleotide sequence length is found
        GeneLength = len(Sequence_0)
        # it is then used to artificially shift the sequence into its +1 frame.
        # The sequence is ended two nucleotides prematurely to make sure that the
        # sequence stays a multiple of 3 in length (the final two nucleotides don't
        # constitute a codon). The same procedure is used to shift the sequence
        # into its +2 frame, except it is truncated one nucleotide from the end instead
        # of two.
        Sequence_1 = Sequence_0[1:GeneLength-2]
        Sequence_2 = Sequence_0[2:GeneLength-1]

        # +1 Frame Counting block
        
        # The length of the sequence in its +1 frame is noted
        Length_1 = len(Sequence_1)
        # Then a list of indices is compiled ranging from 0 to the length of the sequence
        # using a step size of 3. This translates to a list of indices that correspond to
        # the first nucleotide in each codon in the sequence. 
        Index_1 = range(0,len(Sequence_1),3)
        # Scanning the nucleotide sequence by running through the indices previously collected
        for n in Index_1:
            # m is defined to be n+3 so that Sequence[n:m] is the codon starting at index n
            m = n+3
            Codon = Sequence_1[n:m]
            # 1 is added to the total codon counter
            PlusOneCodons +=1
            # If the codon is a stop codon, it is counted by the relevant stop codon counter
            if Codon == 'TAA' or Codon == 'TAG' or Codon == 'TGA':
                PlusOneStops += 1
            # If the codon is a start codon, it is counted by the relevant start codon counter
            if Codon == 'ATG':
                PlusOneStarts += 1



        # +2 Frame Counting - this is the same as the previous block, except is for the gene
        # in its +2 reading frame
        Length_2 = len(Sequence_2)
        Index_2 = range(0,len(Sequence_2),3)
        for n in Index_2:
            m = n+3
            Codon = Sequence_2[n:m]
            PlusTwoCodons += 1
            if Codon == 'TAA' or Codon == 'TAG' or Codon == 'TGA':
                PlusTwoStops += 1
            if Codon == 'ATG':
                PlusTwoStarts += 1

# The results are then printed for the user
print()
print('Plus One Frame:')
print('Total Number of Codons: %s' %PlusOneCodons)
print('Number of Start Codons: %s' %PlusOneStarts)
print('Percentage of Total: %.2f%%' %((PlusOneStarts/PlusOneCodons)*100))
print('1 in %d +1 Codons are Starts\n' %(PlusOneCodons/PlusOneStarts))
print('Number of Stop Codons: %s' %PlusOneStops)
print('Percentage of Total: %.2f%%' %((PlusOneStops/PlusOneCodons)*100))
print('1 in %d +1 Codons are Stops\n' %(PlusOneCodons/PlusOneStops))

print()
print('Plus Two Frame:')
print('Total Number of Codons: %s' %PlusTwoCodons)
print('Number of Start Codons: %s' %PlusTwoStarts)
print('Percentage of Total: %.2f%%' %((PlusTwoStarts/PlusTwoCodons)*100))
print('1 in %d +2 Codons are Starts\n' %(PlusTwoCodons/PlusTwoStarts))
print('Number of Stop Codons: %s' %PlusTwoStops)
print('Percentage of Total: %.2f%%' %((PlusTwoStops/PlusTwoCodons)*100))
print('1 in %d +2 Codons are Stops\n' %(PlusTwoCodons/PlusTwoStops))

