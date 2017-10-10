# The purpose of this program is to take in an gene which overlaps with a second gene
# along with the shared nucleotide sequence and outputs that shared nucleotide sequence
# in the same reading frame as the gene that was read in.
# This is done by removing the minimum number of
# nucleotides from the shared nucleotide sequence so that it is forced into the reading
# frame of the overall sequence. This is done so that the codon usage of the overlapping
# sequence may be analyzed by reading frame. 


# The output of this file is in tab-delimited format, with the unique ID being printed as
# the first output followed by the overlapping sequence in the proper reading frame. 

import csv
import os
from pathlib import Path
import time

# Here the filename of the input file is declared

filename = 'OverlappingVirusGenes_92PairsWithSequenceData.txt'

OutputFilename = Path('./OverlappingSequences_InFrame.txt')

if OutputFilename.is_file():
    os.remove('%s'%OutputFilename)
time.sleep(0.5)


# The csv module is used to read in the file declared earlier. The file should
# be in tab-delimited format. 
with open('%s'%filename, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:

        # The file that's read in should have the following entries:
        # 1 - Unique ID associated with that particular gene
        # 2 - Unique pair ID, used to match overlapping genes. Each pair has a unique ID associated with it
        # 3 - The name of the viral species where the genes came from
        # 4 - The name of the gene (one of the genes in the overlapping pair)
        # 5 - The shared nucleotide sequence associated with both genes involved in the overlap
        # 6 - The nucleotide sequence of the gene (associated with the gene named in entry 4)
        SequenceUID = row[0]
        OverlappingSequence = row[1]
        CodingGene = row[2]

        # If the overlapping sequence and the coding gene are the same, then there isn't anything
        # that needs to be done to put the overlapping sequence into the reading frame of the
        # gene in question. (This happens when one gene in an overlapping pair is entirely contained
        # within its overlapping partner)

        if CodingGene == OverlappingSequence:
            OverlappingSequence_InFrame = OverlappingSequence
        # If the overlapping region and the overall gene aren't the same, then the overlapping region needs to
        # be trimmed so that it's in the reading frame of the gene
        
        else:
            # To determine how to trim the overlapping sequence, the index where the overlapping region starts
            # in the gene is found. The index then is found modulo 3. The indexing works naturally because
            # python starts counting at 0, so if the index is 0 mod 3, then the overlapping region is in the reading
            # frame of the gene and no nucleotides need to be removed from the start of the sequence. If the index is
            # 1 mod 3, then the overlapping sequence is in the +1 reading frame of the gene, and the minimum number of
            # nucleotides that need to be removed from the sequence to force it into the reading frame of the gene is 2.
            # Similarly, if the index is 2 mod 3, one nucletide needs to be removed.  
            StartIndex = CodingGene.index(OverlappingSequence)
            FrameshiftType = StartIndex%3
            if FrameshiftType == 0:
                OverlappingSequence_InFrame = OverlappingSequence
            else:
                if FrameshiftType == 1:
                    # Here, two nucleotides are removed from the start overlapping sequence using string splicing
                    OverlappingSequence_InFrame = OverlappingSequence[2::]
                else:
                    if FrameshiftType == 2:
                        # Similarly, here one nucleotide is removed from the beginning of the overlapping sequence
                        OverlappingSequence_InFrame = OverlappingSequence[1::]
        # Before the final sequence is printed, the overlapping sequence which has been put into the reading frame
        # of the gene in question needs to be shortened so that its length is 0 mod 3. 
        OverlappingLength = len(OverlappingSequence_InFrame)
        OverlappingLength_mod3 = OverlappingLength%3
        OutputOverlappingSequence = OverlappingSequence_InFrame[:OverlappingLength - OverlappingLength_mod3:]

        OutputFile = open('%s' %OutputFilename, 'a')
        output = [SequenceUID, OutputOverlappingSequence]
        OutputFile.write('\t'.join(output[0:]) + '\n')
        OutputFile.close()
