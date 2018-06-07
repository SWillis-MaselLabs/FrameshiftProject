# The purpose of this program is to take in coding sequences and to find all
# open reading frames in the +1 and +2 reading frames. If multiple ORFs are
# found to terminate at a common stop codon, only the longest is returned
# to the user.
import csv
import sys
from Bio.Seq import Seq
import Bio

# The orfs in each coding sequence are saved to lists, distinguished by which
# reading frame they were found in
PlusOneORFs = []
PlusTwoORFs = []

# As a baseline, if an ORF is fewer than 21 nucleotides in length (translating to
# 7 amino acids), it is discarded. Analyses done with these ORFs later implemented
# a stricter length cutoff (25 AA)
lengthCutoff = 21

# The tab-delimited file containing the coding sequences is parsed
filename = 'VirusGenes_Nonoverlapping_Controls_AncestralHomologs.txt'
with open('%s' %filename) as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        UID = row[0]
        PairUID = row[1]
        Genus = row[2]
        Species = row[3]
        GeneName = row[4]
        Accession = row[5]
        ControlDesignation = row[6]
        CodingSeq = row[7]

        # The sequence in the plus one frame is found by removing the first nucleotide
        # and the last two nucleotides in the sequence
        plusOne = CodingSeq[1:len(CodingSeq)-2:]
        plusOneLength = len(plusOne)
        # The indices that correspond to the starting index of each codon are then
        # found
        Index1 = range(0, plusOneLength, 3)

        # Similarly, the sequence in the plus two frame is found by removing the first
        # two and last nucleotide in the sequence
        plusTwo = CodingSeq[2:len(CodingSeq)-1:]
        plusTwoLength = len(plusTwo)
        # And the indices that correspond to the starting index of each codon are declared
        Index2 = range(0, plusTwoLength, 3)

        # Each codon in the plus one frame is then found using the indices declared
        for n in Index1:
            m = n+3
            codon = plusOne[n:m]
            # The ORF is declared as an empty string that will have codons appended to it
            # once a start codon is found
            ORF = ''
            if codon == 'ATG':
                j = n
                # Once a start codon is found, the ORF has all subsequence codons appended
                # to it until a stop codon is found
                ORF += codon
                for j in Index1[j:]:
                    k = j+3
                    newCodon = plusOne[j:k]
                    # If a stop codon is encountered, the stop is appended to the ORF and
                    # the nested for loop is exited so the rest of the sequence can be scanned
                    # for a new start codon and a new ORF is created.
                    if newCodon =='TAG' or newCodon == 'TAA' or newCodon == 'TGA':
                        ORF += newCodon
                        # if the ORF is below the length cutoff, it is discarded
                        if len(ORF) <= lengthCutoff:
                            pass
                        else:
                            # otherwise, it is stored in a list of ORFs with the UID of the gene
                            # that it was found in
                            PlusOneORFs.append([UID,ORF])
                        break
                    else:
                        ORF += newCodon

        # The same procedure is performed for the ORFs in the +2 reading frame.
        for n in Index2:
            m = n+3
            codon = plusTwo[n:m]
            ORF = ''
            if codon == 'ATG':
                j = n
                ORF += codon
                for j in Index2[j:]:
                    k = j+3
                    newCodon = plusTwo[j:k]
                    if newCodon =='TAG' or newCodon == 'TAA' or newCodon == 'TGA':
                        ORF += newCodon
                        if len(ORF) <= lengthCutoff:
                            pass
                        else:
                            PlusTwoORFs.append([UID,ORF])
                        break
                    else:
                        ORF += newCodon

# The ORFs in the plus one frame are then compared to the rest of the ORFs that were found
# If any ORF is a subsequence of another, it is removed
for orf in PlusOneORFs:
    for orf2 in PlusOneORFs:
        if orf[1] in orf2[1]:
            if orf in PlusOneORFs:
                PlusOneORFs.remove(orf)
        elif orf2[1] in orf[1]:
            if orf2 in PlusOneORFs:
                PlusOneORFs.remove(orf2)

# The same is done for the +2 frame
for orf in PlusTwoORFs:
    for orf2 in PlusTwoORFs:
        if orf[1] in orf2[1]:
            if orf in PlusTwoORFs:
                PlusTwoORFs.remove(orf)
        elif orf2[1] in orf[1]:
            if orf2 in PlusTwoORFs:
                PlusTwoORFs.remove(orf2)

                
# The user can choose to have the number of ORFs found in each reading frame
# printed for inspection

#print('Number of +1 ORFs: %s' %len(PlusOneORFs))
#print('Number of +2 ORFs: %s' %len(PlusTwoORFs))

# The results are printed for the user. The files can be saved to a tab-delimited file
# by the python3 GeneratePreexistingORFControls.py >[filename] in the terminal. 
for orf in PlusOneORFs:
    orfSeq = orf[1]
    protein = Seq(orfSeq).translate()
    protein = str(protein)
    noCysProtein = protein.replace('C','')
    noCysProteinLen = len(noCysProtein)
    output = [orf[0],orf[1], protein, noCysProtein, '+1', noCysProteinLen]
    print(*output,sep='\t')


for orf in PlusTwoORFs:
    orfSeq = orf[1]
    protein = Seq(orfSeq).translate()
    protein = str(protein)
    noCysProtein = protein.replace('C','')
    noCysProteinLen = len(noCysProtein)
    output = [orf[0],orf[1], protein, noCysProtein, '+2', noCysProteinLen]
    #print(*output,sep='\t')
