import csv
from pathlib import Path
import time
import os

# Author: Sara Willis

# The purpose of this file is to find the CAI values of gene segments from a 
# tab-delimited file. The CAI, or Codon Adaptation Index, measures the degree of
# to which a gene or gene segment exhibits the same codon bias as the rest of the species

# In the case of this script, the CAI values correspond to sections of overlapping genes
# (in the frame of the gene that contains them). 

# This script takes in two separate files to determine the CAI values. The first is a
# file containing the relative adaptedness values calculated from the raw codon count
# from the species. There is one relative adaptedness value per codon

# The second file contains the gene segments for analysis.

# Here the user may declare the filenames. 

RelativeAdaptednessFile = 'RelativeAdaptednessValues.txt'
GeneSequencesFile = 'GeneSequencesFile.txt'

outputFile = Path('../CodonAdaptationIndexFile.txt')

# The csv module is used to read in the GeneSequencesFile in tab-delimited format
# The entries of the file should be the following:
# 1) The sequence UID corresponding the gene to be analyzed
# 2) The viral species to which the gene belongs
# 3) The name of the gene in question to which the overlapping segment belongs
# 4) The nucleotide sequence of the segment. Note: This sequence is the
# Overlapping sequence, not the entire sequence. The overlapping
# sequence should be trimmed so that it's in the same reading frame as the
# overall gene. 

if outputFile.is_file():
    os.remove('%s'%outputFile)
    time.sleep(0.5)

with open('%s' %GeneSequencesFile, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        UID = row[0]
        Species = row[1]
        GeneName = row[2]
        Sequence = row[3]

        # So that we can access, store and count all the different codons in a sequence
        # dictionaries are created. The number at the beginning of each sub-list is used
        # for identification purposes, i.e. calling a particular sublist and matching it to
        # the same location in another dictionary

        # CodonIndex stores the names of the codons
        # CodonCount is used to store the raw codon counts
        # RelativeAdaptedness stores the relative adaptedness values for each codon. There are
        # unique relative adaptedness values for each codon, and it varies by viral species. These
        # Values will come from the RelativeAdaptednessFile

        
        CodonIndex = [[0,'TTT','TTC'],[1,'TTA','TTG','CTT','CTC','CTA','CTG'],[2,'ATT','ATC','ATA'],[3,'ATG'],[4,'GTT','GTC','GTA','GTG'],[5,'TCT','TCC','TCA','TCG'],[6,'CCT','CCC','CCA','CCG'],[7,'ACT','ACC','ACA','ACG'],[8,'GCT','GCC','GCA','GCG'],[9,'TAT','TAC'],[10,'TAA','TAG','TGA'],[11,'CAT','CAC'],[12,'CAA','CAG'],[13,'AAT','AAC'],[14,'AAA','AAG'],[15,'GAT','GAC'],[16,'GAA','GAG'],[17,'TGT','TGC'],[18,'TGG'],[19,'CGT','CGC','CGA','CGG','AGA','AGG'],[20,'AGT','AGC'],[21,'GGT','GGC','GGA','GGG']]

        CodonCount = [[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]


        RelativeAdaptedness =[[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]

            
        GeneRelativeAdaptedness = []

        Index = ['Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr','Stop','His','Gln','Asn','Lys','Asp','Glu','Cys','Trp','Arg','Ser','Gly']

        # The variable TotalCodonCount is defined to count the total number of codons in each gene segment
        TotalCodonCount = 0

        # Sequence1Length gives the length of the nucleotide segment 
        Sequence1CodonNumber = int(len(Sequence))

        # Next, the number of codons in the sequence is found so that the sequence is properly subdivided
        # so that every codon is correctly accounted for. I.e. we need to splice up the nucleotide segment
        # into blocks of three (the length of a codon). The sequence is then scanned and the number of each
        # type of codon is counted and stored in the correct dictionary
        for n in range(0,Sequence1CodonNumber,3):
            TotalCodonCount += 1
            codon = Sequence[n:n+3]
            for entry in CodonIndex:
                if codon in entry:
                    external_index = CodonIndex.index(entry)
                    internal_index = entry.index(codon)
                    CodonCount[external_index][internal_index]=CodonCount[external_index][internal_index]+1

        # The next step is to open the file that contains the relative adaptedness values
        # The indices are as defined below
        with open('%s' %RelativeAdaptednessFile, 'r') as f:
            reader = csv.reader(f, delimiter = '\t')
            for row in reader:
                RelativeAdaptednessValue_UID = row[0]
                VirusName = row[1]
                RelativeAdapted_Designation = row[2]
                if VirusName != Species:
                    pass
                else:
                    row = [0 if x=='' else x for x in row]
                    
                    UUU = float(row[3])
                    UUC = float(row[4])
                    UUA = float(row[5])
                    UUG = float(row[6])
                    UCU = float(row[7])
                    UCC = float(row[8])
                    UCA = float(row[9])
                    UCG = float(row[10])
                    UAU = float(row[11])
                    UAC = float(row[12])
                    UAA = float(row[13])
                    UAG = float(row[14])
                    UGU = float(row[15])
                    UGC = float(row[16])
                    UGA = float(row[17])
                    UGG = float(row[18])
                    CUU = float(row[19])
                    CUC = float(row[20])
                    CUA = float(row[21])
                    CUG = float(row[22])
                    CCU = float(row[23])
                    CCC = float(row[24])
                    CCA = float(row[25])
                    CCG = float(row[26])
                    CAU = float(row[27])
                    CAC = float(row[28])
                    CAA = float(row[29])
                    CAG = float(row[30])
                    CGU = float(row[31])
                    CGC = float(row[32])
                    CGA = float(row[33])
                    CGG = float(row[34])
                    AUU = float(row[35])
                    AUC = float(row[36])
                    AUA = float(row[37])
                    AUG = float(row[38])
                    ACU = float(row[39])
                    ACC = float(row[40])
                    ACA = float(row[41])
                    ACG = float(row[42])
                    AAU = float(row[43])
                    AAC = float(row[44])
                    AAA = float(row[45])
                    AAG = float(row[46])
                    AGU = float(row[47])
                    AGC = float(row[48])
                    AGA = float(row[49])
                    AGG = float(row[50])
                    GUU = float(row[51])
                    GUC = float(row[52])
                    GUA = float(row[53])
                    GUG = float(row[54])
                    GCU = float(row[55])
                    GCC = float(row[56])
                    GCA = float(row[57])
                    GCG = float(row[58])
                    GAU = float(row[59])
                    GAC = float(row[60])
                    GAA = float(row[61])
                    GAG = float(row[62])
                    GGU = float(row[63])
                    GGC = float(row[64])
                    GGA = float(row[65])
                    GGG = float(row[66])
                    
                    RelativeAdaptedness = [[0,UUU,UUC],[1,UUA,UUG,CUU,CUC,CUA,CUG],[2,AUU,AUC,AUA],[3,AUG],[4,GUU,GUC,GUA,GUG],[5,UCU,UCC,UCA,UCG],[6,CCU,CCC,CCA,CCG],[7,ACU,ACC,ACA,ACG],[8,GCU,GCC,GCA,GCG],[9,UAU,UAC],[10,UAA,UAG,UGA],[11,CAU,CAC],[12,CAA,CAG],[13,AAU,AAC],[14,AAA,AAG],[15,GAU,GAC],[16,GAA,GAG],[17,UGU,UGC],[18,UGG],[19,CGU,CGC,CGA,CGG,AGA,AGG],[20,AGU,AGC],[21,GGU,GGC,GGA,GGG]]

                    # Finally, the raw codon usage data for each of the gene segments
                    # Plus the relative adapteness values are used to calculate the
                    # CAI. The CAI is defined to be the geometric mean of the relative
                    # adaptedness values, i.e. (in LaTex formatting)
                    # CAI = (\prod_{i=1}^{L} w_i)^{1/L}
                    # i.e. the product from i=1 to L, where L is the number of codons
                    # in the segment in question, and w_i is the relative adaptedness
                    # value for the codon i.
                    
                    CodonAdaptationIndex = 1
                    for entry in RelativeAdaptedness:
                        CountIndex = entry[0]
                        entry = entry[1:]
                        for item in entry:
                            numberOfCodons = CodonCount[CountIndex][entry.index(item)+1]
                            #If there are no codons of a particular type, then the relative
                            # adaptedness value is defined to be 0.5 
                            if item == 0.0:
                                item = 0.5
                            CodonAdaptationIndexIntermediate = item**numberOfCodons
                            CodonAdaptationIndex = CodonAdaptationIndex*CodonAdaptationIndexIntermediate
                    # In the even there is not enough information to complete the CAI analysis
                    # for instance, there are no sequence data available, then the CAI is printed as [Null]
                    if TotalCodonCount-CodonCount[3][1]-CodonCount[18][1] <= 0:
                        CodonAdaptationIndex = '[Null]'
                    else:
                        # Otherwise, The CodonAdaptationIndex (the product of all the relative adaptedness values
                        # not yet being raised to the power of 1/L) is raised to the power of the total number of codons,
                        # minutes the nymber of codons with degeneracy (i.e. amino acids which are coded for by only one
                        # amino acid). These are removed as they introduce no new information and so their relative adaptedness
                        # values are always one
                        CodonAdaptationIndex = CodonAdaptationIndex**(1/(TotalCodonCount-CodonCount[3][1]-CodonCount[18][1]))
                        
                    # The relevant information is then printed out in tab-delimited format. The UID first
                    # so the data can be uploaded to the MySQL database, the species the nucleotide belongs
                    # to, the name of the gene the segment was in (and the reading frame it was forced into)
                    # and finally the codon adaptation index

                    # With these data, gene segments are either tentatively classified as ancestral or novel
                    # Higher CAI ==> ancestral (gene more closely resembles overall genome codon usage patterns)
                    # Lower CAI ==> novel (gene less closely resembles overall genome codon usage patterns)

                    output = [UID, Species, GeneName, str(CodonAdaptationIndex)]

                    WriteToOutput = open('%s' %outputFile, 'a')
                    WriteToOutput.write('\t'.join(output[0:])+'\n')
                    WriteToOutput.close()

                    #print('%s\t%s\t%s\t%s' %(UID,Species,GeneName,CodonAdaptationIndex))
