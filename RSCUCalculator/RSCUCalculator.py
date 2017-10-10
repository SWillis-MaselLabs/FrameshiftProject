import numpy
import csv
import os
from pathlib import Path
import sys
import time



# The purpose of this script it to take in codon-usage data in a tab-delimited format
# to calculate both the RSCU values as well as the relative adaptedness.

#The RSCU value of a codon is calculated in the following way:
# -the number of codons associated with a specific amino acid are counted.
# -the total is divided by the degeneracy of that amino acid (i.e. the number of synonymous
# codons which code for that specific amino acid: a number between 1 and 6). Call that total N
# Then, for each codon corresponding to that amino acid, the number of times that codon shows up
# is divided by N. This process gives the RSCU value for every codon.

# Note: if each type of codon corresponding to a specific amino acid is used with equal frequency, then
# the RSCU value for all of those codons will be 1. See Dan Graur's textbook for a more in depth analysis of the
# specific formula.

# Then relative adaptedness for each codon is then calculated. This is found by taking the RSCU value for each
# codon and dividing it by the maximum RSCU value for that amino acid.

# The output of this script is tab-delimited and is in the same format and sequence as what it took in,
# except for instead of the raw codon count, it gives either the RSCU values by codon, or the relative
# adaptedness by codon depending on which output below is selected.


# The two files below are where the RSCU and Relative adaptedness values will
# be printed to
RelativeAdaptednessFile = Path('./RelativeAdaptednessFile.txt')
RSCUFile = Path('./RSCUValues.txt')

# Before new output files are generated, any files that exist from previous runs are removed
# This is done so that old files don't have new entries written to them
if RelativeAdaptednessFile.is_file():
    os.remove('%s'%RelativeAdaptednessFile)
if RSCUFile.is_file():
    os.remove('%s'%RSCUFile)
# The program is given time to remove the files
time.sleep(0.5)

# The name of the input file is declared here
filename = 'RawCodonCount.txt'

with open(filename,'r') as f:
    reader = csv.reader(f,delimiter = '\t')
    for row in reader:
        UID = row[0]
        Virus = row[1]
        GeneDesignation = row[2]
        n = 3
        while n < 67:
            if row[n] == '':
                row[n] = 0
            row[n] = int(row[n])
            n+=1

        UUU = int(row[3])
        UUC = int(row[4])
        UUA = int(row[5])
        UUG = int(row[6])
        UCU = int(row[7])
        UCC = int(row[8])
        UCA = int(row[9])
        UCG = int(row[10])
        UAU = int(row[11])
        UAC = int(row[12])
        UAA = int(row[13])
        UAG = int(row[14])
        UGU = int(row[15])
        UGC = int(row[16])
        UGA = int(row[17])
        UGG = int(row[18])
        CUU = int(row[19])
        CUC = int(row[20])
        CUA = int(row[21])
        CUG = int(row[22])
        CCU = int(row[23])
        CCC = int(row[24])
        CCA = int(row[25])
        CCG = int(row[26])
        CAU = int(row[27])
        CAC = int(row[28])
        CAA = int(row[29])
        CAG = int(row[30])
        CGU = int(row[31])
        CGC = int(row[32])
        CGA = int(row[33])
        CGG = int(row[34])
        AUU = int(row[35])
        AUC = int(row[36])
        AUA = int(row[37])
        AUG = int(row[38])
        ACU = int(row[39])
        ACC = int(row[40])
        ACA = int(row[41])
        ACG = int(row[42])
        AAU = int(row[43])
        AAC = int(row[44])
        AAA = int(row[45])
        AAG = int(row[46])
        AGU = int(row[47])
        AGC = int(row[48])
        AGA = int(row[49])
        AGG = int(row[50])
        GUU = int(row[51])
        GUC = int(row[52])
        GUA = int(row[53])
        GUG = int(row[54])
        GCU = int(row[55])
        GCC = int(row[56])
        GCA = int(row[57])
        GCG = int(row[58])
        GAU = int(row[59])
        GAC = int(row[60])
        GAA = int(row[61])
        GAG = int(row[62])
        GGU = int(row[63])
        GGC = int(row[64])
        GGA = int(row[65])
        GGG = int(row[66])

        CodonIndex = [[0,'TTT','TTC'],[1,'TTA','TTG','CTT','CTC','CTA','CTG'],[2,'ATT','ATC','ATA'],[3,'ATG'],[4,'GTT','GTC','GTA','GTG'],[5,'TCT','TCC','TCA','TCG'],[6,'CCT','CCC','CCA','CCG'],[7,'ACT','ACC','ACA','ACG'],[8,'GCT','GCC','GCA','GCG'],[9,'TAT','TAC'],[10,'TAA','TAG','TGA'],[11,'CAT','CAC'],[12,'CAA','CAG'],[13,'AAT','AAC'],[14,'AAA','AAG'],[15,'GAT','GAC'],[16,'GAA','GAG'],[17,'TGT','TGC'],[18,'TGG'],[19,'CGT','CGC','CGA','CGG','AGA','AGG'],[20,'AGT','AGC'],[21,'GGT','GGC','GGA','GGG']]



        CodonCount = [[0,UUU,UUC],[1,UUA,UUG,CUU,CUC,CUA,CUG],[2,AUU,AUC,AUA],[3,AUG],[4,GUU,GUC,GUA,GUG],[5,UCU,UCC,UCA,UCG],[6,CCU,CCC,CCA,CCG],[7,ACU,ACC,ACA,ACG],[8,GCU,GCC,GCA,GCG],[9,UAU,UAC],[10,UAA,UAG,UGA],[11,CAU,CAC],[12,CAA,CAG],[13,AAU,AAC],[14,AAA,AAG],[15,GAU,GAC],[16,GAA,GAG],[17,UGU,UGC],[18,UGG],[19,CGU,CGC,CGA,CGG,AGA,AGG],[20,AGU,AGC],[21,GGU,GGC,GGA,GGG]]

        RSCUList =[[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]

        RelativeAdaptedness =[[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]


        for list in CodonCount:
            index = list[0]
            newlist = list[1:]
            codoncount = sum(newlist)
            n = len(newlist)
            for entry in newlist :
                if codoncount == 0:
                    RSCUValue = 0
                else:
                    RSCUValue = entry/((1/n)*codoncount)
                RSCUList[index][newlist.index(entry)+1] = RSCUValue

        RSCU_UUU = format(RSCUList[0][1],'.3f')
        RSCU_UUC = format(RSCUList[0][2],'.3f')
        RSCU_UUA = format(RSCUList[1][1],'.3f')
        RSCU_UUG = format(RSCUList[1][2],'.3f')
        RSCU_UCU = format(RSCUList[5][1],'.3f')
        RSCU_UCC = format(RSCUList[5][2],'.3f')
        RSCU_UCA = format(RSCUList[5][3],'.3f')
        RSCU_UCG = format(RSCUList[5][4],'.3f')
        RSCU_UAU = format(RSCUList[9][1],'.3f')
        RSCU_UAC = format(RSCUList[9][2],'.3f')
        RSCU_UAA = format(RSCUList[10][1],'.3f')
        RSCU_UAG = format(RSCUList[10][2],'.3f')
        RSCU_UGU = format(RSCUList[17][1],'.3f')
        RSCU_UGC = format(RSCUList[17][2],'.3f')
        RSCU_UGA = format(RSCUList[10][3],'.3f')
        RSCU_UGG = format(RSCUList[18][1],'.3f')
        RSCU_CUU = format(RSCUList[1][3],'.3f')
        RSCU_CUC = format(RSCUList[1][4],'.3f')
        RSCU_CUA = format(RSCUList[1][5],'.3f')
        RSCU_CUG = format(RSCUList[1][6],'.3f')
        RSCU_CCU = format(RSCUList[6][1],'.3f')
        RSCU_CCC = format(RSCUList[6][2],'.3f')
        RSCU_CCA = format(RSCUList[6][3],'.3f')
        RSCU_CCG = format(RSCUList[6][4],'.3f')
        RSCU_CAU = format(RSCUList[11][1],'.3f')
        RSCU_CAC = format(RSCUList[11][2],'.3f')
        RSCU_CAA = format(RSCUList[12][1],'.3f')
        RSCU_CAG = format(RSCUList[12][2],'.3f')
        RSCU_CGU = format(RSCUList[19][1],'.3f')
        RSCU_CGC = format(RSCUList[19][2],'.3f')
        RSCU_CGA = format(RSCUList[19][3],'.3f')
        RSCU_CGG = format(RSCUList[19][4],'.3f')
        RSCU_AUU = format(RSCUList[2][1],'.3f')
        RSCU_AUC = format(RSCUList[2][2],'.3f')
        RSCU_AUA = format(RSCUList[2][3],'.3f')
        RSCU_AUG = format(RSCUList[3][1],'.3f')
        RSCU_ACU = format(RSCUList[7][1],'.3f')
        RSCU_ACC = format(RSCUList[7][2],'.3f')
        RSCU_ACA = format(RSCUList[7][3],'.3f')
        RSCU_ACG = format(RSCUList[7][4],'.3f')
        RSCU_AAU = format(RSCUList[13][1],'.3f')
        RSCU_AAC = format(RSCUList[13][2],'.3f')
        RSCU_AAA = format(RSCUList[14][1],'.3f')
        RSCU_AAG = format(RSCUList[14][2],'.3f')
        RSCU_AGU = format(RSCUList[20][1],'.3f')
        RSCU_AGC = format(RSCUList[20][2],'.3f')
        RSCU_AGA = format(RSCUList[19][5],'.3f')
        RSCU_AGG = format(RSCUList[19][6],'.3f')
        RSCU_GUU = format(RSCUList[4][1],'.3f')
        RSCU_GUC = format(RSCUList[4][2],'.3f')
        RSCU_GUA = format(RSCUList[4][3],'.3f')
        RSCU_GUG = format(RSCUList[4][4],'.3f')
        RSCU_GCU = format(RSCUList[8][1],'.3f')
        RSCU_GCC = format(RSCUList[8][2],'.3f')
        RSCU_GCA = format(RSCUList[8][3],'.3f')
        RSCU_GCG = format(RSCUList[8][4],'.3f')
        RSCU_GAU = format(RSCUList[15][1],'.3f')
        RSCU_GAC = format(RSCUList[15][2],'.3f')
        RSCU_GAA = format(RSCUList[16][1],'.3f')
        RSCU_GAG = format(RSCUList[16][2],'.3f')
        RSCU_GGU = format(RSCUList[21][1],'.3f')
        RSCU_GGC = format(RSCUList[21][2],'.3f')
        RSCU_GGA = format(RSCUList[21][3],'.3f')
        RSCU_GGG = format(RSCUList[21][4],'.3f')
        
        output = [UID,Virus,RSCU_UUU,RSCU_UUC,RSCU_UUA,RSCU_UUG,RSCU_UCU,RSCU_UCC,RSCU_UCA,RSCU_UCG,RSCU_UAU,RSCU_UAC,RSCU_UAA,RSCU_UAG,RSCU_UGU,RSCU_UGC,RSCU_UGA,RSCU_UGG,RSCU_CUU,RSCU_CUC,RSCU_CUA,RSCU_CUG,RSCU_CCU,RSCU_CCC,RSCU_CCA,RSCU_CCG,RSCU_CAU,RSCU_CAC,RSCU_CAA,RSCU_CAG,RSCU_CGU,RSCU_CGC,RSCU_CGA,RSCU_CGG,RSCU_AUU,RSCU_AUC,RSCU_AUA,RSCU_AUG,RSCU_ACU,RSCU_ACC,RSCU_ACA,RSCU_ACG,RSCU_AAU,RSCU_AAC,RSCU_AAA,RSCU_AAG,RSCU_AGU,RSCU_AGC,RSCU_AGA,RSCU_AGG,RSCU_GUU,RSCU_GUC,RSCU_GUA,RSCU_GUG,RSCU_GCU,RSCU_GCC,RSCU_GCA,RSCU_GCG,RSCU_GAU,RSCU_GAC ,RSCU_GAA,RSCU_GAG ,RSCU_GGU,RSCU_GGC,RSCU_GGA,RSCU_GGG]
        RSCUOutputFile = open('%s'%RSCUFile ,'a')
        RSCUOutputFile.write('\t'.join(output[0:]) + '\n')
        RSCUOutputFile.close()



        for list in RSCUList:
            index = list[0]
            newlist = list[1:]
            maxvalue = max(newlist)
            for item in newlist:
                if maxvalue == 0:
                    RelAdapt = 1
                else:
                    RelAdapt = item/maxvalue
                RelativeAdaptedness[index][newlist.index(item)+1] = RelAdapt

        RelAdapt_UUU = format(RelativeAdaptedness[0][1],'.3f')
        RelAdapt_UUC = format(RelativeAdaptedness[0][2],'.3f')
        RelAdapt_UUA = format(RelativeAdaptedness[1][1],'.3f')
        RelAdapt_UUG = format(RelativeAdaptedness[1][2],'.3f')
        RelAdapt_UCU = format(RelativeAdaptedness[5][1],'.3f')
        RelAdapt_UCC = format(RelativeAdaptedness[5][2],'.3f')
        RelAdapt_UCA = format(RelativeAdaptedness[5][3],'.3f')
        RelAdapt_UCG = format(RelativeAdaptedness[5][4],'.3f')
        RelAdapt_UAU = format(RelativeAdaptedness[9][1],'.3f')
        RelAdapt_UAC = format(RelativeAdaptedness[9][2],'.3f')
        RelAdapt_UAA = format(RelativeAdaptedness[10][1],'.3f')
        RelAdapt_UAG = format(RelativeAdaptedness[10][2],'.3f')
        RelAdapt_UGU = format(RelativeAdaptedness[17][1],'.3f')
        RelAdapt_UGC = format(RelativeAdaptedness[17][2],'.3f')
        RelAdapt_UGA = format(RelativeAdaptedness[10][3],'.3f')
        RelAdapt_UGG = format(RelativeAdaptedness[18][1],'.3f')
        RelAdapt_CUU = format(RelativeAdaptedness[1][3],'.3f')
        RelAdapt_CUC = format(RelativeAdaptedness[1][4],'.3f')
        RelAdapt_CUA = format(RelativeAdaptedness[1][5],'.3f')
        RelAdapt_CUG = format(RelativeAdaptedness[1][6],'.3f')
        RelAdapt_CCU = format(RelativeAdaptedness[6][1],'.3f')
        RelAdapt_CCC = format(RelativeAdaptedness[6][2],'.3f')
        RelAdapt_CCA = format(RelativeAdaptedness[6][3],'.3f')
        RelAdapt_CCG = format(RelativeAdaptedness[6][4],'.3f')
        RelAdapt_CAU = format(RelativeAdaptedness[11][1],'.3f')
        RelAdapt_CAC = format(RelativeAdaptedness[11][2],'.3f')
        RelAdapt_CAA = format(RelativeAdaptedness[12][1],'.3f')
        RelAdapt_CAG = format(RelativeAdaptedness[12][2],'.3f')
        RelAdapt_CGU = format(RelativeAdaptedness[19][1],'.3f')
        RelAdapt_CGC = format(RelativeAdaptedness[19][2],'.3f')
        RelAdapt_CGA = format(RelativeAdaptedness[19][3],'.3f')
        RelAdapt_CGG = format(RelativeAdaptedness[19][4],'.3f')
        RelAdapt_AUU = format(RelativeAdaptedness[2][1],'.3f')
        RelAdapt_AUC = format(RelativeAdaptedness[2][2],'.3f')
        RelAdapt_AUA = format(RelativeAdaptedness[2][3],'.3f')
        RelAdapt_AUG = format(RelativeAdaptedness[3][1],'.3f')
        RelAdapt_ACU = format(RelativeAdaptedness[7][1],'.3f')
        RelAdapt_ACC = format(RelativeAdaptedness[7][2],'.3f')
        RelAdapt_ACA = format(RelativeAdaptedness[7][3],'.3f')
        RelAdapt_ACG = format(RelativeAdaptedness[7][4],'.3f')
        RelAdapt_AAU = format(RelativeAdaptedness[13][1],'.3f')
        RelAdapt_AAC = format(RelativeAdaptedness[13][2],'.3f')
        RelAdapt_AAA = format(RelativeAdaptedness[14][1],'.3f')
        RelAdapt_AAG = format(RelativeAdaptedness[14][2],'.3f')
        RelAdapt_AGU = format(RelativeAdaptedness[20][1],'.3f')
        RelAdapt_AGC = format(RelativeAdaptedness[20][2],'.3f')
        RelAdapt_AGA = format(RelativeAdaptedness[19][5],'.3f')
        RelAdapt_AGG = format(RelativeAdaptedness[19][6],'.3f')
        RelAdapt_GUU = format(RelativeAdaptedness[4][1],'.3f')
        RelAdapt_GUC = format(RelativeAdaptedness[4][2],'.3f')
        RelAdapt_GUA = format(RelativeAdaptedness[4][3],'.3f')
        RelAdapt_GUG = format(RelativeAdaptedness[4][4],'.3f')
        RelAdapt_GCU = format(RelativeAdaptedness[8][1],'.3f')
        RelAdapt_GCC = format(RelativeAdaptedness[8][2],'.3f')
        RelAdapt_GCA = format(RelativeAdaptedness[8][3],'.3f')
        RelAdapt_GCG = format(RelativeAdaptedness[8][4],'.3f')
        RelAdapt_GAU = format(RelativeAdaptedness[15][1],'.3f')
        RelAdapt_GAC = format(RelativeAdaptedness[15][2],'.3f')
        RelAdapt_GAA = format(RelativeAdaptedness[16][1],'.3f')
        RelAdapt_GAG = format(RelativeAdaptedness[16][2],'.3f')
        RelAdapt_GGU = format(RelativeAdaptedness[21][1],'.3f')
        RelAdapt_GGC = format(RelativeAdaptedness[21][2],'.3f')
        RelAdapt_GGA = format(RelativeAdaptedness[21][3],'.3f')
        RelAdapt_GGG = format(RelativeAdaptedness[21][4],'.3f')

        output = [UID,Virus,RelAdapt_UUU,RelAdapt_UUC,RelAdapt_UUA,RelAdapt_UUG,RelAdapt_UCU,RelAdapt_UCC,RelAdapt_UCA,RelAdapt_UCG,RelAdapt_UAU,RelAdapt_UAC,RelAdapt_UAA,RelAdapt_UAG,RelAdapt_UGU,RelAdapt_UGC,RelAdapt_UGA,RelAdapt_UGG,RelAdapt_CUU,RelAdapt_CUC,RelAdapt_CUA,RelAdapt_CUG,RelAdapt_CCU,RelAdapt_CCC,RelAdapt_CCA,RelAdapt_CCG,RelAdapt_CAU,RelAdapt_CAC,RelAdapt_CAA,RelAdapt_CAG,RelAdapt_CGU,RelAdapt_CGC,RelAdapt_CGA,RelAdapt_CGG,RelAdapt_AUU,RelAdapt_AUC,RelAdapt_AUA,RelAdapt_AUG,RelAdapt_ACU,RelAdapt_ACC,RelAdapt_ACA,RelAdapt_ACG,RelAdapt_AAU,RelAdapt_AAC,RelAdapt_AAA,RelAdapt_AAG,RelAdapt_AGU,RelAdapt_AGC,RelAdapt_AGA,RelAdapt_AGG,RelAdapt_GUU,RelAdapt_GUC,RelAdapt_GUA,RelAdapt_GUG,RelAdapt_GCU,RelAdapt_GCC,RelAdapt_GCA,RelAdapt_GCG,RelAdapt_GAU,RelAdapt_GAC ,RelAdapt_GAA,RelAdapt_GAG ,RelAdapt_GGU,RelAdapt_GGC,RelAdapt_GGA,RelAdapt_GGG]

        RelativeAdaptednessOutputFile = open('%s'%RelativeAdaptednessFile, 'a')
        RelativeAdaptednessOutputFile.write('\t'.join(output[0:]) + '\n')
        RelativeAdaptednessOutputFile.close()





        
