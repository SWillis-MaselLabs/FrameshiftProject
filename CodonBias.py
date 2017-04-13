from Bio.Seq import Seq
import numpy

# The purpose of this script is to read in a nucleotide sequence, determine
# the frequency with which each of the 64 codons show up, calculate the GC
# percent content of the overall sequence, as well as the GC percent content
# of each of each position within each codon.
# The data generated here is printed in the same format as that displayed
# in the online codon usage database found at: http://www.kazusa.or.jp/codon/
# The purpose of the script is to generate data to examine codon bias within
# overlapping virus genes to determine which is the ancestral and which is the
# novel gene.

# Here a nucleotide sequence is defined for analysis
NovelSequence = 'GTTCCACAACGGAATGCAGCTACATTTAATCCGGATGCAGGGTATGTGGCATTTATCAGTAAGTATGGGCAGCAGCTCAACTTTACTGTTGCTAGAGTCTTCTTCCTCAACCAGAAGAAGGCCAAGATGGTCTTACATAAGACGCCACAACCAAGTGTCGATCTTACTTTTGCAGGGGTCAAATTTACAGTGGTTAATAACCATTTTCCCCAGTACACTGCAAATCCAGTGTCAGACACTGCCTTTACGCTTCACCGCATCTCGGGCTACTTA'


# To determine whether a sequence can be analyzed, the script checks whether
# or not the sequence appears with a multiple of three nucelotides. If it does
# not, then a warning is issued to the user and the program exits. 
if len(NovelSequence)%3 != 0:
    print('Error: Sequence Not a Multiple of 3')
    exit
else:
    pass


Sequence1Length = int(len(NovelSequence))
Sequence1CodonNumber = int(Sequence1Length)




count = 0

# A list is defined here which contains all codons by which amino acid they
# code for. The first number in each sublist is used for two purposes:
# 1) For easy reference as to which sublist each codon belongs
# 2) Each subsequent list which is presented is in the same formatting
# as this list, so cross-referencing within the program can be done quickly
# and easily by calling on the first entry of each sublist.

CodonIndex = [[0,'TTT','TTC'],[1,'TTA','TTG','CTT','CTC','CTA','CTG'],[2,'ATT','ATC','ATA'],[3,'ATG'],[4,'GTT','GTC','GTA','GTG'],[5,'TCT','TCC','TCA','TCG'],[6,'CCT','CCC','CCA','CCG'],[7,'ACT','ACC','ACA','ACG'],[8,'GCT','GCC','GCA','GCG'],[9,'TAT','TAC'],[10,'TAA','TAG','TGA'],[11,'CAT','CAC'],[12,'CAA','CAG'],[13,'AAT','AAC'],[14,'AAA','AAG'],[15,'GAT','GAC'],[16,'GAA','GAG'],[17,'TGT','TGC'],[18,'TGG'],[19,'CGT','CGC','CGA','CGG','AGA','AGG'],[20,'AGT','AGC'],[21,'GGT','GGC','GGA','GGG']]

# The codon count contains the same reference numbers as the list above, and
# specififies the count of all codons equal to zero. 

CodonCount = [[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]

# This list is in the same formatting as the previous two lists, and
# each reference number of each sublist above corresponds to the amino acid
# list position below

Index = ['Phe','Leu','Ile','Val','Ser','Pro','Thr','Ala','Tyr','Stop','His','Gln','Asn','Lys','Asp','Glu','Cys','Trp','Arg','Ser','Gly']

# The total codon count is specified as zero and will be added to with each
# loop through the nucleotide sequence.
TotalCodonCount = 0

# The following loop counts each type of codon through a simple series of for
# loops. 
for n in range(0,Sequence1CodonNumber,3):
    TotalCodonCount += 1
    codon = NovelSequence[n:n+3]
    for entry in CodonIndex:
        if codon in entry:
            external_index = CodonIndex.index(entry)
            internal_index = entry.index(codon)
            CodonCount[external_index][internal_index]=CodonCount[external_index][internal_index]+1

# Once the prevelence of each codon has been determined, the GC percent content
# is determined

GCContent = [0,0,0]
n=0
for nucleotide in NovelSequence:
    CodonPosition = n%3
    if nucleotide == 'G':
        GCContent[CodonPosition]+=1
        n += 1
    else:
        if nucleotide == 'C':
            GCContent[CodonPosition] += 1
            n += 1
        else:
            n+= 1
            pass
                      
                      

# I've chosen to print out the data in the same format that is found in the codon degeneracy database



print('Total Number of Codons: %d\n' %TotalCodonCount)

print('Format: [triplet][frequency: per thousand]([number])\n\n')

# Codons UUU,UCU,UAU,UGU
# UUU = CodonCount[0][1]
# UCU = CodonCount[5][1]
# UAU = CodonCount[9][1]
# UGU = CodonCount[17][1]

print('UUU %.1f (%d)\tUCU %.1f (%d)\tUAU %.1f (%d)\tUGU %.1f (%d)' %(numpy.divide(CodonCount[0][1],TotalCodonCount)*1000,CodonCount[0][1],numpy.divide(CodonCount[5][1],TotalCodonCount)*1000,CodonCount[5][1],numpy.divide(CodonCount[9][1],TotalCodonCount)*1000,CodonCount[9][1],numpy.divide(CodonCount[17][1],TotalCodonCount)*1000,CodonCount[17][1]))

# Codons UUC,UCC,UAC,UGC
# UUC = CodonCount[0][2]
# UCC = CodonCount[5][2]
# UAC = CodonCount[9][2]
# UGC = CodonCount[17][2]

print('UUC %.1f (%d)\tUCC %.1f (%d)\tUAC %.1f (%d)\tUGC %.1f (%d)' %(numpy.divide(CodonCount[0][2],TotalCodonCount)*1000,CodonCount[0][2],numpy.divide(CodonCount[5][2],TotalCodonCount)*1000,CodonCount[5][2],numpy.divide(CodonCount[9][2],TotalCodonCount)*1000,CodonCount[9][2],numpy.divide(CodonCount[17][2],TotalCodonCount)*1000,CodonCount[17][2]))

# Codons UUA,UCA,UAA,UGA
# UUA = CodonCount[1][1]
# UCA = CodonCount[5][3]
# UAA = CodonCount[10][1]
# UGA = CodonCount[10][3]

print('UUA %.1f (%d)\tUCA %.1f (%d)\tUAA %.1f (%d)\tUGA %.1f (%d)' %(numpy.divide(CodonCount[1][1],TotalCodonCount)*1000,CodonCount[1][1],numpy.divide(CodonCount[5][3],TotalCodonCount)*1000,CodonCount[5][3],numpy.divide(CodonCount[10][1],TotalCodonCount)*1000,CodonCount[10][1],numpy.divide(CodonCount[10][3],TotalCodonCount)*1000,CodonCount[10][3]))

# Codons UUG,UCG,UAG,UGG
# UUG = CodonCount[1][2]
# UCG = CodonCount[5][4]
# UAG = CodonCount[10][2]
# UGG = CodonCount[18][1]

print('UUG %.1f (%d)\tUCG %.1f (%d)\tUAG %.1f (%d)\tUGG %.1f (%d)\n' %(numpy.divide(CodonCount[1][2],TotalCodonCount)*1000,CodonCount[1][2],numpy.divide(CodonCount[5][4],TotalCodonCount)*1000,CodonCount[5][4],numpy.divide(CodonCount[10][2],TotalCodonCount)*1000,CodonCount[10][2],numpy.divide(CodonCount[18][1],TotalCodonCount)*1000,CodonCount[18][1]))

# Codons CUU,CCU,CAU,CGU
# CUU = CodonCount[1][3]
# CCU = CodonCount[6][1]
# CAU = CodonCount[11][1]
# CGU = CodonCount[19][1]

print('CUU %.1f (%d)\tCCU %.1f (%d)\tCAU %.1f (%d)\tCGU %.1f (%d)' %(numpy.divide(CodonCount[1][3],TotalCodonCount)*1000,CodonCount[1][3],numpy.divide(CodonCount[6][1],TotalCodonCount)*1000,CodonCount[6][1],numpy.divide(CodonCount[11][1],TotalCodonCount)*1000,CodonCount[11][1],numpy.divide(CodonCount[19][1],TotalCodonCount)*1000,CodonCount[19][1]))

# Codons CUC,CCC,CAC,CGC
# CUC = CodonCount[1][4]
# CCC = CodonCount[6][2]
# CAC = CodonCount[11][2]
# CGC = CodonCount[19][2]

print('CUC %.1f (%d)\tCCC %.1f (%d)\tCAC %.1f (%d)\tCGC %.1f (%d)' %(numpy.divide(CodonCount[1][4],TotalCodonCount)*1000,CodonCount[1][4],numpy.divide(CodonCount[6][2],TotalCodonCount)*1000,CodonCount[6][2],numpy.divide(CodonCount[11][2],TotalCodonCount)*1000,CodonCount[11][2],numpy.divide(CodonCount[19][2],TotalCodonCount)*1000,CodonCount[19][2]))

# Codons CUA,CCA,CAA,CGA
# CUA = CodonCount[1][5]
# CCA = CodonCount[6][3]
# CAA = CodonCount[12][1]
# CGA = CodonCount[19][3]

print('CUA %.1f (%d)\tCCA %.1f (%d)\tCAA %.1f (%d)\tCGA %.1f (%d)' %(numpy.divide(CodonCount[1][5],TotalCodonCount)*1000,CodonCount[1][5],numpy.divide(CodonCount[6][3],TotalCodonCount)*1000,CodonCount[6][3],numpy.divide(CodonCount[12][1],TotalCodonCount)*1000,CodonCount[12][1],numpy.divide(CodonCount[19][3],TotalCodonCount)*1000,CodonCount[19][3]))

# Codons CUG,CCG,CAG,CGG
# CUG = CodonCount[1][6]
# CCG = CodonCount[6][4]
# CAG = CodonCount[12][2]
# CGG = CodonCount[19][4]

print('CUG %.1f (%d)\tCCG %.1f (%d)\tCAG %.1f (%d)\tCGG %.1f (%d)\n' %(numpy.divide(CodonCount[1][6],TotalCodonCount)*1000,CodonCount[1][6],numpy.divide(CodonCount[6][4],TotalCodonCount)*1000,CodonCount[6][4],numpy.divide(CodonCount[12][2],TotalCodonCount)*1000,CodonCount[12][2],numpy.divide(CodonCount[19][4],TotalCodonCount)*1000,CodonCount[19][4]))

# Codons AUU,ACU,AAU,AGU
# AUU = CodonCount[2][1]
# ACU = CodonCount[7][1]
# AAU = CodonCount[13][1]
# AGU = CodonCount[20][1]

print('AUU %.1f (%d)\tACU %.1f (%d)\tAAU %.1f (%d)\tAGU %.1f (%d)' %(numpy.divide(CodonCount[2][1],TotalCodonCount)*1000,CodonCount[2][1],numpy.divide(CodonCount[7][1],TotalCodonCount)*1000,CodonCount[7][1],numpy.divide(CodonCount[13][1],TotalCodonCount)*1000,CodonCount[13][1],numpy.divide(CodonCount[20][1],TotalCodonCount)*1000,CodonCount[20][1]))

# Codons AUC,ACC,AAC,AGC
# AUC = CodonCount[2][2]
# ACC = CodonCount[7][2]
# AAC = CodonCount[13][2]
# AGC = CodonCount[20][2]


print('AUC %.1f (%d)\tACC %.1f (%d)\tAAC %.1f (%d)\tAGC %.1f (%d)' %(numpy.divide(CodonCount[2][2],TotalCodonCount)*1000,CodonCount[2][2],numpy.divide(CodonCount[7][2],TotalCodonCount)*1000,CodonCount[7][2],numpy.divide(CodonCount[13][2],TotalCodonCount)*1000,CodonCount[13][2],numpy.divide(CodonCount[20][2],TotalCodonCount)*1000,CodonCount[20][2]))

# Codons AUA,ACA,AAA,AGA
# AUA = CodonCount[2][3]
# ACA = CodonCount[7][3]
# AAA = CodonCount[14][1]
# AGA = CodonCount[19][5]

print('AUA %.1f (%d)\tACA %.1f (%d)\tAAA %.1f (%d)\tAGA %.1f (%d)' %(numpy.divide(CodonCount[2][3],TotalCodonCount)*1000,CodonCount[2][3],numpy.divide(CodonCount[7][3],TotalCodonCount)*1000,CodonCount[7][3],numpy.divide(CodonCount[14][1],TotalCodonCount)*1000,CodonCount[14][1],numpy.divide(CodonCount[19][5],TotalCodonCount)*1000,CodonCount[19][5]))

# Codons AUG,ACG,AAG,AGG
# AUG = CodonCount[3][1]
# ACG = CodonCount[7][4]
# AAG = CodonCount[14][2]
# AGG = CodonCount[19][6]

print('AUG %.1f (%d)\tACG %.1f (%d)\tAAG %.1f (%d)\tAGG %.1f (%d)\n' %(numpy.divide(CodonCount[3][1],TotalCodonCount)*1000,CodonCount[3][1],numpy.divide(CodonCount[7][4],TotalCodonCount)*1000,CodonCount[7][4],numpy.divide(CodonCount[14][2],TotalCodonCount)*1000,CodonCount[14][2],numpy.divide(CodonCount[19][6],TotalCodonCount)*1000,CodonCount[19][6]))

# Codons GUU,GCU,GAU,GGU
# GUU = CodonCount[4][1]
# GCU = CodonCount[8][1]
# GAU = CodonCount[15][1]
# GGU = CodonCount[21][1]

print('GUU %.1f (%d)\tGCU %.1f (%d)\tGAU %.1f (%d)\tGGU %.1f (%d)' %(numpy.divide(CodonCount[4][1],TotalCodonCount)*1000,CodonCount[4][1],numpy.divide(CodonCount[8][1],TotalCodonCount)*1000,CodonCount[8][1],numpy.divide(CodonCount[15][1],TotalCodonCount)*1000,CodonCount[15][1],numpy.divide(CodonCount[21][1],TotalCodonCount)*1000,CodonCount[21][1]))

# Codons GUC,GCC,GAC,GGC
# GUC = CodonCount[4][2]
# GCC = CodonCount[8][2]
# GAC = CodonCount[15][2]
# GGC = CodonCount[21][2]

print('GUC %.1f (%d)\tGCC %.1f (%d)\tGAC %.1f (%d)\tGGC %.1f (%d)' %(numpy.divide(CodonCount[4][2],TotalCodonCount)*1000,CodonCount[4][2],numpy.divide(CodonCount[8][2],TotalCodonCount)*1000,CodonCount[8][2],numpy.divide(CodonCount[15][2],TotalCodonCount)*1000,CodonCount[15][2],numpy.divide(CodonCount[21][2],TotalCodonCount)*1000,CodonCount[21][2]))

# Codons GUA,GCA,GAA,GGA
# GUA = CodonCount[4][3]
# GCA = CodonCount[8][3]
# GAA = CodonCount[16][1]
# GGA = CodonCount[21][3]

print('GUA %.1f (%d)\tGCA %.1f (%d)\tGAA %.1f (%d)\tGGA %.1f (%d)' %(numpy.divide(CodonCount[4][3],TotalCodonCount)*1000,CodonCount[4][3],numpy.divide(CodonCount[8][3],TotalCodonCount)*1000,CodonCount[8][3],numpy.divide(CodonCount[16][1],TotalCodonCount)*1000,CodonCount[16][1],numpy.divide(CodonCount[21][3],TotalCodonCount)*1000,CodonCount[21][3]))

#Codons GUG,GCG,GAG,GGG
# GUG = CodonCount[4][4]
# GCG = CodonCount[8][4]
# GAG = CodonCount[16][2]
# GGG = CodonCount[21][4]

print('GUG %.1f (%d)\tGCG %.1f (%d)\tGAG %.1f (%d)\tGGG %.1f (%d)\n\n' %(numpy.divide(CodonCount[4][4],TotalCodonCount)*1000,CodonCount[4][4],numpy.divide(CodonCount[8][4],TotalCodonCount)*1000,CodonCount[8][4],numpy.divide(CodonCount[16][2],TotalCodonCount)*1000,CodonCount[16][2],numpy.divide(CodonCount[21][4],TotalCodonCount)*1000,CodonCount[21][4]))

FirstLetter = GCContent[0]
SecondLetter = GCContent[1]
ThirdLetter = GCContent[2]
CodingGC = FirstLetter+SecondLetter+ThirdLetter

percentGC = numpy.divide(CodingGC,len(NovelSequence))*100

FirstLetterPercent = numpy.divide(FirstLetter,TotalCodonCount)*100
SecondLetterPercent = numpy.divide(SecondLetter,TotalCodonCount)*100
ThirdLetterPercent = numpy.divide(ThirdLetter,TotalCodonCount)*100

print('_____________________________________________________________\n')

print('Coding GC %.2f \t 1st letter GC %.2f \t 2nd letter GC %.2f \t 3d letter GC %.2f' %(percentGC,FirstLetterPercent,SecondLetterPercent,ThirdLetterPercent))

# Because in some cases there are very few codons in a strand (<100), it can be difficult to do
# straight comparissons. In order to facilitate codon bias analysis, I have included a count of
# the highly degenerate codons which code for Leucine and Arginine.

Leu = CodonCount[1][1:]
LeuPercent = []
Arg = CodonCount[19][1:]
ArgPercent = []
TotalLeuCount = 0
TotalArgCount = 0

for n in Leu:
    number = int(n)
    TotalLeuCount += number

for n in Arg:
    number = int(n)
    TotalArgCount += number

for n in Leu:
    number = int(n)
    percentLeu = numpy.divide(number,TotalLeuCount)*100
    percentLeu = format(percentLeu,'.2f')
    LeuPercent.append(percentLeu)

for n in Arg:
    number = int(n)
    percentArg = format(numpy.divide(number,TotalArgCount)*100,'.2f')
    ArgPercent.append(percentArg)


print('\nLeucine')
print('Leu Codons: %s' %str(CodonIndex[1][1:]).replace('T','U'))
print('Leu: %s' %CodonCount[1][1:])
print('Leu Percentage: %s' %LeuPercent)
print('\nArginine')
print('Arg Codons: %s' %str(CodonIndex[19][1:]).replace('T','U'))
print('Arg: %s' %CodonCount[19][1:])
print('Arg Percentage: %s' %ArgPercent)
