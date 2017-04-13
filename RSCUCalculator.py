from Bio.Seq import Seq
import numpy
import csv


with open('GeneSegmentsForCodonAnalysis.txt','r') as f:
    reader = csv.reader(f,delimiter = '\t')
    for row in reader:
        UID = row[0]
        Species = row[1]
        GeneName = row[2]
        Verification = row[3]
        Sequence = row[4]

        file = open("%s_%s.txt" %(Species,GeneName),"w")

        if len(Sequence)%3 != 0:
            print('Error: Sequence Not a Multiple of 3')
            exit
        
        else:
            pass


        Sequence1Length = int(len(Sequence))
        Sequence1CodonNumber = int(Sequence1Length)




        count = 0


        CodonIndex = [[0,'TTT','TTC'],[1,'TTA','TTG','CTT','CTC','CTA','CTG'],[2,'ATT','ATC','ATA'],[3,'ATG'],[4,'GTT','GTC','GTA','GTG'],[5,'TCT','TCC','TCA','TCG'],[6,'CCT','CCC','CCA','CCG'],[7,'ACT','ACC','ACA','ACG'],[8,'GCT','GCC','GCA','GCG'],[9,'TAT','TAC'],[10,'TAA','TAG','TGA'],[11,'CAT','CAC'],[12,'CAA','CAG'],[13,'AAT','AAC'],[14,'AAA','AAG'],[15,'GAT','GAC'],[16,'GAA','GAG'],[17,'TGT','TGC'],[18,'TGG'],[19,'CGT','CGC','CGA','CGG','AGA','AGG'],[20,'AGT','AGC'],[21,'GGT','GGC','GGA','GGG']]



        CodonCount = [[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]

        RSCUList =[[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]

        RelativeAdaptedness =[[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]

        GeneRelativeAdaptedness = []

        Index = ['Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr','Stop','His','Gln','Asn','Lys','Asp','Glu','Cys','Trp','Arg','Ser','Gly']

        TotalCodonCount = 0

        for n in range(0,Sequence1CodonNumber,3):
            TotalCodonCount += 1
            codon = Sequence[n:n+3]
            for entry in CodonIndex:
                if codon in entry:
                    external_index = CodonIndex.index(entry)
                    internal_index = entry.index(codon)
                    CodonCount[external_index][internal_index]=CodonCount[external_index][internal_index]+1

        file.write("Total Codon Count: %d\n\n" %TotalCodonCount)
        #print('AminoAcidProduct\tCodon\tRSCU')
        for count in CodonCount:
            AminoAcidReference = count[0]
            ListOfCount = count[1::]
            total = sum(ListOfCount)
            NumberOfSynonymous = len(ListOfCount)
            n=1
            for codon in ListOfCount:
                codon_index = n
                if total == 0:
                    RSCU = 0
                else:
                    RSCU = codon/((1/NumberOfSynonymous)*total)
        
                RSCUList[AminoAcidReference][codon_index] = format(RSCU,'.3f')
                n+=1

        for entry in RSCUList:
            AminoAcidReference = entry[0]
            entry = entry[1::]
            maximum = float(max(entry))
            m=1
            for item in entry:
                codon_index = m
                if maximum == 0:
                    w_i = 0
                else:
                    w_i = float(item)/maximum
                RelativeAdaptedness[AminoAcidReference][codon_index]=w_i
                GeneRelativeAdaptedness.append(w_i)
                m+=1
        file.write('AminoAcid\tCodon\tRSCU\tRelativeAdaptedness\n')
        for entry in CodonIndex:
            AminoAcidReference = entry[0]
            entry = entry[1::]
            m = 1
            for item in entry:
                codon_index = m
                codon = entry[m-1]
                AminoAcid = Index[AminoAcidReference]
                RSCU_codon = RSCUList[AminoAcidReference][codon_index]
                RelativeAdaptedness_codon = RelativeAdaptedness[AminoAcidReference][codon_index]
                file.write('%s\t%s\t%s\t%s\n' %(AminoAcid,codon,RSCU_codon,RelativeAdaptedness_codon))
        file.close()
 

        
