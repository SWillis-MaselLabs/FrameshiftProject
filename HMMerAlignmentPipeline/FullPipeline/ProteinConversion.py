import Bio
from Bio.Seq import Seq
import csv


def ProteinConversion(OverlappingGenesFilename,ControlsFilename):

    File = open('AllSequencesFile.txt','w')
    with open('%s' %OverlappingGenesFilename, 'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            UID = row[0]
            GeneName = row[1]
            Sequence = Seq(row[2])
            Protein = Sequence.translate()
            output = [str(UID), str(GeneName), str(Protein)]

            File.write('\t'.join(output[0:]) + '\n')


    with open('%s' %ControlsFilename, 'r') as f:
        reader = csv.reader(f,delimiter = '\t')
        for row in reader:
            UID = row[0]
            NewUID = 'c'+UID
            GeneDesignation = row[1]
            GeneName = row[2]
            ProteinSequence = row[3]
            if GeneDesignation != 'ScrambledByNucleotideControl':
                output = [str(NewUID),str(GeneName), str(ProteinSequence)]
                File.write('\t'.join(output[0:]) + '\n')

    File.close()

        
