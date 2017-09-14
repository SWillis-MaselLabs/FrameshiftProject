from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import csv
import sys

def Fasta_File_Maker():

    filename = 'AllSequencesFile.txt'
    with open('%s' %filename, 'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:

            UID = row[0]
            GeneName = row[1]
            Sequence = row[2]
            record=SeqRecord(Seq(Sequence, generic_protein), id = '%s' %UID, name = '%s' %GeneName)
            sys.stdout=open('FastaFile1.txt','a')
            print(record.format('fasta'))
            sys.stdout = open('FastaFile2.txt','a')
            print(record.format('fasta'))
            
