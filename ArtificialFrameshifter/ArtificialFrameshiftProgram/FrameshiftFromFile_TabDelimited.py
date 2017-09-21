from Bio import SeqIO
from Bio.Seq import Seq
from Remove_Start_Stop_Function import removeStartStop
from Stop_Codon_Remover_Function import stopCodonRemover
from Excess_Nucleotide_Trimmer import frontTrimmer, backTrimmer
import csv
import os, os.path
from pathlib import Path
import time

filename = 'VirusGenes_Nonoverlapping_Controls_Complete.txt'


# The purpose of this program is to take in data from a tab-delimited file,
# and return two different nucleotide sequences
# in tab-delimited format. The two sequences will be the original sequence
# artificially frameshifted into the +1 and +2 frame. The results will
# be printed to two tab-delimited files .

# All shifting and manipulating is done with the stop and start codons
# removed from the body of the sequence. The start and stop codons are
# reattached and the end of the program to each sequence.
# All sequences will have any stop codons appearing in the new reading frame
# removed.


#First the sequences is read from a tab-delimited file.

# OutputList is a list created for each entry and will be where the output
# data will be stored and will eventually be printed from
# UID is the unique ID associated with the entry from the MySQL database where it's stored
# Family is the viral family the gene is from
# Species is the name of the viral species where the gene is from
# GenomeAccession is the ncbi accession number where the gene can be accessed
# GeneDesignation is the control designation for the gene in question. Unshifted genes
# should have a designation of ControlGene. The output frameshifted sequences will be given a
# different GeneDesignation, either ShiftedReadingFrameControlPlusOne or
# ShiftedReadingFrameControlPlusOne depending on the frameshift type
# Gene is the name of the gene
# NucleotideSequence is the input nucleotide sequence which will be used to generate the controls
# ProteinSequence is the translated NucleotideSequence
# NoCysProtein is the translated NucleotideSequence with all cysteines removed

# Format of Output of this file:
# 1 - UID
# 2 - Family
# 3 - Species
# 4 - Genome Accession Number
# 5 - Protein Accession Number
# 6 - Gene Designation
# 7 - Gene Name
# 8 - Nucleotide Sequence
# 9 - Protein Sequence
# 10 - Protein Sequence (No Cys)



# For the purposes of keeping track of entries, a counter is introduced here, and will count the
# number of entries that are successfully manipulated.

# The program starts by checking to see if old results already exist
# If so, they are deleted (this is necessary to prevent writing new entries
# to old files from previous runs)
# The filenames for both of the output files are declared here.
# a while loop instead of an if statement is used to make sure the
# files are actually removed. This is to prevent the program from
# running so quickly that the files are not actually erased from
# the directory leading to old files having new results appended to them
PlusOneFrameshiftedControlsFilename = Path('../PlusOneFrameshiftedControls.txt')
PlusTwoFrameshiftedControlsFilename = Path('../PlusTwoFrameshiftedControls.txt')
while PlusOneFrameshiftedControlsFilename.is_file():
    os.remove('%s'%PlusOneFrameshiftedControlsFilename)

while PlusTwoFrameshiftedControlsFilename.is_file():
    os.remove('%s'%PlusTwoFrameshiftedControlsFilename)

# A number of filters will be implemented to make sure the sequences are:
# 1) Of the correct length (they must be a multiple of 3 in length)
# 2) Are actually protein coding sequences
# 3) Don't contain an N (denoting an unknown section)
# These filters will be printed to stdout at the end of the script to notify
# the user of how many genes were lost to these filters
EntryCounter = 0
NumberWithWrongDesignation = 0
NumberWithNInSequence = 0
NumberWithWrongLength = 0

# The file with the coding genes is opened and the entries
# are given variable names
with open('%s'%filename, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        PlusOneOutputList = []
        PlusTwoOutputList = []
        
        UID = row[0]
        Family = row[1]
        Species = row[2]
        GenomeAccession = row[3]
        GeneAccession = row[4]
        GeneDesignation = row[5]
        Gene = row[6]
        NucleotideSequence = row[7]
        ProteinSequence = row[8]
        NoCysProtein = row[9]
        
        # The output list is set up by adding the data that will not change
        PlusOneOutputList.extend((UID,Family,Species,GenomeAccession,GeneAccession))
        PlusTwoOutputList.extend((UID,Family,Species,GenomeAccession,GeneAccession))

        

        # Filters are implemented to make sure the nucleotide sequences that are CodingGene, and
        # that there are no unknown portions of the sequence and that the length is a multiple of
        # 3. If not, the sequence is ignored
        if GeneDesignation != 'CodingGene':
            NumberWithWrongDesignation += 1
        else:
            if 'N' in str(NucleotideSequence):
                NumberWithNInSequence += 1
            else:
                if len(str(NucleotideSequence))%3 != 0:
                    NumberWithWrongLength +=1
                
                # The nucleotide sequence has its start and stop codons removed using the submodule removeStartStop
                # and the start and stop codons associated with the sequence are stored for later use.
                else:
                    EntryCounter += 1
                    
                    Unaltered_Sequence = str(NucleotideSequence)
                    Trimmed_Sequence = removeStartStop(Unaltered_Sequence)
                    Trimmed_Unshifted_Sequence = str(Trimmed_Sequence[0])
                    start_codon = str(Trimmed_Sequence[1])
                    stop_codon = str(Trimmed_Sequence[2])

                    # Now the manipulated sequences are generated. The +1 frameshift removes the first nucleotide
                    # while the +2 frameshift removes the first and second nucleotide
                    Plus_1_Frameshifted_Sequence  = Trimmed_Unshifted_Sequence[1::]
                    Plus_2_Frameshifted_Sequence = Trimmed_Unshifted_Sequence[2::]


                    # Now the frameshifted sequences are trimmed so they contain a multiple of 3 nucleotides:
                    # This is done using the submodule backTrimmer, which will remove the minimum number of
                    # nucleotides to keep the sequence a multiple of 3 in length
                    Plus_1_Frameshifted_Sequence = str(backTrimmer(Plus_1_Frameshifted_Sequence))
                    Plus_2_Frameshifted_Sequence = str(backTrimmer(Plus_2_Frameshifted_Sequence))


                    # The next step removes all stop codons from the sequences
                    Trimmed_Unshifted_Sequence_Stops_Removed = Trimmed_Unshifted_Sequence
                    Plus_1_Frameshifted_Stops_Removed = str(stopCodonRemover(Plus_1_Frameshifted_Sequence))
                    Plus_2_Frameshifted_Stops_Removed = str(stopCodonRemover(Plus_2_Frameshifted_Sequence))
     

                    # The start and stop codons are now reattached (and if they were non-existent
                    # and replaced with '***', those start/stop codons will be removed
                    Trimmed_Unshifted_Sequence = start_codon + Trimmed_Unshifted_Sequence_Stops_Removed + stop_codon
                    Trimmed_Unshifted_Sequence = Trimmed_Unshifted_Sequence.replace("*","")
                    
                    Plus_1_Frameshifted = start_codon + Plus_1_Frameshifted_Stops_Removed + stop_codon
                    Plus_1_Frameshifted = Plus_1_Frameshifted.replace("*","")
        
                    Plus_2_Frameshifted = start_codon + Plus_2_Frameshifted_Stops_Removed + stop_codon
                    Plus_2_Frameshifted = Plus_2_Frameshifted.replace("*","")
        
                    # The cysteines are removed from the resulting translated sequences. This is because
                    # These sequences will be run through the disorder-predicting software IUPred (see Materials and
                    # Methods)
                    Plus1_Protein = Seq(Plus_1_Frameshifted).translate()
                    Plus1_Protein_NoCys = str(Plus1_Protein).replace('C','')
                    
                    Plus2_Protein = Seq(Plus_2_Frameshifted).translate()
                    Plus2_Protein_NoCys = str(Plus2_Protein).replace('C','')
                    
                    # Finally, because of the way the MySQL database was set up, it was easier to only have one control group
                    # printed at a time.
                    # The output is printed in the same format as the input file, only with the nucleotide sequence, protein
                    # sequence and no-cys protein sequences changed, with the GeneDesignation changed to match the control.

                    
                    # The finalized Plus One Frameshifts are printed to an output file
                    PlusOneOutputList.extend(('ShiftedReadingFrameControlPlusOne',Gene,str(Plus_1_Frameshifted),str(Plus1_Protein),str(Plus1_Protein_NoCys)))
                    PlusOneFrameshiftedControls = open('%s'%PlusOneFrameshiftedControlsFilename,'a')
                    PlusOneFrameshiftedControls.write('\t'.join(PlusOneOutputList[0:]) + '\n')
                    PlusOneFrameshiftedControls.close()

                    
                    # The finalized Plus Two Frameshifts are printed to a separate file
                    # both files are stored one directory up from this script.
                    PlusTwoOutputList.extend(('ShiftedReadingFrameControlPlusTwo',Gene,str(Plus_2_Frameshifted),str(Plus2_Protein),str(Plus2_Protein_NoCys)))
                    PlusTwoFrameshiftedControls = open('%s'%PlusTwoFrameshiftedControlsFilename,'a')
                    PlusTwoFrameshiftedControls.write('\t'.join(PlusTwoOutputList[0:]) + '\n')
                    PlusTwoFrameshiftedControls.close()


                    

# The program finally notifies the user how many genes were lost to each of the
# filters and how many genes were successfully manipulated
print('Number of entries successfully manipulated: %s' %EntryCounter)
print('Number lost because they did not have the designation CodingGene: %s' %NumberWithWrongDesignation)
print('Number lost due to presence of "N" (denotes unknown section): %s' %NumberWithNInSequence)
print('Number lost because nucleotide sequence not a multiple of 3: %s' %NumberWithWrongLength)




   
