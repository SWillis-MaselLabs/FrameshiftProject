from Bio import SeqIO
from Bio.Seq import Seq
from Remove_Start_Stop_Function import removeStartStop
from Stop_Codon_Remover_Function import stopCodonRemover
from Excess_Nucleotide_Trimmer import frontTrimmer, backTrimmer
import csv

# The purpose of this program is to take in data from a tab-delimited file,
#read in the genetic sequences,and return six different nucleotide sequences
# in tab-delimited fomat. The six sequences will be:
# original sequence, +1 frameshift, +2 frameshift, reverse complement,
# +1 reverse complement frameshift, +2 reverse complement frameshift.
# All shifting and manipulating is done with the stop and start codons
# removed from the body of the sequence. The start and stop codons are
# reattached and the end of the program to each sequence.
# All sequences will have any stop codons removed from their main bodies,
# except the main sequence (which should not have any stop codons in it)


#First we need to read in sequences to work with from a tab-delimited file.

with open('GallusGallus_NoDuplicates.tab', 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:

        # Because I will need to use this program with a few different files
        # with different formats (all tab-delimited, though), the various
        # entries will be defined here for ease of use, and so if the program
        # needs to be modified, the indices will only need to be changed here
        
        EnsemblUID = row[0]
        TranscriptUID = row[1]
        NucleotideSequence = row[2]
        ProteinSequence = row[3]
        NoCysProteinSequence = row[4]
        GeneDesignation = row[5]
        
        if GeneDesignation != 'CodingGene':
            pass
        else:
            if 'N' in str(NucleotideSequence):
                pass
            else:
                if len(str(NucleotideSequence))%3 != 0:
                    pass
                
                else:
                    Unaltered_Sequence = str(NucleotideSequence)
                    Trimmed_Sequence = removeStartStop(Unaltered_Sequence)
                    Trimmed_Unshifted_Sequence = str(Trimmed_Sequence[0])
                    start_codon = str(Trimmed_Sequence[1])
                    stop_codon = str(Trimmed_Sequence[2])

                    # Now we can generate the strings we will need to work with. First, we produce the
                    # +1, +2, -1, -2 frameshifted sequences
                    
                    Plus_1_Frameshifted_Sequence  = Trimmed_Unshifted_Sequence[1::]
                    Plus_2_Frameshifted_Sequence = Trimmed_Unshifted_Sequence[2::]
                    Minus_1_Frameshifted_Sequence = Trimmed_Unshifted_Sequence[: len(Trimmed_Unshifted_Sequence)-1 :]
                    Minus_2_Frameshifted_Sequence = Trimmed_Unshifted_Sequence[: len(Trimmed_Unshifted_Sequence)-2 :]

                    # Now the frameshifted sequences will be trimmed so they contain a multiple of 3 nucleotides:
                    # We will use the frontTrimmer on the -1 and -2 frameshifted sequences as this will translate
                    # to being the same thing as applying the backTrimmer to the +1, +2 framshifted reverse
                    # complement sequences
                    # Also, while the input sequence SHOULD be a multiple of 3, I will run it through the
                    # backTrimmer as well.

                    Trimmed_Unshifted_Sequence = str(backTrimmer(Trimmed_Unshifted_Sequence))
                    Plus_1_Frameshifted_Sequence = str(backTrimmer(Plus_1_Frameshifted_Sequence))
                    Plus_2_Frameshifted_Sequence = str(backTrimmer(Plus_2_Frameshifted_Sequence))
                    Minus_1_Frameshifted_Sequence = str(frontTrimmer(Minus_1_Frameshifted_Sequence))
                    Minus_2_Frameshifted_Sequence = str(frontTrimmer(Minus_2_Frameshifted_Sequence))
        

                    # Here the reverse complement and the +1, +2 frameshifted sequences of the
                    # The reverse complement will be generated as strings:

                    Reverse_Complement = str(Seq(Trimmed_Unshifted_Sequence).reverse_complement())
                    Plus_1_Frameshifted_Reverse_Complement = str(Seq(Minus_1_Frameshifted_Sequence).reverse_complement())
                    Plus_2_Frameshifted_Reverse_Complement = str(Seq(Minus_2_Frameshifted_Sequence).reverse_complement())
        

                    # The next step removes the stop codons from all of our sequences

                    Trimmed_Unshifted_Sequence_Stops_Removed = Trimmed_Unshifted_Sequence
                    Plus_1_Frameshifted_Stops_Removed = str(stopCodonRemover(Plus_1_Frameshifted_Sequence))
                    Plus_2_Frameshifted_Stops_Removed = str(stopCodonRemover(Plus_2_Frameshifted_Sequence))
                    Reverse_Complement_Stops_Removed = str(stopCodonRemover(Reverse_Complement))
                    Plus_1_Frameshifted_Reverse_Complement_Stops_Removed = str(stopCodonRemover(Plus_1_Frameshifted_Reverse_Complement))
                    Plus_2_Frameshifted_Reverse_Complement_Stops_Removed = str(stopCodonRemover(Plus_2_Frameshifted_Reverse_Complement))
        

                    # The start and stop codons will now be reattached (and if they were non-existent
                    # and replaced with '***', those start/stop codons will be removed

                    Trimmed_Unshifted_Sequence = start_codon + Trimmed_Unshifted_Sequence_Stops_Removed + stop_codon
                    Trimmed_Unshifted_Sequence = Trimmed_Unshifted_Sequence.replace("*","")
                    
                    Plus_1_Frameshifted = start_codon + Plus_1_Frameshifted_Stops_Removed + stop_codon
                    Plus_1_Frameshifted = Plus_1_Frameshifted.replace("*","")
        
                    Plus_2_Frameshifted = start_codon + Plus_2_Frameshifted_Stops_Removed + stop_codon
                    Plus_2_Frameshifted = Plus_2_Frameshifted.replace("*","")
        
                    Reverse_Complement = start_codon + Reverse_Complement_Stops_Removed + stop_codon
                    Reverse_Complement = Reverse_Complement.replace("*","")
        
                    Plus_1_Frameshifted_Reverse_Complement = start_codon + Plus_1_Frameshifted_Reverse_Complement_Stops_Removed +stop_codon
                    Plus_1_Frameshifted_Reverse_Complement = Plus_1_Frameshifted_Reverse_Complement.replace("*","")
            
                    Plus_2_Frameshifted_Reverse_Complement = start_codon + Plus_2_Frameshifted_Reverse_Complement_Stops_Removed + stop_codon
                    Plus_2_Frameshifted_Reverse_Complement = Plus_2_Frameshifted_Reverse_Complement.replace("*","")



                    # I've eliminated the code to convert the nucleotide
                    # sequences into amino acid sequences
                    # because I decided to generate the amino
                    # acid sequences in a different script in a different location
                    #
                    # The script which generates the output of the amino acid
                    # translations can be found in the folder
                    # Desktop/Research/Gene_data_and_programs/Data_For_Upload_to_MySQL
                    # The file is codon_to_amino_acid_converter.py

                    # The following code will calculate the length-differences between the
                    # original nucleotide sequence and the frameshifted sequences

                    Plus_1_Frameshifted_Protein = Seq(Plus_1_Frameshifted).translate()
                    Plus_1_Frameshifted_Protein_NoCys = str(Plus_1_Frameshifted_Protein).replace("C","")

                    Plus_2_Frameshifted_Protein = Seq(Plus_2_Frameshifted).translate()
                    Plus_2_Frameshifted_Protein_NoCys = str(Plus_2_Frameshifted_Protein).replace("C","")

                    Reverse_Complement_Protein = Seq(Reverse_Complement).translate()
                    Reverse_Complement_Protein_NoCys = str(Reverse_Complement_Protein).replace("C","")

                    Plus_1_Frameshifted_Reverse_Complement_Protein = Seq(Plus_1_Frameshifted_Reverse_Complement).translate()
                    Plus_1_Frameshifted_Reverse_Complement_Protein_NoCys = str(Plus_1_Frameshifted_Reverse_Complement_Protein).replace("C","")

                    Plus_2_Frameshifted_Reverse_Complement_Protein = Seq(Plus_2_Frameshifted_Reverse_Complement).translate()
                    Plus_2_Frameshifted_Reverse_Complement_Protein_NoCys = str(Plus_2_Frameshifted_Reverse_Complement_Protein).replace("C","")

                    
                    # The output of this file will be given in the following order. All entries will be tab-delimited
                    # But for clarity, the entries have been broken up into groups
                    # Entry[1] = EnsemblUID, Entry[2] = Transcript UID
                    #
                    # Entry[3] = +1 Frameshifted Nucleotide Sequence, Entry[4] = +1 Frameshifted Protein Sequence
                    # Entry[5] = +1 Frameshifted NoCys Protein Sequence, Entry[6] = Control Designation
                    #
                    # Entry[7] = +2 Frameshifted Nucleotide Sequence, Entry[8] = +2 Frameshifted Protein Sequence
                    # Entry[9] = +2 Frameshifted NoCys Protein Sequence, Entry[10] = Control Designation
                    #
                    # Entry[11] = Reverse Complement Nucleotide Sequence, Entry[12] = Reverse Complement Protein Sequence
                    # Entry[13] = Reverse Complement NoCys Protein Sequence, Entry[14] = Control Designation
                    #
                    # Entry[15] = +1 Reverse Complement Nucleotide Sequence, Entry[16] = +1 Reverse Complement Protein Sequence
                    # Entry[17] = +1 Reverse Complement NoCys Protein Sequence, Entry[18] = Control Designation
                    #
                    # Entry[19] = +2 Reverse Complement Nucleotide Sequence, Entry[20] = +2 Reverse Complement Protein Sequence
                    # Entry[21] = +2 Reverse Complement NoCys Protein Sequence, Entry[22] = Control Designation
                    
                    print('%s\t%s\t%s\t%s\t%s\tShiftedReadingFrameControlPlusOne\t%s\t%s\t%s\tShiftedReadingFrameControlPlusTwo\t%s\t%s\t%s\tShiftedReadingFrameControlRevComp\t%s\t%s\t%s\tShiftedReadingFrameControlRevCompPlusOne\t%s\t%s\t%s\tShiftedReadingFrameControlRevCompPlusTwo' %  (EnsemblUID,TranscriptUID,Plus_1_Frameshifted,Plus_1_Frameshifted_Protein, Plus_1_Frameshifted_Protein_NoCys, Plus_2_Frameshifted,Plus_2_Frameshifted_Protein,Plus_2_Frameshifted_Protein_NoCys, Reverse_Complement,Reverse_Complement_Protein,Reverse_Complement_Protein_NoCys, Plus_1_Frameshifted_Reverse_Complement, Plus_1_Frameshifted_Reverse_Complement_Protein, Plus_1_Frameshifted_Reverse_Complement_Protein_NoCys, Plus_2_Frameshifted_Reverse_Complement,Plus_2_Frameshifted_Reverse_Complement_Protein, Plus_2_Frameshifted_Reverse_Complement_Protein_NoCys))


   
