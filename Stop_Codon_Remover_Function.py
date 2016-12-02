#The purpose of this file is to create a subroutine which will read
#in a nucleotide sequence, will identify whether there
#are stop codons within the body of the sequence, and if so will remove them.

def stopCodonRemover(nucleotide_sequence):
    for n in range(len(nucleotide_sequence)-1):

        # If any stop codons are found, they are replaced with the string '***'
        # This is to keep the sequence the same length so that the entire sequence
        # can be searched and the stop codons can be identified
        # before any removal takes place
        
        if nucleotide_sequence[n] == 'T' and n%3 ==0:
            if nucleotide_sequence[n+1] =='A' and nucleotide_sequence[n+2] == 'A':
                nucleotide_sequence = nucleotide_sequence[:n]+str('***') + nucleotide_sequence[n+3 : ]
            else:
                if nucleotide_sequence[n+1] == 'A' and nucleotide_sequence[n+2] == 'G':
                    nucleotide_sequence = nucleotide_sequence[:n]+str('***') + nucleotide_sequence[n+3 : ]
                else:
                    
                    if nucleotide_sequence[n+1] == 'G' and nucleotide_sequence[n+2] == 'A':
                        nucleotide_sequence = nucleotide_sequence[:n]+str('***') + nucleotide_sequence[n+3 : ]
                    else:
                        nucleotide_sequence = nucleotide_sequence
        else:
            nucleotide_sequence = nucleotide_sequence

        #The last step is to take the sequence that was generated above and to replace all of the
        # '*' characters with nothing so they are removed.
    nucleotide_sequence_with_stop_codons_removed = nucleotide_sequence.replace("*","")
    return nucleotide_sequence_with_stop_codons_removed



