# The purpose of this program is to create a function which takes a nucleotide sequence in as an
# input, then which checks the sequence to determine whether or not it starts with a
# start codon and ends with a stop codon. In both cases, if it does, the function will strip the
# nucleotide sequence of these codons and will produce as an output the stripped nucleotide sequence
# the start codon and the stop codon as three separate strings
# If no stop/start codons are found, the program will return the
# string '***' in their place, which is easily removable by the user down the road






def removeStartStop(Input_nucleotide_sequence):
    start_codon = str('ATG')
    
    # The first if statement checks to see whether the input nucleotide sequence
    # starts with a start codon. If it does the sequence is spliced to remove the start
    # codon. 
    if Input_nucleotide_sequence[0]=='A' and Input_nucleotide_sequence[1]=='T' and Input_nucleotide_sequence[2]=='G':
        nucleotide_sequence_without_start = Input_nucleotide_sequence[3 : :]

    # If it does not start with a start codon, the original sequence is returned and the start_codon
    # variable is returned as an empty string.
    else:
        nucleotide_sequence_without_start = Input_nucleotide_sequence
        start_codon = "***"

    trimmed_nucleotide_sequence = nucleotide_sequence_without_start

    # The next set of conditional statements determine whether or not the sequence ends
    # with a stop codon. If it does, it returns the stop codon associated with the sequence
    # as well as the nucleotide sequence (with the start codon removed) with its stop codon removed
    # If it does not, it returns the strung '***' for the stop_codon variable and returns the
    # original sequence.

    # When this function is used in a program, if the user wishes to attach the start and stop
    # codons to the body of the sequence, but does not with to have the string
    # '***' in their sequence, then use the following command on the resulting
    # sequence: sequence_string.replace("*","")

    if len(nucleotide_sequence_without_start) == 0:
        stop_codon = str('***')
    else:
        
        if nucleotide_sequence_without_start[len(nucleotide_sequence_without_start)-3] == 'T':
    
            if nucleotide_sequence_without_start[len(nucleotide_sequence_without_start)-2] == 'A' and nucleotide_sequence_without_start[len(nucleotide_sequence_without_start)-1]=='G':
        
                trimmed_nucleotide_sequence = nucleotide_sequence_without_start[: len(nucleotide_sequence_without_start)-3 :]
                stop_codon = str('TAG')
        
            else:
                if nucleotide_sequence_without_start[len(nucleotide_sequence_without_start)-2]=='A' and nucleotide_sequence_without_start[len(nucleotide_sequence_without_start)-1] == 'A':
        
                    trimmed_nucleotide_sequence = nucleotide_sequence_without_start[: len(nucleotide_sequence_without_start)-3 :]
                    stop_codon = str('TAA')
        
                else:
                    if nucleotide_sequence_without_start[len(nucleotide_sequence_without_start)-2] == 'G' and  nucleotide_sequence_without_start[len(nucleotide_sequence_without_start)-1] =='A':
        
                        trimmed_nucleotide_sequence = nucleotide_sequence_without_start[: len(nucleotide_sequence_without_start)-3 :]
                        stop_codon = str('TGA')
                    else:
                        stop_codon = str('***')
        else:
            trimmed_nucleotide_sequence = nucleotide_sequence_without_start
            stop_codon = str('***')

    return trimmed_nucleotide_sequence, start_codon, stop_codon
