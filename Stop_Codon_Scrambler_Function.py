import random

# The purpose of this program is to look for stop codons
# within the body of a nucleotide sequence. If there are, then we select
# a nucleotide corresponding to the stop codon at random and distribute it,
# at random, somewhere else in the body of the sequence. We do this until
# there are no more stop codons in the body of our sequence.




def stopCodonScrambler(scrambled_and_trimmed_nucleotide_sequence):
    true = 1
    false = 0
    stop_codons_in_body = true


    while stop_codons_in_body == true:
        
        # count is used to determine whether or not a stop codon has been found
        # If nothing is added to the count
        # then no stop codon was found and at the end of this loop
        # stop_codons_in_body is set to false.
        count = 0

        # 'n' will be used to look at the idices of our sequence.
        for n in range(len(scrambled_and_trimmed_nucleotide_sequence)-1):

            
            # For there to be a start codon, the codon must start with T, which means
            # there must be a T at an index which is a multiple of three. If this is true
            # we move on to further investigations of the codon. If not, then we continue
            # scanning the sequence
            if scrambled_and_trimmed_nucleotide_sequence[n] == 'T' and n%3 ==0:
                
                # The next if statements look to see what follows the T at the start of a codon. If any
                # of the three if statements below are true, a stop codon exists. If not, there is no
                # stop codon.
                if scrambled_and_trimmed_nucleotide_sequence[n+1] == 'A' and scrambled_and_trimmed_nucleotide_sequence[n+2] == 'A':
                    # If a stop codon is found, a random nucleotide within the codon is selected by its index 
                    random_integer_to_select_nucleotide = random.randint(n,n+2)
                    
                    # The nucleotide selected is removed from the sequence by concatenating the sequences
                    # on either side of it. 
                    scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide = scrambled_and_trimmed_nucleotide_sequence[: random_integer_to_select_nucleotide :] + scrambled_and_trimmed_nucleotide_sequence[random_integer_to_select_nucleotide+1 ::]
                    # A random integer is then chosen to select an index within the shortened sequence
                    random_integer_to_splice_sequence = random.randint(0, len(scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide)-1)
                    
                    # The nucleotide sequence is spliced at the index which was selected by the
                    # random variable above, the removed nucleotide is inserted in this location
                    # and the whole sequence is concatenated.
                    scrambled_and_trimmed_nucleotide_sequence= scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide[: random_integer_to_splice_sequence :] +scrambled_and_trimmed_nucleotide_sequence[random_integer_to_select_nucleotide]+ scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide[random_integer_to_splice_sequence: :]
                    count = count + 1

                    
                else:
                    if scrambled_and_trimmed_nucleotide_sequence[n+1] == 'A' and scrambled_and_trimmed_nucleotide_sequence[n+2] == 'G':
                        random_integer_to_select_nucleotide = random.randint(n,n+2)
                        scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide = scrambled_and_trimmed_nucleotide_sequence[: random_integer_to_select_nucleotide :] + scrambled_and_trimmed_nucleotide_sequence[random_integer_to_select_nucleotide+1 ::]
                        random_integer_to_splice_sequence = random.randint(0, len(scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide)-1)
                        scrambled_and_trimmed_nucleotide_sequence= scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide[: random_integer_to_splice_sequence :] +scrambled_and_trimmed_nucleotide_sequence[random_integer_to_select_nucleotide]+ scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide[random_integer_to_splice_sequence: :]
                        count = count + 1
                
                
                    else:
                        if scrambled_and_trimmed_nucleotide_sequence[n+1] == 'G' and scrambled_and_trimmed_nucleotide_sequence[n+2] == 'A':
                            random_integer_to_select_nucleotide = random.randint(n,n+2)
                            scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide = scrambled_and_trimmed_nucleotide_sequence[: random_integer_to_select_nucleotide :] + scrambled_and_trimmed_nucleotide_sequence[random_integer_to_select_nucleotide+1 ::]
                            random_integer_to_splice_sequence = random.randint(0, len(scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide)-1)
                            scrambled_and_trimmed_nucleotide_sequence= scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide[: random_integer_to_splice_sequence :] +scrambled_and_trimmed_nucleotide_sequence[random_integer_to_select_nucleotide]+ scrambled_and_trimmed_nucleotide_sequence_missing_nucleotide[random_integer_to_splice_sequence: :]
                            count = count + 1
                        else:
                            pass

        #If the count is zero, then no stop codons were found
        if count == 0:
            stop_codons_in_body =false

    return scrambled_and_trimmed_nucleotide_sequence
