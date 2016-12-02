import random

#The purpose of this file is to create a function which takes in a nucleotide
#sequence which has had its start and stop codons removed. It then randomly
#reassigns indices to the nucleotides within the sequence and returns a sequence
#which contains the same nucleotides as the original sequence, but in random
#order



def sequenceScrambler(trimmed_nucleotide_sequence):

    # An empty string is generated to append strings of
    # nucleotides to later.
    scrambled_and_trimmed_nucleotide_sequence = ""

    # An empty list is created to keep track of integers
    # so we do not generate duplicates
    Used_integers = []

    
    # So long as our list of used integers is less than the
    # length of the nucleotide sequence, we continue generating
    # new random integers. The random integers will correspond
    # to the indices of the nucelotides within the sequence
    while len(Used_integers)<len(trimmed_nucleotide_sequence):

        # The random integer must correspond to the index of a
        # nucleotide so it must be between 0 and 1 less than the length
        # of the sequence (since counting starts at 0 for indices)
        random_integer = random.randint(0,len(trimmed_nucleotide_sequence)-1)

        # If the integer we have generated has not yet been used
        # we keep it and add it to the list of used integers.
        if random_integer not in Used_integers:
            Used_integers.append(random_integer)
            
        # The number is then used to select the nucleotide corresponding
        # to that index value and concatenates it onto the
        # nucleotide sequence.
            scrambled_and_trimmed_nucleotide_sequence = scrambled_and_trimmed_nucleotide_sequence + trimmed_nucleotide_sequence[random_integer]
            
        # If it is already in the list of used integers, we
        # do nothing and generate a new number
        else:
            pass

    return scrambled_and_trimmed_nucleotide_sequence
