# The purpose of this program is to write two functions which take in a nucleotide
# sequence and determine whether its length is a multiple of three. If it is not,
# frontTrimmer will take off enough nucleotides from the start of the sequence so
# that the sequence's length is a multiple of 3.

# backTrimmer will also determine whether or not the length of the sequence
# is a multiple of 3, but it will remove the excess nucleotides from the end of the
# sequence instead of the start

def frontTrimmer(nucleotide_sequence):
    if len(nucleotide_sequence)%3 == 0:
        return nucleotide_sequence

    else:
        return nucleotide_sequence[len(nucleotide_sequence)%3 ::]

def backTrimmer(nucleotide_sequence):
    if len(nucleotide_sequence)%3 == 0:
        return nucleotide_sequence

    else:
        return nucleotide_sequence[: len(nucleotide_sequence)-len(nucleotide_sequence)%3 :]

