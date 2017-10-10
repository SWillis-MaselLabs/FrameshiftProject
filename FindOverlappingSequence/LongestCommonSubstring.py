# This submodule contains the function longest_common_substring
# which takes in two strings as arguments and returns the longest
# substring that is shared by both of them. It will be used to find the
# nucleotide sequence shared by two overlapping genes.

def longest_common_substring(s1,s2):
    m = [[0]*(1+len(s2)) for i in range(1 + len(s1))]
    longest, x_longest =0,0
    for x in range(1,1+len(s1)):
        for y in range(1,1+len(s2)):
            if s1[x-1]==s2[y-1]:
                m[x][y]=m[x-1][y-1]+1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest-longest: x_longest]
