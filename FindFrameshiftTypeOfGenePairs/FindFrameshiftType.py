# The purpose of this program is to look at a pair of same-sense overlapping genes 
# and to give the respective reading frames of the genes with respect to its overlapping
# partner.
# The program then prints the unique ID and unique pair ID associated with the overlapping
# pair along with the frameshift type in tab-delimited format.


# First we import the csv module so we can read in our file
import csv


# The filename of the desired file is declared
# This file should be in tab-delimited format
filename = 'OverlappingVirusGenes_93PairsWithSequenceData.txt'

# So that the program doesn't double-count overlapping pairs, a list of unique IDs is kept
# and used to skip genes/gene pairs which have already been analyzed
UsedSequenceUIDs = []

# The csv module is used to read in the file declared earlier. The file should
# be in tab-delimited format
with open('%s'%filename, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:

        # The file that's read in should have the following entries:
        # 1 - Unique ID associated with that particular gene
        # 2 - Unique pair ID, used to match overlapping genes. Each pair has a unique ID associated with it
        # 3 - The name of the viral species where the genes came from
        # 4 - The name of the gene (one of the genes in the overlapping pair)
        # 5 - The shared nucleotide sequence associated with both genes involved in the overlap
        # 6 - The nucleotide sequence of the gene (associated with the gene named in entry 4)
        Sequence1UID = row[0]
        PairUID_1 = row[1]
        Species_1 = row[2]
        GeneName_1 = row[3]
        OverlappingSequence_1 = row[4]
        CodingGene1 = row[5]

        # First the program makes sure that the UID hasn't been used yet. If it has, then it
        # skips it so no double-counting occurs
        if Sequence1UID in UsedSequenceUIDs:
            pass
        # If it hasn't been used yet, then the UID is added to the list of used UIDs and the file is
        # opened for a second time so that it can be scanned for its overlapping pair
        else:
            UsedSequenceUIDs.append(Sequence1UID)
            with open('%s' %filename, 'r') as f2:
                reader2 = csv.reader(f2, delimiter = '\t')
                for item in reader2:
                    # Some of the entries will be the same as initially declared, but they are given
                    # different variable names so that they can be compared with one another
                    Sequence2UID = item[0]
                    PairUID_2 = item[1]
                    Species_2 = item[2]
                    GeneName_2 = item[3]
                    OverlappingSequence_2 = item[4]
                    CodingGene2 = item[5]
                    
                    # If the pair IDs match up, then we may proceed. However, so that we don't match a gene
                    # to itself, we need to check that the two sequence UIDs don't match. If they do match, then that entry
                    # is skipped. If they don't match, then they are overlapping pairs and the second UID is added
                    # to the list of used UIDs
                    # Once we have the overlapping pair, a single OverlappingSequence variable is declared for clarity
                    # The two entries should have been the same to start with)
                    if PairUID_2 == PairUID_1:
                        if Sequence2UID not in UsedSequenceUIDs:
                            UsedSequenceUIDs.append(Sequence2UID)
                            OverlappingSequence = OverlappingSequence_2


                            # First the program looks to see if the overlap is an internal overlap
                            # i.e. the case when one gene is entirely contained within another gene
                            # If this is the case, then we find the starting index modulo 3 for where
                            # that internal gene begins in the external gene
                            # This gives us the reading frame the internal gene is found in
                            # with respect to the external gene. Since overlaps are symmetrical,
                            # it's overlapping pair's respective reading frame is 3 minus the
                            # starting index mod 3.

                            if CodingGene1 == OverlappingSequence:

                                Sequence_StartIndex = CodingGene2.index(CodingGene1)
                                FrameshiftType = SequenceStartIndex%3
                                print('%s\t%s\t+%s' %(Sequence1UID,PairUID_1,3-FrameshiftType))
                                print('%s\t%s\t+%s' %(Sequence2UID,PairUID_1,FrameshiftType))

                            # if the first gene is not entirely contained within the second, then
                            # the second is checked to see whether it is completely contained in the first
                            # if it is, a mirror image of what was done above is done to find the respective
                            # reading frames of the two genes
                            else:
                                if CodingGene2 == OverlappingSequence:
                                    Sequence_StartIndex = CodingGene1.index(CodingGene2)
                                    FrameshiftType = Sequence_StartIndex%3
                                    print('%s\t%s\t+%s' %(Sequence1UID,PairUID_1,3-FrameshiftType))
                                    print('%s\t%s\t+%s' %(Sequence2UID,PairUID_1,FrameshiftType))
                                    
                                # If the overlap is not an internal overlap, then it is a terminal overlap
                                # meaning that the ends of the two genes overlap with one another without
                                # either being completely contained within the other.
                                # The overlapping sequence is then used to find the frameshift type in
                                # the same ways as above. 
                                else:
                                    StartIndex = CodingGene1.index(OverlappingSequence)
                                    # If the starting index is 0, that means that the front of the gene being
                                    # looked at is within the overlapping region, and so the other gene is
                                    # looked to find the reading frame (this can be done because of the
                                    # symmetry of reading frames. If a gene is in the +1/+2 frame of its partner,
                                    # Then its partner is in the +2/+1 frame of that gene)
                                    if StartIndex == 0:
                                        New_StartIndex = CodingGene2.index(OverlappingSequence)
                                        Frameshift_1_from_2 = New_StartIndex%3
                                        print('%s\t%s\t+%s' %(Sequence2UID,PairUID_2,3-Frameshift_1_from_2))
                                        print('%s\t%s\t+%s' %(Sequence1UID,PairUID_1,Frameshift_1_from_2))
                                        
                                    else:
                                        Frameshift_2_from_1 = StartIndex%3
                                        print('%s\t%s\t+%s' %(Sequence1UID,PairUID_1,3-Frameshift_2_from_1))
                                        print('%s\t%s\t+%s' %(Sequence2UID,PairUID_2,Frameshift_2_from_1))

                                        # In all cases, the final output is of the form
                                        # [SEQUENCE UID] [TAB] [PAIR UID] [TAB] [FRAMESHIFT FROM OVERLAPPING PAIR]
                                        

    
