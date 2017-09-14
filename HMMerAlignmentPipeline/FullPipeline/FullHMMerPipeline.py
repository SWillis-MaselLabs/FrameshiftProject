# The purpose of this file is to read in two tab-delimited files
# that have been downloaded from a MySQL database. The first
# is a file containing overlapping genes, the second is a file
# containing non-overlapping controls.

# The Overlapping genes file will have the following entries:
# [1] UID
# [2] Gene Name
# [3] Overlapping Nucleotide Sequence (in the same frame as the
# overall gene)

# The controls file will have the following entries:
# [1] UID
# [2] Gene Designation (either Coding Gene, Scrambled, or
# Artificially Frameshifted)
# [3] Protein Sequence

# We want to run the overlapping sections of the overlapping genes,
# the frameshifted controls and the coding genes through the program
# HMMer to see if there are any homologous genes. This program will utilize
# phmmer to assign homology group IDs. It will then output two files:
# [1] The UIDs and homology group IDs for the overlapping genes
# in tab-delimited format
# [2] The same as the above, but for the non-overlapping controls.

# First, declare the two files you wish to use:

OverlappingGenesFilename = 'OverlappingVirusGenes_93PairsWithSequenceData.txt'
ControlsFilename = 'VirusGenes_Nonoverlapping_Controls_Complete.txt'

############################################################
############################################################
############################################################
# The necessary modules and submodules are imported
import os
import sys
import csv
import os.path
from ProteinConversion import ProteinConversion
from Fasta_File_Maker import Fasta_File_Maker
from SetConsolidator import SetConsolidator

# ProteinConversion converts the overlapping nucleotide
# sequences into protein sequences, and then combines the overlapping genes
# file and the non-overlapping controls file into one large tab-delimited
# file
ProteinConversion(OverlappingGenesFilename,ControlsFilename)
#os.remove('./%s'%OverlappingGenesFilename)
#os.remove('./%s'%ControlsFilename)

# Fasta_File_Maker converts the large tab-delimited file into Fasta format
# so that HMMer can be used. Two separate fasta files are made. They are
# identical. One is used as the query file, the other serves as the
# database. 
Fasta_File_Maker()
# phmmer is then used to compare each sequence in the query file to the
# database file. The output is saved as pHMMerOutput.txt
os.system('phmmer --tblout pHMMerOutput.txt FastaFile1.txt FastaFile2.txt')
# Once HMMer is run, the fasta files are removed from the directory
os.remove('./FastaFile1.txt')
os.remove('./FastaFile2.txt')

# We now begin the process of assigning homology group IDs
HMMerOutputFile = 'pHMMerOutput.txt'
FileForComparison = 'AllSequencesFile.txt'

# Because of the formatting of the tblout file from hmmer, we need to specify
# unwanted characters so we can get rid of them. In this case, the unwanted
# characters are spaces, of which there are many

# The list homologousPairs will be used to store genes that hit one another
# pairwise (as output from HMMer)

# TotalList will be used to keep track of the unique IDs of the genes
# that had a hit in HMMer. This way, we can later assign unique IDs
# to all the sequences that only hit themselves 

unwanted = ['']
homologousPairs = []
TotalList = []

# Using the csv module, we open the first file

with open('%s' %HMMerOutputFile, 'r') as f:
    reader = csv.reader(f, delimiter = ' ')
    for row in reader:
        # We declare an empty list so we can reformat the hmmer output file
        newRow = []
        # we look at each element line by line in the file and if it's not
        # a blank space, then we add it to the newly reformatted row
        for element in row:
            if element not in unwanted:
                newRow.append(element)
        # If a line begins with a hash, then it's either a column heading
        # notation in the hmmer output and we ignore it
        if newRow[0][0] != '#':
            # The first element in the new row is the UID of the query gene
            # The second is the gene that it hit
            # We combine these into a list containing both genes
 
            gene1 = newRow[0]
            gene2 = newRow[2]
            pair = [gene1,gene2]
            # We sort the genes so that we can compare lists (since lists
            # are order-dependent)
            pair = sorted(pair)
            # If the genes are the same, we discard them (since every gene
            # hits itself)
            if gene1 != gene2:
                # The TotalList of gene IDs that had a hit then has
                # these genes appended to it
                TotalList = TotalList + pair
                # Since if A hits B, then B hits A, if we count all
                # non-identical entries that hit, we double-count, so
                # before we add the pair to the list of homologous pairs,
                # we check to see whether they're already in there.
                if pair not in homologousPairs:
                    homologousPairs.append(pair)


# The output from SetConsolidator will be a list of disjoint
# lists containing UIDs that hit eachother in the HMMer run
homologyGroups_NoDuplicates = SetConsolidator(homologousPairs)

# We define homologyGroupID to be 0. We will sequentially add 1
# to it for each sublist in homologyGroups_NoDuplicates
homologyGroupID = 0
# We will also keep track of which UIDs have a homologyGroupID
# assigned to it. Any UIDs that do not have a homologyGroupID at
# the end of this program will be assigned a unique one
GenesWithHomologyGroupIDs = []
Genes_NoHomologs = []

# The output files are declared here. These will be where the desired
# outputs will be printed. One for the non-overlapping controls, one
# for the overlapping genes. The files will be saved one directory up
# so they're easy to locate

save_path = '../'
ControlFileFullName = os.path.join(save_path, 'ControlsHomologyGroups.txt')
OverlappingFileFullName = os.path.join(save_path, 'OverlappingHomologyGroups.txt')
ControlsFile = open('%s'%ControlFileFullName,'w')
OverlappingFile = open('%s'%OverlappingFileFullName,'w')

# Each sublist in homologyGroups_NoDuplicates counts as a homology
# group. Each sublist will have a unique homologyGroupID assigned to it
# as an integer in the last position of each sublist
for group in homologyGroups_NoDuplicates:
    homologyGroupID+=1
    group.append(homologyGroupID)

# Next, we print the UIDs and their homology groups in tab-delimited
# format to output files
for group in homologyGroups_NoDuplicates:
    for element in group:
        if group.index(element) != len(group)-1:
            # Because we had to combine the UIDs into one large file
            # to run HMMer, we added a 'c' to each of the control's
            # UIDs to distinguish them from the overlapping genes'
            # UIDs. As we're printing our output file, we use this
            # c to distinguish which gene is an overlapping gene
            # and which is a control. The c is removed for the
            # output file so the homology UIDs can be easily
            # uploaded into the MySQL database
            if element[0] != 'c':
                output = [str(element),str(group[len(group)-1])]
                OverlappingFile.write('\t'.join(output[0:]) + '\n')
            else:
                output = [str(element[1:]),str(group[len(group)-1])]
                ControlsFile.write('\t'.join(output[0:]) + '\n')
            # Every UID that has a HomologyUID is appended to this
            # list so that any that doesn't appear will have a
            # unique Homology ID given to it
            GenesWithHomologyGroupIDs.append(element)

# finally, the large file with all the sequences is opened and
# the UIDs are looked at to see which ones need a unique homology group ID. 

with open('%s'%FileForComparison,'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        UID = row[0]
        GeneName = row[1]
        Sequence = row[2]

        if UID not in GenesWithHomologyGroupIDs:
            homologyGroupID+=1
            if UID[0] != 'c':
                output = [str(UID),str(homologyGroupID)]
                OverlappingFile.write('\t'.join(output[0:]) + '\n')
            else:
                output = [ str(UID[1:]),str(homologyGroupID)]
                ControlsFile.write('\t'.join(output[0:]) + '\n')


# Finally, we close the output files.
OverlappingFile.close()
ControlsFile.close()



