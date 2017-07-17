Read me for the file: VirusGenes_Nonoverlapping_Controls.txt 

This file contains 92 non-overlapping viral genes used to generate artificially frameshifted controls. The file is in tab-delimited format with the following headings:

SequenceUID - The unique sequence ID used in the MySQL database for reference


Family	- The viral family to which the gene belongs


Species	- The viral species to which the gene belongs


GenomeAccessionNumber - The accession number that the viral species’ genome can be accessed with in the NCBI database


ProteinAccessionNumber	- The accession number that the gene’s protein product can be accessed with in the NCBI database


GeneDesignation - There are three possible designations. The first is “CodingGene”, which indicates that the listed nucleotide and protein sequence are the original sequences corresponding to that gene, pulled from the NCBI database. The second is “ShiftedReadingFrameControlPlusOne”, which indicates that the listed nucleotide and protein sequences have been artificially shifted into the +1 reading frame of the original coding gene. Similarly, the third is “ShiftedReadingFrameControlPlusTwo” which is the coding gene artificially shifted into the +2 reading frame.


GeneName - The name of the gene corresponding to the given entry



CodingGene - The nucleotide sequence of the listed gene either in its native state or artificially frameshifted, depending on GeneDesignation 



ProteinSequence	 - The translated amino acid sequence corresponding to the CodingGene


NoCysProteinSequence - The ProteinSequence with all cysteine’s removed (see materials and methods)


NoCysIUPredMeanISD - The mean ISD of the gene, calculated by IUPred using the NoCysProteinSequence



NoCysOrigDataIUPredAminoAcidISDScores - The raw data given for each amino acid sequence in the NoCysProteinSequence as given by IUPred



NoCysProteinSeqLength - The number of amino acids in NoCysProteinSequence 