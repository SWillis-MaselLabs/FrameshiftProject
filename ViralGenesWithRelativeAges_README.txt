The file ViralGenesWithRelativeAges contains all the data for the 34 overlapping viral genes for which the relative ages could be determined (either by phylostratigraphy or codon adaptation analysis). Each entry contains the following information:



SequenceUID - The unique ID ascribed to each entry in the original database



OrderOrGroup - The taxonomic order or group of the viral species in question



Family - The taxonomic family of the viral species in question 



Genus - The taxonomic genus of the viral species in question



Species	- The viral species that the overlapping gene came from



PValue	- If the codon degeneracy data of the viral species were available, then the gene’s relative adaptedness values were analyzed and compared with the overlapping partner’s relative adaptedness values. A Mann-Whitney U-Test was then run comparing the two distributions 



GenomeAccessionNumber - The accession number for the viral species’ genome in the NCBI database



ProteinAccessionNumber	- The accession number for the protein associated with this gene in the NCBI database



GeneName - The name of this entry’s gene 



OverlappingGene	 - The name of the gene that overlaps with this entry



GeneDesignation - This is either ancestral or novel. Ancestral corresponds to the overprinted protein, while novel is the gene which is born either wholly or in-part de novo



Verification - Corresponds to the relative ages of the genes in question. The value is either Tentative or Verified. Verified corresponds to genes which have been classified as either ancestral or novel by phylostratigraphy analysis in the literature in which they were found. Tentative corresponds to genes which were analyzed by this work using codon degeneracy. They were classified as either ancestral or novel by analyzing their codon adaptation index (CAI), and later a p-value cutoff was chosen (of .02) and all tentatively-classified genes with p-values above the cutoff were discarded. 



OverlapType - This entry is either “Terminal” or “Internal”. Terminal signifies that the overlap in question is an end-to-end overlap, in which neither gene is completely contained within the other. “Internal” signifies that one of the genes involved in the overlap is completely contained within its overlapping partner. 



LengthOfOverlap	 - The length of the overlap is the number of nucleotides in length of the overlapping sequence. 



OverlappingSequence - The overlapping sequence is the nucleotide sequence which is completely contained within the overlapping gene and its partner. It can either be the whole gene, or only a section depending on the type of overlap.

	

SequenceUsedForCodonAnalysis - The sequence used for codon analysis is related to the overlapping sequence. The sequence used for codon analysis is the overlapping sequence with zero, one or two nucleotides removed from the beginning of the sequence in order to put the overlapping sequence into the reading frame of the gene in question. A sufficient number of nucleotides were always removed from the tail of the sequence to ensure that the sequence length remained a multiple of three. 



CodingGene - The coding gene is the nucleotide sequence of the gene in the entry collected from the NCBI database

	

ProteinSequence	 - The entry’s protein sequence 



NoCysProteinSequence - The entry’s protein sequence with all cysteines removed. This is the sequence which was analyzed by IUPred



FrameshiftFromAncestral - This is either +0, +1, or +2, and is the frameshift of the gene in question with respect to its overlapping ancestral counterpart. If the entry is classified as “Ancestral”, then its frameshift will be +0. If the entry is classified as “Novel”, then its frameshift will be +1 or +2. 

	

FrameshiftFromNovel - 	This is either +0, +1, or +2, and is the frameshift of the gene in question with respect to its overlapping novel counterpart. If the entry is classified as “Novel”, then its frameshift will be +0. If the entry is classified as “Ancestral”, then its frameshift will be +1 or +2. 



OverlappingGeneAccessionNumber	- The accession number of the gene that overlaps with the entry



NoCysIUPredMeanISD - The mean ISD of the gene as predicted by the computer program IUPred



NoCysOrigDataIUPredAminoAcidISDScores	- The raw data from the IUPred run, specifically the individual values scored by each amino acid in the NyCysProteinSequence 



NoCysProteinSeqLength - the number of amino acids in the NoCysProteinSequence entry. 

