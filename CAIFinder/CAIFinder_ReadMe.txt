The purpose of this file is to calculate the CAI value of nucleotide sequences.

The CAI is used to determine the degree to which the codon usage pattern of a gene or gene segment matches the overall genome of the species. 

The program takes in two tab-delimited files. 

1) A file that contains all of the relative adaptedness values calculated for each codon in the species. The file should have the following entries:
	1) UID
	2) Species
	3) Designation (this says whether the values are raw codon count or relative adaptedness. Only entries with Designation = “RelativeAdapted” will be read.
	4-68) The relative adaptedness values for codons in the following order:
	UUU, UUC, UUA, UUG, UCU, UCC, UCA, UCG, UAU, UAC, UAA, UAG, UGU, UGC, UGA
	UGG, CUU, CUC, CUA, CUG, CCU, CCC, CCA, CCG, CAU, CAC, CAA, CAG, CGU, CGC, 
	CGA, CGG, AUU, AUG, AUA, AUG, ACU, ACC, AGA, ACG, AAU, AAC, AAA, AAG, AGU, AGC, AGA, 	AGG, GUU, GUC, GUA, GUG, GCU, GCC, GCA, GCG, GAU, GAC, GAA, GAG, GGU, GGC, GGA, GGG

2) A file that contains the nucleotide sequences from the genes that need CAI values ascribed to them. This file should have the following entries:
	1) UID
	2) Species
	3) Gene Name
	4) Nucleotide Sequence

For the overlapping genes in this project, only the overlapping nucleotide sequence (in the frame of the gene in question) was used. 

The program then goes through the nucleotide sequence in question and counts the occurrences of each codon. It then uses the relative adaptedness values pulled from the first file to calculate the CAI of the gene or gene segment by taking the geometric mean of the relative adaptedness values. 