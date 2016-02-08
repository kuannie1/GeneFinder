# -*- coding: utf-8 -*-
"""
gene_finder.py: compilation of functions

@author: Anne Ku

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if (nucleotide == 'A'):
        return 'T'
    elif (nucleotide == 'T'):
        return 'A'
    elif (nucleotide == 'C'):
        return 'G'
    elif (nucleotide == 'G'):
	    return 'C'
    else:
		return "Try again with an actual nucleotide"


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("TGAACGCGGTAGCGTAC")
    'GTACGCTACCGCGTTCA'
    """
    #Uses get_complement to return the complimentary string
    direct_complement = ""
    need_way_to_reverse = ""
    for char in dna:
        t = get_complement(char)
        direct_complement += t
    return direct_complement[::-1]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTATAAAA")
    'ATGTATAAAA'
    >>> rest_of_ORF("ATGAGATGAG")
    'ATGAGA'
    """
    stop_codons = ['TAG', 'TAA', 'TGA']
    group = []
    dna_sequence = []
    last = (-1 * (len(dna) %3)) 


    no_stop_codon = True
#turned dna string into groups of 3 within a list
    for i in range(int(len(dna)/3)):
    	nucleotide = dna[3*i] + dna[3*i + 1] + dna[3*i + 2]
    	group.append(nucleotide)


#rules out the elements from the stop-codon until the end
    for codon in group:
    	if codon not in stop_codons:
    		dna_sequence.append(codon)
    	else:
    	    no_stop_codon = False
            break
#combines all the happy nucleotides together for a great result:
    dna_sequence = "".join(dna_sequence)
    if ((len(dna) % 3) != 0) and no_stop_codon:
    	return dna
    else:
    	return dna_sequence


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATATGTGTAGATAGATGTGCCC")
    ['ATGCATATGTGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATATGTGTAGAATGATGTGCCC")
    ['ATGCATATGTGTAGAATGATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATTGAATGTGTAGATAAATG")
    ['ATGCAT', 'ATGTGTAGA', 'ATG']
    """
    i = 0 
    group = [] 
    while i < len(dna): #looping through dna
    	codon = dna[i:i+3]
    	if codon == "ATG": #checks to see if codon starts
    		ORF_string = rest_of_ORF(dna[i:])
    		group.append(ORF_string)
    		i += len(ORF_string)
    	else: #updates so it can start again, looks at the next codon
    		i += 3
    		
    return group


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCGAATGTAGCATGAATGA")
    ['ATGCGAATG', 'ATGAATGA', 'ATGA']
    >>> find_all_ORFs("ATGCGATAAATGTAGCATGAATGA")
    ['ATGCGA', 'ATG', 'ATGAATGA', 'ATGA']
    """
    group = find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])
    return group


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']

    """
    complement = get_reverse_complement(dna)
    group = find_all_ORFs(dna) + find_all_ORFs(complement)
    return group

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_ORFs = find_all_ORFs_both_strands(dna)
    ORF_lengths = []
    #Finds the lengths of all ORFs and returns their lengths in new list
    for ORF in all_ORFs:
    	one_length = len(ORF)
    	ORF_lengths.append(one_length)
    #locates longest ORF (self-explanatory)
    index_location = ORF_lengths.index(max(ORF_lengths))

    return all_ORFs[index_location]


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    #shuffling the right # of times
    incrementing_trial = 0
    output_ORF = []
    while incrementing_trial < num_trials:
    	dna = shuffle_string(dna)
    	ORF = longest_ORF(dna)
    	output_ORF.append(ORF) 

    	incrementing_trial += 1

    #like before, this is for finding the longest ORF
    ORF_lengths = []    
    for ORF in output_ORF:
    	one_length = len(ORF)
    	ORF_lengths.append(one_length)

    index_location = ORF_lengths.index(max(ORF_lengths))

    return output_ORF[index_location]


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    group = []
    protein = ""
    #Groups strand into codons
    for i in range(int(len(dna)/3)):
    	nucleotide = dna[3*i] + dna[3*i + 1] + dna[3*i + 2]
    	group.append(nucleotide)
    #Conversion
    for codon in group:
    	amino_acid = aa_table[codon]
    	protein += amino_acid

    return protein


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    threshold = longest_ORF_noncoding(dna, 1500)
    amino_acid_sequences = []
    #next, find all ORFs on both strands, and then return a list 
    #containing the amino acid sequence encoded by any ORFs that 
    #are longer than the threshold computed above using 
    #longest_ORF_noncoding

    ORFs_both_strands = find_all_ORFs_both_strands(dna)
    longer_than_threshold = []
    shorter_than_threshold = []

    for ORF in ORFs_both_strands:
		if ORF > threshold:
			longer_than_threshold.append(ORF) #The desired list!



    for ORF in longer_than_threshold:
		amino_acid_sequences.append(coding_strand_to_AA(ORF))

    return amino_acid_sequences


if __name__ == "__main__":
    import doctest
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print gene_finder(dna)
    doctest.testmod()
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose = True)