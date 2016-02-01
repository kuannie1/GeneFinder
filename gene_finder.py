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
    else:
	    return 'C'

    # TODO: implement this, remember to ctrl+B!!!
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    direct_complement = ""
    need_way_to_reverse = ""
    for char in dna:
        t = get_complement(char)
        direct_complement = direct_complement + t
    #HOW DO I DO REVERSE???
    for i in range(0, (len(direct_complement))):
        need_way_to_reverse = need_way_to_reverse + direct_complement[-1 - i]	
    return need_way_to_reverse
    # TODO: implement this
    pass


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
    #Stop codons: TAG, TAA, and TGA
    group = []
    dna_sequence = []
    last = (-1 * (len(dna) %3)) 
    no_stop_codon = True
#turned dna string into groups of 3 within a list
    for i in range(int(len(dna)/3)):
    	nucleotide = dna[3*i] + dna[3*i + 1] + dna[3*i + 2]
    	group.append(nucleotide)
#rules out the elements from the stop-codon until the end

    for j in group:
    	if (j != 'TGA') and (j != 'TAA') and (j != 'TAG'):
    		dna_sequence.append(j)
    	else:
    	    no_stop_codon = False
            break
#combines all the happy nucleotides together for a great result:
    dna_sequence = "".join(dna_sequence)
    if ((len(dna) % 3) != 0) and no_stop_codon:
    	return dna
    
#    elif ((len(dna) % 3) != 0) and not(no_stop_codon):
 #       return dna    
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

    """
    i = 0 #initial index condition
    group = [] 
    while i < len(dna): #looping through dna
    	codon = dna[i:i+3]
    	if codon == "ATG": #checks to see if codon starts
    		ORF_string = rest_of_ORF(dna[i:])
    		group.append(ORF_string)
    		i = len(ORF_string) + i
    	else: #updates so it can start again, looks at the next codon
    		i = i + 3
    		
    return group


    # TODO: implement this
    pass


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
    """
   
    group = [] 
    group = group + find_all_ORFs_oneframe(dna)
    dna = dna[1:]
    group = group + find_all_ORFs_oneframe(dna)
    dna = dna[1:]
    group = group + find_all_ORFs_oneframe(dna)
    return group
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> 
    """
    original = dna
    complement = get_reverse_complement(dna)
    group = []
    group = group + find_all_ORFs(original)
    group = group + find_all_ORFs(complement)
    return group
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose = True)
