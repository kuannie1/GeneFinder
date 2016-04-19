from load import load_nitrogenase_seq
from load import load_metagenome
#import pypy
from gene_finder import *


nitrogenase = str(load_nitrogenase_seq())
metagenome = load_metagenome()
one_metagenome = metagenome[0] #a tuple

def longestSubstringMetagenomeSnippet(enzyme, metagenome_snippet):
    m = [[0] * (1 + len(metagenome_snippet)) for i in range(1 + len(enzyme))] #creates a list of lists
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(enzyme)): #x starts with 1, goes up to (1 + length of 1st string)
        for y in range(1, 1 + len(metagenome_snippet)): #y starts with 1, goes up to (1 + length of 2nd string)
            if enzyme[x - 1] == metagenome_snippet[y - 1]: 
                m[x][y] = m[x - 1][y - 1] + 1 #if the upper-left diagonal characters are equal, 
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return enzyme[x_longest - longest: x_longest]


def longestSubstringAllMetagenomeSnippets(enzyme, metagenome):
	i = 0
	ORF_lengths = []
	genome_names = [] #I think that's what you call a part of a metagenome?
	for name, metagenome_snippet in metagenome:
		ORF = longestSubstringMetagenomeSnippet(enzyme, metagenome_snippet)
		ORF_lengths.append(len(ORF))
		print i
		i += 1
	largest_size = max(ORF_lengths)
	index_that_solves_everything = a.index(max(ORF_lengths))

	return metagenome[index_that_solves_everything]

longestSubstringAllMetagenomeSnippets(nitrogenase, metagenome)