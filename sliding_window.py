#!/usr/bin/env python
# sliding.py
'''
This script generates all k-mers of a specified length from a given string.

Usage: 
python3 sliding.py <kmer_size> <string>

Arguments:
kmer_size: The length of each k-mer.
string: The string from which k-mers are to be generated.

Example:
python3 sliding.py 3 ATTGTCT
'''

import sys

def get_kmers(k, string):
    ''' Returns a list of all k-mers in the given string. '''
    kmers = []
    end = len(string) - k + 1
    for start in range(0, end):
        kmers.append(string[start:start + k])

    return kmers

if __name__ == "__main__":

    # Check to make sure there are at least two arguments
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires 2 arguments: a kmer size and a string")
    
    k = int(sys.argv[1])
    string = sys.argv[2]

    if k <= 0:
        raise ValueError("kmer size must be a positive integer")
    if k > len(string):
        raise ValueError("kmer size cannot be larger than the length of the string")

    kmers = get_kmers(k, string)
    for i, kmer in enumerate(kmers):
        print(f"{i}\t{kmer}")

