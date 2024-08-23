#!/usr/bin/env python
# sliding.py
'''usage python3 sliding.py 3 ATTGTCT
'''

import re
import sys

def get_kmers(k, string):
    ''' Returns a list of all k-mers in the given string
    '''
    kmers = []
    end = len(string) - k + 1
    for start in range(0, end):
        kmers.append(string[start:start + k])

    return kmers\

if __name__ == "__main__":

    # Check to make sure there are at least two arguments
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires 2 arguments: a kmer size and then a string")
    k = int(sys.argv[1])
    string = sys.argv[2]

    kmers = get_kmers(k, string)
    for i in range(len(kmers)):
        print("{}\t{}".format(i, kmers[i]))

