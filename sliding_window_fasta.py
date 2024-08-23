#!/usr/bin/env python
# sliding_window_fasta.py

import sys
from Bio import SeqIO
import argparse

def sliding_window(k, string):
    '''Returns a list of all k-mers in a given string.'''
    return [string[start:start + k] for start in range(len(string) - k + 1)]

def gc_content(seq):
    """Returns [0,1], the fraction of GCs in the input string."""
    gc_count = seq.lower().count('g') + seq.lower().count('c')
    return gc_count / float(len(seq))

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate GC content for k-mers in a FASTA file.")
    parser.add_argument("kmer_length", type=int, help="Length of the k-mers.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    for record in SeqIO.parse(args.fasta_file, "fasta"):
        print(f">{record.description}")
        for kmer in sliding_window(args.kmer_length, str(record.seq)):
            print(f"{kmer}\t{gc_content(kmer):.2f}")

if __name__ == "__main__":
    main()

