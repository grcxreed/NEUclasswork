#!/usr/bin/env python
#BioPython_seqio.py
import sys
from Bio import SeqIO

def reverse_complement_sequence(infast):
    """this will return a SeqRecord object containing a Seq that is the reverse 
    complement to the input fasta file"""
    records = SeqIO.parse(infast,"fasta")
    reverse_records = list()

    for seqrecord in records:
        #get the reverse complement and create a new seq object
        new_seq = seqrecord.seq.reverse_complement()
        seqrecord.seq = new_seq
        reverse_records.append(seqrecord)
    return reverse_records

def main():
    SeqIO.write(revers_complement_sequence(sys.argv[1]),sys.argv[2], "fasta")

if __name__=="__main__":
    main()
