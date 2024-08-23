#!/usr/bin/env python3
#BioPython_seqio.py
from Bio.Seq import Seq
from Bio import SeqIO

#reads a fasta file with SeqIO
def transcribeYeast():
    records = [rec.reverse_complement(id="rc_"+rec.id, description="reverse complement") \
            for rec in SeqIO.parse("yeast.fasta","fasta")]
   # outputs a new fasta file whose contents 
   # are the reverse compliment of the sequences 
   # from the original fasta file
    SeqIO.write(records, "reverse_complement_yeast.fasta","fasta")

if __name__=="__main__":
    transcribeYeast()

