#!/usr/bin/env python3
#translate_APOE.py
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translateAPOE():
        records = [rec.translate(id="aa_"+rec.id, description = "translated sequence") \
            for rec in SeqIO.parse("APOE_refseq_transcript.fasta","fasta")]
        SeqIO.write(records,"apoe_aa.fasta","fasta")

if __name__=="__main__":
    translateAPOE()
