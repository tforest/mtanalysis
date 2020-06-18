#!/usr/bin/env python3
from Bio import SeqIO
from sys import argv
inputf=argv[1]

records = SeqIO.parse(inputf, "fasta")
count = SeqIO.write(records, inputf+'.stockholm', "stockholm")
print(inputf,": Converted %i records" % count)
