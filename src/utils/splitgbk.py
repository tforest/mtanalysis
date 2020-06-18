#python3 splitgbk.py < input.gbk

from Bio import SeqIO
import sys

for rec in SeqIO.parse(sys.stdin, "genbank"):
   SeqIO.write([rec], open(rec.id + ".gbk", "w"), "genbank")
