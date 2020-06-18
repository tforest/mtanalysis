#!/usr/bin/env python3
"""
Gene by Gene : GenBank to FASTA Amino Acids (*.gbk to *.faa)
"""

from Bio import GenBank
from Bio import SeqIO
import sys

gbk_filename = sys.argv[1]
faa_filename = gbk_filename+"_converted.faa"

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for k, seq_feature in enumerate(seq_record.features) :
        if seq_feature.type=="CDS" :
            if "gene" not in seq_feature.qualifiers.keys():
                key = "product"
            else:
                key = "gene"
            assert len(seq_feature.qualifiers['translation'])==1
            #output_handle.write(">"+seq_feature.qualifiers['gene'][0]+"_"+
            output_handle.write(">"+"_".join(seq_feature.qualifiers[key][0].split(" "))+"_"+
                                gbk_filename.split('/')[-1]+
                                "_feature_"+str(k)+
                                "_from_"+seq_record.name+", ORGANISM "+
                                seq_record.annotations['organism']+","+
                                "\n"+seq_feature.qualifiers['translation'][0].strip('-')+"\n")

output_handle.close()
input_handle.close()
print("Done")
