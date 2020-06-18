#!/usr/bin/env python3

import sys
import re

name = sys.argv[1].split("/")[-1].split(".")[0]
#print(name)
def parse_taxono_line(tab):
    schema = ['mycota',
              'mycotina',
              'mycetes',
              'ales',
              'aceae']
    taxo = {0:'',
            1:'',
            2:'',
            3:'',
            4:''}
    for i, elem in enumerate(schema):
       for ele in tab:
           if ele.endswith(elem):
               taxo[schema.index(elem)] = ele
               
    return [tab[2]]+list(taxo.values())

with open(sys.argv[1], "r") as inputf:
    lines = inputf.readlines()
    for i, line in enumerate(lines):
        if line.startswith("  ORGANISM"):
            organisms = lines[i:i+4]
            break
tax = []
for elem in organisms[1:]:
    string = elem.strip().split("; ")
    tax+=([i.strip('.;') for i in string])

taxo = parse_taxono_line(tax)
for taxa in taxo:
    print(taxa, end=", ")
print(organisms[0].split('ORGANISM')[1].strip()+", "+name)
#print(name.strip())
