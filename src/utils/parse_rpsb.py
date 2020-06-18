#!/usr/bin/env python3

""" Parse RPSB
A direct output parser for RPSBPROC (CD-SEARCH like)
ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/
"""

from sys import argv

with open(argv[1]) as f:
    for line in f.readlines():
        if line.startswith("#") or not line.startswith("QUERY"):
            continue
        else:
            print(line.split("\t")[-1].strip().split(",")[0])
