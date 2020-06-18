#!/usr/bin/env python3

import sys
import os

KMTC_OUTPUT = sys.argv[1]
 
with open(KMTC_OUTPUT, "r") as input:
    for line in input.readlines():
        line = line.strip()
        if line.startswith("Cluster #"):
            cluster = "CLUST_"+''.join(line.split(' ')[2:])
            os.mkdir(cluster)
        if line.startswith("Tree #"):
            tree_name = line.split("\t")[1].strip()[2:]
            tree = ''.join(line.split("\t")[2:])
            with open(cluster+"/"+tree_name+".nw", "w") as out:
                out.write(tree)
