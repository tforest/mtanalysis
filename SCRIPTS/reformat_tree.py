#!/usr/bin/env python3
# Copyright 2016, Felix Thalen

"""
Reads a Newick tree and converts the bootstrap value for each and every node in
that tree into the desired format.
"""

import sys
import argparse

import hashlib

try:
    from ete3 import Tree
except:
    print >>sys.stderr, 'ETE Toolkit not found.'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--replace',
            help='replace the input file with the output', action='store_true')
    parser.add_argument('-d', '--decimal', help ='integer -> decimal, \
            default: decimal -> integer', action='store_true')
    parser.add_argument('-e', '--extension',
            help='specify extension, default: .mod.tre', action='store_true')
    parser.add_argument('-s', '--serialize',
            help='serialize leaf names and stores names correspondance in a file to recover them later', action='store_true')
    parser.add_argument('--unserialize',
            help='recover leaf names from previous correspondance file', action='store_true')
    parser.add_argument('--serfile',
            help='serialization file (optional. By default same as in input with .ser suffix)', action='store')
    parser.add_argument('tree', help='file containing the tree')
    return parser.parse_args()


def main():
    args = parse_args()

    # Use the extension specified by the user if present
    if args.extension:
        ext = args.extension
    else:
        ext = '.mod.tre'

    # Load the tree
    t = Tree(args.tree)

    if args.serialize and args.unserialize:
            print("Error! You cannot serialize and unserialize at the same time.")
            exit(0)
    
    if args.serialize:
            serial = {}

    if args.unserialize:
        unserial = {}
        if args.serfile:
            serfile = args.serfile
        else:
            serfile = args.tree+".ser"

        
        with open(serfile, "r") as unserfile:
            for unser in unserfile.readlines():
                line = unser.split(",")
                unserial[line[0]] = line[1].strip()
    for i, leaf in enumerate(t):

        hashed = str(hashlib.md5(leaf.name.encode()).hexdigest())[:10]
        if args.serialize:
            serial["Spec_"+hashed] = leaf.name
      
        if len(leaf.name.split("~")) == 3:
            # remove dots from names
            leaf.name = leaf.name.replace('.', '')
            # replace leaf name without gene name (remove cox1~ etc...)
            # also replace ~ with _ underscores 
            leaf.name = "_".join(leaf.name.split("~")[-2:])
            
        if args.unserialize:
            leaf.name = unserial[leaf.name]
            
        if args.serialize:
            leaf.name = "Spec_"+hashed
    # Iterate over the nodes, convert to desired value
    for node in t.iter_search_nodes():
        node.support = round(node.support*100)
        node.dist = round(node.dist, 5)

    # If the replace flag is set, replace the input file with the output file.
    # Otherwise create a new file with the '.mod.tre' extension
    if args.replace:
        out = args.tree
    else:
        out = args.tree + ext
    t.write(format=0, outfile=out)

    if args.serialize:
        with open(out+".ser", "w") as serfile:
            for ser in serial:
                serfile.write(str(ser)+","+serial[ser]+"\n")
    
if __name__ == '__main__':
    main()


