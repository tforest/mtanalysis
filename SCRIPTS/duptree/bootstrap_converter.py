#!/usr/bin/env python
# Copyright 2016, Felix Thalen

"""
Reads a Newick tree and converts the bootstrap value for each and every node in
that tree into the desired format.
"""

import sys
import argparse
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

if __name__ == '__main__':
    main()
