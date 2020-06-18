#!/usr/bin/env python3

from Bio import AlignIO
from Bio import Alphabet
import sys

input_clustal = sys.argv[1]
AlignIO.convert(input_clustal, "clustal", input_clustal+".nexus", "nexus", alphabet=Alphabet.generic_dna)
