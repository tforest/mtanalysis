#!/usr/bin/env python3

"""
args : input output
"""

from bs4 import BeautifulSoup

import sys
import re

input_file = sys.argv[1]
out_file = sys.argv[2]
files_list = []
base_url = "https://genome.jgi.doe.gov/portal"
# with open (input_file, "r") as f:
#     for line in f:
#         res = re.findall(r'url=(.*)(gz)', line)
#         try:
#             splice = res[0][0]
#             kept = splice[2:]
#             files_list.append(kept)
#         except:
#             continue

infile = open(input_file, "r")
contents = infile.read()

soup = BeautifulSoup(contents, 'xml')



for line in soup.find_all('file'):
    filetype = line.get("fileType")
    if filetype == "Assembly":
        try:
            url = line.get("url")
            res = re.findall(r'url=(.*)(gz)', url)
            files_list.append(base_url+res[0][0]+res[0][1])
        except:
            continue

with open (out_file, "w") as f:
    for filename in files_list:
        f.write(filename+'\n')
