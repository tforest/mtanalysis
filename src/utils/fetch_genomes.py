#!/usr/bin/env python3

""" Fetch genomes from NCBI Genome DB
Forest Thomas
Licence GPLv3
"""
from bs4 import BeautifulSoup
import sys
import re
import requests

def read_data(filename):
    """ Returns a BS object from xml file
    """

    if filename.startswith("http"):
        contents = requests.get(filename).content
    else:
        infile = open(filename, "r")
        contents = infile.read()
    soup = BeautifulSoup(contents, "xml")
    #print(soup)
    return soup

def fetch(soup):
    """
    """
    print(soup)
    parsed = soup.find("div", {"class":"sequence"}).text
        
    parsed=re.sub('if(.*)(;)', '', parsed)
    parsed = parsed.split("\n",2)[2];
    return parsed

def write_file(output, content):
    with open (output, "w") as f:
        f.write(content)
            
def main():
    """ 
    """
    # parse data
    soup = read_data(filename=sys.argv[1])
    content = fetch(soup)
    write_file(output=sys.argv[2], content=content)

if __name__ == '__main__':
    main()
