#!/usr/bin/env python3

""" Parse JGI browser
Forest Thomas
MIT Licence
"""

from bs4 import BeautifulSoup
import sys
import re

def read_data(filename):
    """ Returns a BS object from xml file
    """

    infile = open(filename, "r")
    contents = infile.read()
    soup = BeautifulSoup(contents, 'xml')
    return soup

def fetch_urls(soup, keywords, types_list):
    """
    """
    base_url = "https://genome.jgi.doe.gov/portal"
    files_list = []
    for line in soup.find_all('file'):
        filetype = line.get("fileType")
        if filetype in types_list:       
            try:
                url = line.get("url")
                res = re.findall(r'url=(.*)(gz)', url)
                kept_url = base_url+res[0][0]+res[0][1]
                if len(keywords) > 0:
                    for key in keywords:
                        if key in kept_url:
                            files_list.append(kept_url)
                else:
                    files_list.append(kept_url)       
            except:
                continue
    return files_list

def write_list(output, files_list):
    with open (output, "w") as f:
        for filename in files_list:
            f.write(filename+'\n')

def build_urls_list(soup, out_file, keywords=[], types_list=['Assembly']):
    """ Build a file with all URLs of selected organisms
    """
    # get urls list
    urls = fetch_urls(soup, keywords, types_list)
    # write urls
    write_list(out_file, urls)

def stats(soup, keyword=""):
    """ Compute stats on file DB
    """
    filetypes = {}
    for line in soup.find_all('file'):
        filename = line.get("filename")
        if keyword.lower() in filename.lower():
            filetype = line.get("fileType")
            if filetype not in filetypes:
                filetypes[filetype] = [line]
            else:
                filetypes[filetype].append(line)
    # print types stats        
    for filetype in filetypes:
        print(filetype, len(filetypes[filetype]))

def load_ids(filename):
    ids = []
    with open(filename, 'r') as idsfile:
        for line in idsfile.readlines():
            ids.append(line.strip())
    print(ids)
    return ids
def main():
    """ 
    """
    # parse data
    soup = read_data(filename=sys.argv[1])
    # stats(soup, keyword="mito")
    # stats(soup, keyword="")
    try:
        ids = load_ids(sys.argv[3])
        build_urls_list(soup, sys.argv[2], keywords=ids)
    except:
        print("without ids")
        build_urls_list(soup, sys.argv[2])

if __name__ == '__main__':
    main()
