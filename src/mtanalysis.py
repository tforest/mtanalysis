#!/usr/bin/env python3

import sys
import argparse
from os import listdir
from os.path import isfile, join
import traceback
import numpy as np

import matplotlib.pyplot as plt

def gen_files_list(args):
    lst = []
    if args.f:
        lst += args.f
    if args.r:
        lst += [args.r+"/"+f for f in listdir(args.r) if isfile(join(args.r, f))]
    #store unique filenames in list to be returned
    lst = list(set(lst))
    lst.sort()
    return lst
def print_entry(entry):
    for line in entry:
        print(line)
def load_data(args):
    """
    """
    data = {}
    files_lst = gen_files_list(args)
    for k, file in enumerate(files_lst):
        print("Reading file: {}/{} ({}%)".format(1+k, len(files_lst),
                                               100*(k+1)/len(files_lst)))
        with open(file, "r") as f:
            data[len(data.keys())] = f.read().splitlines()
    return data, files_lst
def split_entries(f):
    k = 0
    entries = {}
    for line in f:
        if line.startswith("//"):
            k+=1
            # do not parse line if it is //
            continue
        if k not in entries.keys():
            entries[k] = []
        entries[k].append(line)
    return entries
def check_orientation(entry, eps=0.05):
    reverse = 0
    forward = 0
   
    for line in entry:
        if line.startswith("     gene"):
            if "complement" in line:
                reverse += 1
            else:
                forward += 1

    if reverse == 0:
        prop = 1
    else:
        prop = forward/reverse
    if prop > 0 + eps and prop < 1 - eps:
        qualif = 'diff'
    else:
        qualif = 'same'
    stats = {'f':forward,
             'r':reverse,
             'prop':prop,
             'q':qualif
    }
    return stats
def analyse_file(f, fname):
    entries = split_entries(f)
    prop = []
    for entry in entries.values():
        orient = check_orientation(entry, eps=0.2)
        prop.append(orient['prop'])       
        # if(orient['q']=='same'):
        print("Orientation(f/r):{}/{}({}): {} < {}".format(orient['f'],
                                                      orient['r'],
                                                      orient['prop'],
                                                      orient['q'],
                                                      fname))
    print("Median:", np.median(prop))
    
def add_annotation(args,f_lst):
    ids = []
    annot = []
    with open(args.ids, "r") as metadata:
        # jump first line
        metadata.readline()
        for line in metadata.readlines():
            ids.append(line.strip().split(',')[-1])
            annot.append(line.strip().split(','))

    for f in f_lst:
        with open(f, "r") as fi:
            newfile = []
            newname = False
            oflag = False
            for line in fi:
                if line.startswith("  ORGANISM"):
                    oflag = True
                    name = f.split('/')[-1].split('.')[0].split('_Mito')[0]
                    for id in ids:
                        if id.strip() in name.strip():
                            newname = annot[ids.index(id)][-3]
                            annotation = annot[ids.index(id)]
                            break
                    if newname!=False:
                        newfile.append("  ORGANISM  "+newname.strip()+"\n"+"            "+
                                       "; ".join(annotation[0:5]).strip()+"\n")
                    else:
                        newfile.append("  ORGANISM  "+name+"\n")
                    newname = False
                elif line.startswith("            ") and oflag:
                    continue
                else:    
                    newfile.append(line)
                    oflag = False
        with open(f, "w") as fout:
            for line in newfile:
                fout.write(line)
def main():
    """ 
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",
                        help="input file(s)",
                        action="store",
                        nargs="*")
    parser.add_argument("-r",
                        help="repertory with input files",
                        action="store")
    parser.add_argument("-ids",
                        help="ids file",
                        action="store")
    try:
        args = parser.parse_args()
        data, f_lst = load_data(args)
        for f in data.values():
            fname = f_lst[list(data.values()).index(f)].split('/')[-1]
            analyse_file(f, fname)
        if args.ids:
            add_annotation(args, f_lst)
    except:
        # if optionnal arguments are not specified
        traceback.print_exc()
        parser.print_help()
        sys.exit(0)
if __name__ == "__main__":
    main()
