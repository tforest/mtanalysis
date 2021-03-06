#!/usr/bin/env python3

import sys
import argparse
from os import listdir
from os.path import isfile, join, isdir
import traceback
import numpy as np

import matplotlib.pyplot as plt

# custom modules
import rpsb_parser
import tools

def gen_files_list(input_source, source_type):
    lst = []
    if source_type == "FOLDER":
        lst += [input_source+"/"+f for f in listdir(input_source) if isfile(join(input_source, f))]
    elif source_type == "FILES_LIST":
        lst += input_source
    #store unique filenames in list to be returned
    lst = list(set(lst))
    lst.sort()
    return lst
def ambiguous_open(input_source):
    lst = []
    if type(input_source) == list:
        for elem in input_source:
            if isfile(elem):
                lst += gen_files_list(input_source, "FILES_LIST")
            elif isdir(elem):
                lst += gen_files_list(elem, "FOLDER")
            else:
                return False
    elif isdir(input_source):
        lst += gen_files_list(input_source, "FOLDER")
    elif isfile(input_source):
        lst += gen_files_list(input_source, "FILES_LIST")
    else:
        return False
    return list(np.unique(lst))
    
def input_source_selector(args):
    if args.f:
        lst = gen_files_list(args.f, "FILES_LIST")
    if args.r:
        lst = gen_files_list(args.r, "FOLDER")
    return lst
def print_entry(entry):
    for line in entry:
        print(line)
def load_data(args):
    """
    """
    data = {}
    files_lst = input_source_selector(args)
    print(files_lst)
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

def synteny_print(entry):
    gene = []
    for line in entry:
        if line.startswith("     gene"):
            if "complement" in line:
                sign = "-"
            else:
                sign = "+"
        if "/gene" in line:
            if not "trn" in line and not "orf" in line:
                word = line.strip().split('"')[1]
                if word not in gene:
                    gene.append(sign+word)
    return gene
def analyse_file(f, fname):
    entries = split_entries(f)
    prop = []
    synteny = {}
    # print header
    for entry in entries.values():
        orient = check_orientation(entry, eps=0.2)
        prop.append(orient['prop'])       
        # if(orient['q']=='same'):
        # print("Orientation(f/r):{}/{}({}): {} < {}".format(orient['f'],
        #                                               orient['r'],
        #                                               orient['prop'],
        #                                               orient['q'],
        #                                               fname))
        if "_Mito" in fname:
            spec_name = fname.split("_Mito")[0]
        else:
            spec_name = fname.split(".")[0]
        print(spec_name+","+str(orient['f'])+","+
              str(orient['r'])+","+
              str(orient['prop']))
        synteny[spec_name] = ",".join(synteny_print(entry))
    synt = spec_name, ",".join(synteny.values())
    return synt
    #print("Median:", np.median(prop))
    
def most_frequent(List): 
    return max(set(List), key = List.count) 
  
def align_synteny(synteny):
    synt = []
    ordered_synt = []
    for elem in synteny:
        for elem in elem[1].split(','):
            synt.append(elem)
    origin = most_frequent(synt)
    for elem in synteny:
        spec_name = elem[0]
        synt = elem[1].split(",")
        if origin in synt:
            ind = synt.index(origin)
            ordered_synt.append(','.join([spec_name]+synt[ind:]+synt[:ind]))
    return ordered_synt
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
                        newfile.append("  ORGANISM  "+newname.strip()+"\n"+
                                       "            "+
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
    parser.add_argument("-a",
                        help="Perform analysis of genomes",
                        action="store_true")
    parser.add_argument("-f",
                        help="input file(s)",
                        action="store",
                        nargs="*")
    parser.add_argument("-r",
                        help="repertory with input files",
                        action="store")
    parser.add_argument("-ids",
                        help="ids file for annotation (usually in CSV format)",
                        action="store")
    parser.add_argument("-rpsb",
                        help="parse RPS-BLAST results got from rpsbproc",
                        action="store")
    parser.add_argument("-comp",
                        help="Compare FASTA files to Domains identified by rpsb",
                        action="store")
    parser.add_argument("-tokenize",
                        help="Tokenize fasta headers",
                        action="store", nargs="*")
    parser.add_argument("-postannot",
                        help="Post annotation treatment (needs a list of fasta" \
                        "files to be parsed)",
                        action="store", nargs="*")
    parser.add_argument("-headermode",
                        help="Specify header mode for postannot" \
                        "(0, 1 or 2)",
                        action="store", nargs="*")
    parser.add_argument("-clusters",
                        help="Clusters of genes from fasta files",
                        action="store", nargs="*")
    parser.add_argument("-intersect",
                        help="Get intersection of orthogroups given a set of " \
                        "orthologous sequences groups as fasta files",
                        action="store", nargs="*")
    parser.add_argument("-selectgroup",
                        help="Select specific species from taxonomic characteristics",
                        action="store", nargs="*")
    parser.add_argument("-concat",
                        help="Concatenate all fasta sequences that have all genes in common in a group of files",
                        action="store", nargs="*")
    parser.add_argument("-renamefromid",
                        help="Rename fasta header accordingly to ID spec in taxonomy file",
                        action="store", nargs="*")
    # parser.add_argument("-ids",
    #                     help="Accession IDs file",
    #                     action="store")
    try:
        args = parser.parse_args()
        if args.f or args.r:
            # if it is a basic analysis
            data, f_lst = load_data(args)
            if args.a:
                # if analysis is chosen
                synteny = []
                print("spcies_id,nb_fwd,nb_rev,prop")
                for f in data.values():
                    fname = f_lst[list(data.values()).index(f)].split('/')[-1]
                    synteny.append(analyse_file(f, fname))
                synteny = align_synteny(synteny)
                print("---- SYNTENY ----")
                for synt in synteny:
                    print(synt)
            if args.ids:
                # add annotation mode
                add_annotation(args, f_lst)
        if args.rpsb:      
            domains = rpsb_parser.parse(args.rpsb)
            if (args.f or args.r) and args.comp:
                # if input masta files then, compare to parsing
                rpsb_parser.compare(data, f_lst, domains, out_folder = args.comp)
        if args.tokenize:
            tools.tokenize(args.tokenize)
        if args.postannot:
            fasta = ambiguous_open(args.postannot)
            #tools.post_annot(fasta)
            if args.clusters and args.ids:
                tools.remove_spurious_clusters(fasta, ambiguous_open(args.clusters),
                                               ids=args.ids)
            elif args.clusters:
                #tools.remove_spurious_clusters(fasta, ambiguous_open(args.clusters))
                tools.split_clusters(fasta, ambiguous_open(args.clusters))
            elif args.ids:
                if args.headermode:
                    tools.rename_mfannot_proteome(args.postannot, args.ids, headermode=args.headermode)
                else:
                    tools.rename_mfannot_proteome(args.postannot, args.ids)    
            else:
                tools.detect_bad_genomes(fasta)
        if args.intersect:
            tools.ortho_intersect(args.intersect)
        if args.selectgroup:
            tools.subset_seq_from_tax(args.selectgroup[0], args.selectgroup[1],
                                      args.selectgroup[2:])
        if args.concat:
            tools.list_common_genes(args.concat)
        if args.renamefromid and args.ids:
            tools.rename_from_id(args.ids, args.renamefromid)
    except:
        # if mandatory arguments are not specified
        traceback.print_exc()
        parser.print_help()
        sys.exit(0)
if __name__ == "__main__":
    main()
