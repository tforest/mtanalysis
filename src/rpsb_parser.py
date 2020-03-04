
""" Parse RPSB
A direct output parser for RPSBPROC (CD-SEARCH like)
ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/
"""
import mtanalysis
from sys import argv
import fasta
DEBUG = False

def parse_file(file):
    filename = file.split('/')[-1]
    # the hits list found by rpsb
    hits = []
    print("Parsing", filename, "...")
    with open(file) as f:
        for line in f.readlines():
            if line.startswith("#") or not line.startswith("QUERY"):
                # pass if it's not a hit line description
                continue
            else:
                hits.append(">"+line.split("\t")[-1].strip())
    return filename, hits
def parse(files_folder):
    domains = {}
    files = mtanalysis.gen_files_list(input_source=files_folder, source_type="FOLDER")
    for elem in files:
        key, domains[key] = parse_file(elem)
    return domains
def write_res(fasta_kept, fasta_unk, out_file):
    with open(out_file+'.kept', 'w') as fasta_kept_f:
        for line in fasta_kept:
            fasta_kept_f.write(line+'\n')
    with open(out_file+'.unk', 'w') as fasta_unk_f:
        for line in fasta_unk:
            fasta_unk_f.write(line+'\n')
    if DEBUG:
        for line in unkown:
            if line.startswith('>'):
                print(line)
def compare(fastas, files_list, domains, out_folder):
    for file in files_list:
        fasta_kept, fasta_unk = fasta.parse(file, domains, extension = '.asn.out')
        out_file = out_folder + '/' + file.split('/')[-1]
        print("Separating file", file, "...")
        write_res(fasta_kept, fasta_unk, out_file)
