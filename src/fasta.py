DEBUG = True

def parse_domains(filename, domains, extension='.asn.out'):
    """
    args:
    extension is the suffit added by the rpsb 
    """    
    keep = False
    kept = []
    unknown = []
    base_filename = filename.split("/")[-1]
    with open(filename) as fasta:
        for line in fasta.readlines():
            line = line.strip()
            if line.startswith('>'):
                keep = False
                if line in domains[base_filename+extension]:
                    if DEBUG: print(line,
                                    domains[base_filename+extension].index(line))
                    keep = True
            if keep: kept.append(line)
            if line not in kept: unknown.append(line)
    return kept, unknown

def parse(filename):
    """
    Parse a regular fasta file and split in list of entries inside a dict 
    """
    entries = {}
    base_filename = filename.split("/")[-1]
    with open(filename) as fasta:
        for line in fasta.readlines():
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
            else:
                if key in entries.keys():
                    entries[key] += line
                else:
                    entries[key] = line

    return entries

def set_header(stream, mode, names=None):
    new_stream = {}
    spec_name = ""
    spec_id = ""
    for head, seq in stream.items():
        for id, name in names.items():
            if id in head:
                spec_id = id
                spec_name = name
                if mode == 0:
                    new_head = head.split(";")[0]+"~"+'_'.join(name.split(" "))+"~"+id
                elif mode == 1:
                    new_head = head.split("_")[0]+"~"+'_'.join(name.split(" "))+"~"+id
                else:
                    print("Error ! Please select a correct mode for renaming headers!")
                    exit(0)
                new_stream[new_head] = seq 
    return new_stream, spec_name, spec_id
