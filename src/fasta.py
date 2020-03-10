DEBUG = False

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
