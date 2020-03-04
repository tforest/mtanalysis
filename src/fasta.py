DEBUG = False

def parse(filename, domains, extension='.asn.out'):
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
