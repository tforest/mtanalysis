from Bio import Entrez
import sys
Entrez.email = "Your.Name.Here@example.org"

with open(sys.argv[1], "r") as inputf:
    for ids in inputf.readlines():
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
        the_id = ids.strip()
        with open(the_id+".gb", "w") as f:
            for line in handle.readlines():
                f.write(line)
        handle.close()
