import mtanalysis
import fasta
from Bio import Entrez

def add_name_from_nc_id(ids_file):
    Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are

    ids = []
    names = {}
    with open(ids_file, 'r') as idsinput:
        for line in idsinput.readlines():
            ids.append(line.strip())
    with open(ids_file) as ids_input:
        ids_lines = ids_input.readlines()
        for i, id in enumerate(ids_lines):
            try:
                handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
                print(i, len(ids_lines), i/len(ids_lines)*100)
                names[id] = (handle.readlines()[0].split(',')[0])
            except:
                print("This entry was not found in Nucleotide DB:", id)
    return names
def tokenize(files_input):
    files = mtanalysis.ambiguous_open(input_source=files_input)
def post_annot(files):
    print(files)
def load_clusters(clusters):
    clusts = {}
    for cluster in clusters:
        with open(cluster) as clust:
            for line in clust.readlines():
                id = line.split()[0]
                val = line.split()[1]
                if id in clusts.keys():
                    clusts[id].append(val)
                else:
                    clusts[id] = [val]
    return clusts
def remove_spurious_clusters(fastas, clusters, ids=None):
    clusts = load_clusters(clusters)
    max_key, max_value = max(clusts.items(), key = lambda x: len(set(x[1])))
    if ids:
        # add a dict with IDS-name corresp
        ids_names = add_name_from_nc_id(ids)
    for fas in fastas:        
        fasta_stream = fasta.parse(fas)
        with open(fas+'.fasta', 'w') as fas_out:
            for key in fasta_stream:
                print('_'.join(key.split('_')[-3:-1]))
                if key in max_value:
                    fas_out.write(">"+key+'\n')
                    fas_out.write(fasta_stream[key]+'\n')

def detect_bad_genomes(fastas):
    duplicates = []
    species_list = []
    genes = {}
    for fas in fastas:        
        fasta_stream = fasta.parse(fas)
        for gene in list(fasta_stream.keys()):
            #print(gene.split('_')[0], gene.split('_')[1])
            gene_name = gene.split('_')[0]
            if len(gene.split('_')[1]) == 1:
                species_id = gene.split('_')[2]
            else:
                species_id = gene.split('_')[1]
            print(species_id, gene.split('_'))
            if gene_name in genes.keys():
                genes[gene_name].append(species_id)
            else:
                genes[gene_name] = [species_id]
            
    for gene in genes:            
            for species in genes[gene]:
                if species not in species_list:
                    species_list.append(species)
                if genes[gene].count(species) > 1:
                    #print(gene, species)
                    #if gene in ['cox1', 'cox2', 'cob']:
                    if species not in duplicates:
                        duplicates.append(species)
    print(duplicates, len(duplicates), len(species_list))
    for species in species_list:
        if species not in duplicates:
            #print(species)
            continue
