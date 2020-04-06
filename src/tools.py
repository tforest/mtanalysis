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

def split_clusters(fastas, clusters, ids=None, filter=0.1):
    clusts = load_clusters(clusters)
    all_fastas = {}
    print("Number of species fasta files =", len(fastas),
          "\nNumber of clusters =", len(clusts))
    #max_key, max_value = max(clusts.items(), key = lambda x: len(set(x[1])))
    if ids:
        # add a dict with IDS-name corresp
        ids_names = add_name_from_nc_id(ids)
    for fas in fastas:        
        fasta_stream = fasta.parse(fas)
        all_fastas.update(fasta_stream)
    for i, clust in enumerate(clusts):
        clustname = "OG"+format(int(clust.split("_")[-1]), '06d')
        print("[", round(i/len(clusts)*100),"% ]", "Cluster",
              clustname, "size =", len(clusts[clust]))
        if(len(clusts[clust])>len(fastas)-len(fastas)*0.1):
            print("Keeping", clustname, len(clusts[clust]), ">", len(fastas)-len(fastas)*filter)
            with open(clustname+'.fa', 'w') as clust_out:
                for key in clusts[clust]:
                    for fasta_key in all_fastas.keys():
                        if key in fasta_key:
                            clust_out.write(">"+fasta_key+"\n"+all_fastas[fasta_key]+"\n")
    #print(all_fastas.keys())
        # with open(fas+'.clust', 'w') as fas_out:
        #     for key in fasta_stream:
        #         #print('_'.join(key.split('_')[-3:-1]))
        #         fas_out.write(">"+key+'\n')
        #         fas_out.write(fasta_stream[key]+'\n')
                    
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
def rename_mfannot_proteome(fastas, ids_file, column_id = -1, column_name = -3):
    ids = []
    names = {}
    with open(ids_file, 'r') as idsinput:
        for line in idsinput.readlines():
            id = line.strip().split(',')[column_id].strip()
            name = line.strip().split(',')[column_name].strip()
            if name == "":
                name = "UNKNOWN"
            names[id] = name
    for fas in fastas:        
        fasta_stream = fasta.parse(fas)
        # update fasta headers with new names
        fasta_stream, spec_name, spec_id = fasta.set_header(fasta_stream,
                                                            names)
        new_filename = spec_name+"_"+spec_id+".fasta"
        new_filename = new_filename.replace("/", "")
        with open(new_filename, "w") as out:
            print("Generating", new_filename)
            for head, seq in fasta_stream.items():
                out.write(">"+head+"\n"+seq+"\n")

def ortho_intersect(fastas):
    group = {}
    all_fastas = {}
    # parse fastas and build all groups
    for k, fas in enumerate(fastas):        
        fasta_stream = fasta.parse(fas)
        all_fastas.update(fasta_stream)
        group[k] = []
        for key in fasta_stream.keys():
            if len(key.split("~")) != 3:
                new_head = [key.split("_")[0], "_".join(key.split("_")[1].split('-')), '_'.join(key.split("_")[2:])]
                new_head = '~'.join(new_head)
                group[k].append(new_head)
            else:
                group[k].append(key)
    new_group = {}
    for grp_a in group:
        new_group[grp_a] = []
        for grp_b in group:
            if grp_b == grp_a:
                break
            for key_a in group[grp_a]:
                if key_a in group[grp_b]:
                    new_group[grp_a].append(key_a)
        print(round(grp_a/len(group)*100), "%")
    print("-----COMP-----")        
    # for k in new_group:
    #     print(len(new_group[k]))
        
    for i, k in enumerate(sorted(new_group, key=lambda k: len(new_group[k]), reverse=True)):
        if i>20:
            break
        print(k, len(new_group[k]))
        with open("GROUP_"+str(i), "w") as output:    
            for entry in new_group[k]:
                if entry.startswith("orf"):
                    continue
                else:
                    output.write(">"+entry+"\n"+all_fastas[entry]+"\n")
