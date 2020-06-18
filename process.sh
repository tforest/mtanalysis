#!/bin/bash

# This script aims to run an overall analysis of genomes using orthologous groups construction 
# It takes genbank files in input with extension 'gbf' by default.

# ARGS
GENOMES_PATH=$1

# PROGRAMS PATHS
#GBK_TO_FAA="$HOME/Documents/Tools/gbk_to_faa.py
GBK_TO_FAA="SCRIPTS/gbk_to_faa.py"
### CONSTANTS

GENBANK_EXTENSION=".gbf"

#Eurotiales (Microsporum canis)
#SELECTED_TAXO_GRP="eurotiales;NC_012832"

# BOLETALES
#SELECTED_TAXO_GRP="boletales;Aspbisp1"


#Sordariales
SELECTED_TAXO_GRP="sordariales;CodFL1790_1;NC_027509"

# OUTPUT
#OUT_FOLD="EUROTIALES_2"
#OUT_FOLD="BOLETALES_2"
#OUT_FOLD="SORDARIALES_2"
OUT_FOLD="SORDARIALES_3"

########### MAIN

## Convert GBK to fasta files
for filename in $GENOMES_PATH/*$GENBANK_EXTENSION; do
    echo "python3 ~/Documents/Tools/gbk_to_faa.py $filename"
    python3 $GBK_TO_FAA $filename
done

mkdir $OUT_FOLD
mkdir $OUT_FOLD/PROT_FASTAS
mv $GENOMES_PATH/*_converted.faa $OUT_FOLD/PROT_FASTAS/


mtanalysis -postannot $OUT_FOLD/PROT_FASTAS/* -ids DATA/metadata.csv
mkdir $OUT_FOLD/RAW_PROTEOME
mv $OUT_FOLD/PROT_FASTAS/* $OUT_FOLD/RAW_PROTEOME/
mv *.fasta $OUT_FOLD/PROT_FASTAS/
rm $OUT_FOLD/PROT_FASTAS/_.fasta

mtanalysis -selectgroup $SELECTED_TAXO_GRP DATA/metadata.csv $OUT_FOLD/PROT_FASTAS/*.fasta
find $OUT_FOLD/PROT_FASTAS/*.ret -size  0 -print -delete

mkdir $OUT_FOLD/REDUCED_PROTEOMES

mv $OUT_FOLD/PROT_FASTAS/*.ret $OUT_FOLD/REDUCED_PROTEOMES

mkdir $OUT_FOLD/ORTHOGROUPS/
rename 's/(.*)/$1.fasta/' $OUT_FOLD/REDUCED_PROTEOMES/*.ret



# ORTHOFINDER
orthofinder -f $OUT_FOLD/REDUCED_PROTEOMES

mkdir $OUT_FOLD/ORTHOGROUPS/ORTHOFINDER

mv $OUT_FOLD/REDUCED_PROTEOMES/Ortho*/Res*/Orthogroup_Sequences/* $OUT_FOLD/ORTHOGROUPS/ORTHOFINDER

mkdir $OUT_FOLD/SILIX
cat $OUT_FOLD/REDUCED_PROTEOMES/*.fasta > $OUT_FOLD/SILIX/CONCAT.fasta
mkdir $OUT_FOLD/SILIX/BLAST
makeblastdb -in $OUT_FOLD/SILIX/CONCAT.fasta -dbtype prot -out $OUT_FOLD/SILIX/BLAST/PROTEOME_DB

blastp -db $OUT_FOLD/SILIX/BLAST/PROTEOME_DB -query $OUT_FOLD/SILIX/CONCAT.fasta -outfmt 6 -out $OUT_FOLD/SILIX/BLAST/blastall.out -num_threads 4

silix $OUT_FOLD/SILIX/CONCAT.fasta $OUT_FOLD/SILIX/BLAST/blastall.out -n > $OUT_FOLD/SILIX/node.silix
# HIFIX
mkdir $OUT_FOLD/HIFIX
hifix -t 4 $OUT_FOLD/SILIX/CONCAT.fasta $OUT_FOLD/SILIX/BLAST/blastall.net $OUT_FOLD/SILIX/node.silix > $OUT_FOLD/HIFIX/clusters.hfx

mkdir $OUT_FOLD/ORTHOGROUPS/HIFIX

mtanalysis -postannot $OUT_FOLD/REDUCED_PROTEOMES/*.fasta -clusters $OUT_FOLD/HIFIX/clusters.hfx

mv OG* $OUT_FOLD/ORTHOGROUPS/HIFIX

mkdir $OUT_FOLD/ORTHOGROUPS/INTERSECT

mtanalysis -intersect $OUT_FOLD/ORTHOGROUPS/ORTHOFINDER/* $OUT_FOLD/ORTHOGROUPS/HIFIX/*

mv GROUP_* $OUT_FOLD/ORTHOGROUPS/INTERSECT 

find $OUT_FOLD/ORTHOGROUPS/INTERSECT/* -size 0 -print -delete

for filename in $OUT_FOLD/ORTHOGROUPS/INTERSECT/*; do
    echo "Building tree for $(basename $filename)"
    xvfb-run ete3 build -w standard_raxml -a $filename -o $OUT_FOLD/ETE/INTERSECT/$(basename $filename) --tools-dir ~/.etetoolkit/ext_apps-latest/ -C4 > $(basename $filename).log
    sleep 5s
done

mv GROUP_*.log $OUT_FOLD/ETE/INTERSECT/

mkdir $OUT_FOLD/TREES
find $OUT_FOLD/ETE/ -type f -name "*.nw" -exec cp {} $OUT_FOLD/TREES \;

for filename in $OUT_FOLD/TREES/*.nw; do
    echo "Reformat Tree $(basename $filename)"
    SCRIPTS/reformat_tree.py $filename
done

#concat all trees 
for f in $OUT_FOLD/TREES/*.mod.tre;do cat $f; echo; done> $OUT_FOLD/TREES/TREES.tre
# add the unrooted tag for duptree
sed -i 's/^/[\&U]/' $OUT_FOLD/TREES/TREES.tre


#Building consensus tree
./SCRIPTS/duptree/duptree -i $OUT_FOLD/TREES/TREES.tre -o $OUT_FOLD/TREES/CONSENSUS.out --nogenetree

#Serialize for Kmedoid Tree Clustering
for filename in $OUT_FOLD/TREES/*.nw; do
    echo "Reformat Tree $(basename $filename)"
    SCRIPTS/reformat_tree.py -s $filename 
done

#concat all trees 
for f in $OUT_FOLD/TREES/*.mod.tre;do cat $f; echo; done> $OUT_FOLD/TREES/TREES_ser.tre

cat $OUT_FOLD/TREES/*.ser > $OUT_FOLD/TREES/TREES_ser.tre.ser

# Clustering phase

KMTC -tree $OUT_FOLD/TREES/TREES_ser.tre 3

mkdir $OUT_FOLD/CLUSTERING

mv output.txt ' input_.txt' stat.csv $OUT_FOLD/CLUSTERING/

SCRIPTS/KMTC_parser.py $OUT_FOLD/CLUSTERING/output.txt

mv CLUST_* $OUT_FOLD/CLUSTERING/

#Reconstruct Trees with original leaf names

for cluster in $OUT_FOLD/CLUSTERING/CLUST_*; do
   for tree in $cluster/*; do
       echo "Unserialize $(basename $tree)"
       SCRIPTS/reformat_tree.py --unserialize $tree --serfile $OUT_FOLD/TREES/TREES_ser.tre.ser -r
   done
done
