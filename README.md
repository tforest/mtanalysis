# MTanalysis

A custom toolkit to analyse fungal mitogenomes.

### Prerequisites

- argparse
- numpy
- matplotlib
- biopython

### Running

You should be ready to run the program, by calling, for example
```
python3 mtanalysis.py

### Help

  -h, --help            show this help message and exit
  -a                    Perform analysis of genomes
  -f [F [F ...]]        input file(s)
  -r R                  repertory with input files
  -ids IDS              ids file for annotation (usually in CSV format)
  -rpsb RPSB            parse RPS-BLAST results got from rpsbproc
  -comp COMP            Compare FASTA files to Domains identified by rpsb
  -tokenize [TOKENIZE [TOKENIZE ...]]
                        Tokenize fasta headers
  -postannot [POSTANNOT [POSTANNOT ...]]
                        Post annotation treatment (needs a list of fastafiles to be parsed)
  -headermode [HEADERMODE [HEADERMODE ...]]
                        Specify header mode for postannot(0, 1 or 2)
  -clusters [CLUSTERS [CLUSTERS ...]]
                        Clusters of genes from fasta files
  -intersect [INTERSECT [INTERSECT ...]]
                        Get intersection of orthogroups given a set of orthologous sequences groups as fasta files
  -selectgroup [SELECTGROUP [SELECTGROUP ...]]
                        Select specific species from taxonomic characteristics
  -concat [CONCAT [CONCAT ...]]
                        Concatenate all fasta sequences that have all genes in common in a group of files
  -renamefromid [RENAMEFROMID [RENAMEFROMID ...]]
                        Rename fasta header accordingly to ID spec in taxonomy file


## Author

* **FOREST Thomas** - *M2BI* - [Univ-Paris-Diderot.fr](https://www.univ-paris-diderot.fr/) - thomas.forest@etu.univ-paris-diderot.fr

## License

This project is licensed under the CeCILL-C License - see [eCILL-B FREE SOFTWARE LICENSE AGREEMENT](http://cecill.info/licences/Licence_CeCILL-B_V1-en.html) for details
