# TIGER version 2
Target / Integrative Genetic Element Retriever, version 2

# Citation
Mageeney CM, Lau BY, Wagner JW, Hudson CM, Schoeniger JS, Krishnakumar R and Williams KP. 2020. New candidates for regulated gene integrity revealed through precise mapping of integrative genetic elements. Nucleic Acids Research 48(8):4052-4065 (doi.org/10.1093/nar/gkaa156)
Mageeney CM, Trubl G, Williams KP. 2022. Improved mobilome delineation in fragmented genomes. Front Bioinform doi: 10.3389/fbinf.2022.866850

# INSTALLATION
git clone -b TIGER2 https://github.com/sandialabs/TIGER.git

User should have a reference genome blast database available, such as refseq_genomic

User should download Pfam-A.hmm from pfam and place it or a symbolic link to it in the TIGERPATH/db directory

```ln -s /ABSOLUTE_PATH/Pfam-A.hmm TIGERPATH/db/Pfam-A.hmm```

The following programs (with suggested versions) must be properly installed and in the user's path:
* blastn, blastdbcmd, makeblastdb 2.6.0+
* prokka 1.11
* bedtools 2.27.1
* tRNAscan-SE 2.0.2
* cmscan 1.1.2
* Perl Core: List::Util, File::Spec, Cwd, Getopt::Long
* Perl Noncore: IPC::Run3

# RUNNING
Requires a .fa file and .tax file in the same folder with the same prefix, as at TIGERPATH/testdata
The tab-separated fields of the one-line .tax file are: 
 1. taxid
 2. organism
 3. Division;Phylum;Class;Order;Family;Genus;Species
 4. Genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
 5. Nickname (short name for organism, eg Eco837 for the 837th E. coli genome)

## Sample calls to try within /testdata (PATH: to TIGER installation; DB: to reference genome blast database)
```perl PATH/bin/islander.pl -verbose genome.fa &> islander.log```

```perl PATH/bin/tiger.pl -verbose -db DB -fasta genome.fa &> tiger.log```

```perl PATH/bin/typing.pl genome.island.nonoverlap.gff &> typing.log```

```perl PATH/bin/resolve.pl mixed > resolved.gff 2> resolved.log```  # Not yet updated for use with TIGER2

```perl PATH/bin/typing.pl resolved.gff &> typing.log```  # Use only when resolve.pl is updated

## Notes:
before rerunning islander.pl: ```rm genome.stats```

before rerunning tiger.pl: ```rm genome.island.nonoverlap.gff```
