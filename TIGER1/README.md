# TIGER
Target / Integrative Genetic Element Retriever

# Citation
Mageeney CM, Lau BY, Wagner JW, Hudson CM, Schoeniger JS, Krishnakumar R and Williams KP. 2000. New candidates for regulated gene integrity revealed through precise mapping of integrative genetic elements. bioRxiv 2020.01.24.918748 (doi:10.1101/2020.01.24.918748)

# INSTALLATION
git clone https://github.com/sandialabs/TIGER.git

User should have a reference genome blast database available, such as refseq_genomic

User should download Pfam-A.hmm from pfam (we currently use the version from July 2010) and place it or a symbolic link to it in the comparator/db directory

```ln -s /ABSOLUTE_PATH/Pfam-A.hmm comparator/db/Pfam-A.hmm```

The following programs (with suggested versions) must be properly installed and in the user's path:
* blastn, blastdbcmd, makeblastdb 2.6.0+
* prokka 1.11
* bedtools 2.27.1
* tRNAscan-SE 2.0.2
* cmscan 1.1.2
* Perl Core: List::Util, File::Spec, Cwd, Getopt::Long
* Perl Noncore: IPC::Run3

# RUNNING
Requires a .fa file and .tax file in the same folder with the same prefix, as at comparator/testdata
The tab-separated fields of the one-line .tax file are: 
 1. taxid
 2. organism
 3. Division;Phylum;Class;Order;Family;Genus;Species
 4. Genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
 5. Nickname (short name for organism, eg Eco837 for the 837th E. coli genome)

## Sample calls to try within /testdata (PATH: to TIGER installation; DB: to reference genome blast database)
```perl PATH/islander.pl -verbose genome.fa &> islander.log```

```perl /projects/islands/tiger/TIGER/tiger.pl -verbose -db DB -fasta genome.fa &> tiger.log```

```perl PATH/resolve.pl mixed > resolved.gff 2> resolved.log```

```perl PATH/bin/typing.pl resolved.gff &> typing.log```

## Notes:
before rerunning islander.pl: ```rm genome.stats```

before rerunning comparator.pl: ```rm genome.island.nonoverlap.gff```

