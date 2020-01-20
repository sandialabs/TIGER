# TIGER
The Integrative Genetic Elements Retriever: Regulated gene integrity revealed through precise mapping of genomic islands

# INSTALLATION
git clone https://github.com/sandialabs/TIGER.git

User should have a reference genome blast database available, such as refseq_genomic

User should download Pfam-A.hmm from pfam (we currently use the version from July 2010) and place it or a symbolic link to it in the comparator/db directory

```ln -s /ABSOLUTE_PATH/Pfam-A.hmm comparator/db/Pfam-A.hmm```

The following programs (with suggested versions) must properly installed and in the user's path:
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

## Sample calls to try with /testdata
```perl PATH/islander.pl -verbose genome.fa```

```perl PATH/comparator.pl -verbose -db /databases/blast/refseq_genomic -fasta genome.fa  # switch in your own path to refseq_genomic```

```perl PATH/resolve.pl mixed```

## Notes:
before rerunning islander.pl: ```rm genome.stats```

before rerunning comparator.pl: ```rm genome.island.nonoverlap.gff```

