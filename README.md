# TIGER: Targeted Integrative Genetic Element Retriever

**Version:** 2  
**Authors:** Catherine Mageeney, Gary Trubl, Kelly Williams


## Table of Contents
- [Citations](#Citations)
- [Software Description](#Description)
- [Software Dependencies](#Dependencies)
- [Installation](#Installation)
- [Usage Guide](#Usage)
- [Contact](#contact)


## Citations
If you're using this software in a publication, please cite:
1. Mageeney CM, Trubl G, Williams KP. 2022. Improved mobilome delineation in fragmented genomes. *Front Bioinform* doi: [10.3389/fbinf.2022.866850](https://doi.org/10.3389/fbinf.2022.866850)
2. Mageeney CM, Lau BY, Wagner JW, Hudson CM, Schoeniger JS, Krishnakumar R, Williams KP. 2020. New candidates for regulated gene integrity revealed through precise mapping of integrative genetic elements. *Nucleic Acids Research* 48(8):4052-4065 (doi: [10.1093/nar/gkaa156](https://doi.org/10.1093/nar/gkaa156))
3. Hudson CM, Lau BY, Williams KP. 2015. Islander: a database of precisely mapped genomic islands in tRNA and tmRNA genes. *Nucleic Acids Research* V43(D1):D48–D53 (doi: [10.1093/nar/gku1072](https://doi.org/10.1093/nar/gku1072))

## Software Description
Within this repository are two programs avilable for use: Islander and Tiger. These two programs use different methodology for identifiying putative integrative genetic elements in genomic contigs.
### Islander
To identify genomic islands the Islander software looks for tyrosine integrase genes located within either a tmRNA or tRNA. It first identifies tRNA and tmRNAs using tRNAscan-SE (tRNA), BRUCE (tmRNA), and ARAGORN (both). It then identifies nearby integrases (excluding Xer and integron subclasses) using a integrase specific HMM with HMMER3. Using the sequence of the identified tRNA or tmRNA with a nearby integrase, a Blast search is then conducted to search for the cognate end which corresponds to the end of the putative integrated genomic island. 
<p align="center">
  <img src=image-1.png />
</p>

<sub>Figure 1. Islander algorithm. (A and B) Population Phase: tRNA and tmRNA genes (tDNAs), tDNA fragments and integrase genes are placed on the chromosome, and each interval between a tDNA and its cognate fragments is considered a candidate island. (C) Filtering Phase: Candidates pass through several filters, including tests for an integrase gene, correct fragment/tDNA orientation and length. (D) Resolution Phase: Multiple candidates at the same tDNA are resolved, identifying tandem arrays when each island in the array has its own tDNA fragment and integrase gene. (c) Hudson CM et. al (2015)<sub>


## Software Dependencies
TIGER requires the following programs to be available as system-wide executables and has been tested with the following software versions. We reccomend installing these packages and their dependencies using conda (instructions: [Installation](#installation)).

- Prokka v1.11 (https://github.com/tseemann/prokka)
- PfTools v3 (https://github.com/sib-swiss/pftools3)
- tRNAscan-SE 2.0 (https://github.com/UCSC-LoweLab/tRNAscan-SE)
- Perl IPC::Run3 

## Installation
This tutorial recomends that you have a working version of Anaconda or Miniconda (Miniconda is suggested: https://docs.anaconda.com/miniconda/) to properly install all dependencies.

```bash
#Create a conda environment for Tiger
conda create --name Tiger

#load your new environment
conda activate Tiger

# install prokka. Note, prokka has a known bug with conda that you must address by editing its code 
# base. See instructions below.

    conda install -c conda-forge -c bioconda -c defaults prokka

#setup/update prokka databases

    prokka --setupdb

#Install pftools

    conda install bioconda::pftools

#Install tRNA-Scan

    conda install bioconda::trnascan-se

#Install Perl modules

    conda install bioconda::perl-ipc-run3

#Install Git (if not already available on your system)

    conda install anaconda::git


# download Tiger repository from GitHub
    git clone -b TIGER https://github.com/sandialabs/TIGER.git

# Optional: We recomend adding the path of TIGER/bin to your bash profile to avoid having to call the
# full program during use. We always advise making a copy of your bash profile before making any 
# edits as a typo may have severe consequences.

# ensure the execultables in the TIGER repository have executable permission on your system
    
    cd <path to TIGER folder>/bin
    chmod +x *

#install wget if not available on your system
    
    conda install anaconda::wget

# get PFAM_A hmm database from NCBI: this is an OLD PFAM database (v35.0) if you want another 
# version, you can change it in the ftp path below (note, this is the version that is suggested to
# use with TIGER per https://doi.org/10.3389/fbinf.2022.866850). This database must be installed to 
# "Path to Tiger Folder"/db as PFAM-A.hmm

    wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz -o <path to TIGER folder>/db/Pfam-A.hmm.gz

# install gzip if not already available on your system
    
    conda install conda-forge::gzip

# uncompress PFAM database:

    gunzip <path to TIGER folder>/db/Pfam-A.hmm.gz

#IMPORTANT: Fix a bug in Prokka:
# I'm not an author of this program but I've found this bug reported in several forums without an
# easy solve given. What's happening is that the prokka main script is unable to parse the version
# number of bioperl that it retrieves with conda correctly (giving the error: Argument "1.7.8" isn't 
# numeric in numeric lt (<)) and thus it exits the program immediately. My "quick and ugly" fix is to
# hash out the lines in their code that are doing this check. So, in the prokka script (this script 
# should be found at "path to your conda distribution"/envs/"conda env name"/bin/prokka) add a hash
# (#) to the beginning of lines 256-259 like so:

    "#my $minbpver = "1.006002"; # for Bio::SearchIO::hmmer3"
    "#my $bpver = $Bio::Root::Version::VERSION;"
    "#msg("You have BioPerl $bpver");"
    "#err("Please install BioPerl $minbpver or higher") if $bpver < #$minbpver;"

# We recommend using either nano or vim to make these changes

```


## Usage Guide


# TIGER Requirments:



![alt text](image.png)

Requires a .fa file in the same folder with the same prefix.
You may add a .tax file to specify what the genome is.
The tab-separated fields of the one-line .tax file are: 
 1. taxid
 2. organism
 3. Division;Phylum;Class;Order;Family;Genus;Species
 4. Genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
 5. Nickname (short name for organism, eg Eco837 for the 837th E. coli genome)

## Sample calls to try within /testdata (PATH: to TIGER installation; DB: to reference genome blast database)

```perl PATH/bin/islander.pl -verbose genome.fa &> islander.log```

```perl PATH/bin/tiger.pl -verbose -db DB -cross simple -fasta genome.fa &> tiger.log```

```perl PATH/bin/resolve.pl mixed lenient genome genome.island.merge.gff islander.gff &> resolve.log```

```perl PATH/bin/typing.pl resolved.gff &> typing.log```

```perl PATH/bin/typing.pl genome.island.nonoverlap.gff &> typing.log```

## Notes:
before rerunning islander.pl: ```rm genome.stats```

before rerunning tiger.pl: ```rm genome.island.nonoverlap.gff```

## Contact
While we're working on a more complete usage guide, please reach out to eltorra@sandia.gov for help if this software is necessary for your research