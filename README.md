# TIGER: Targeted Integrative Genetic Element Retriever

**Version:** 2  
**Authors:** Catherine Mageeney, Gary Trubl, Kelly Williams


## Table of Contents
- [Citations](#Citations)
- [Software Description](#Description)
- [Software Dependencies](#Dependencies)
- [Installation](#Installation)
- [Usage Guide](#Guide)
- [Contact](#contact)


## Citations
If you're using this software in a publication, please cite:
1. Mageeney CM, Trubl G, Williams KP. 2022. Improved mobilome delineation in fragmented genomes. *Front Bioinform* doi: [10.3389/fbinf.2022.866850](https://doi.org/10.3389/fbinf.2022.866850)
2. Mageeney CM, Lau BY, Wagner JW, Hudson CM, Schoeniger JS, Krishnakumar R, Williams KP. 2020. New candidates for regulated gene integrity revealed through precise mapping of integrative genetic elements. *Nucleic Acids Research* 48(8):4052-4065 (doi: [10.1093/nar/gkaa156](https://doi.org/10.1093/nar/gkaa156))
3. Hudson CM, Lau BY, Williams KP. 2015. Islander: a database of precisely mapped genomic islands in tRNA and tmRNA genes. *Nucleic Acids Research* V43(D1):D48–D53 (doi: [10.1093/nar/gku1072](https://doi.org/10.1093/nar/gku1072))

## Software Description
Within this repository are two programs available for use: Islander and TIGER. These two programs use different methodology for identifiying putative integrative genetic elements in genomic contigs.
### Islander
To identify genomic islands the Islander software looks for tyrosine integrase genes located within either a tmRNA or tRNA. It first identifies tRNA and tmRNAs using tRNAscan-SE (tRNA), BRUCE (tmRNA), and ARAGORN (both). It then identifies nearby integrases (excluding Xer and integron subclasses) using a integrase specific HMM with HMMER3. Using the sequence of the identified tRNA or tmRNA with a nearby integrase, a BLAST search is then conducted to search for the cognate end which corresponds to the end of the putative integrated genomic island. All candidate islands are then subjected to a series of tests to determine whether they represent true integrate genomic islands (see image below).
<p align="center">
  <img src=image-1.png />
</p>

<sub> Graphical Description of the Islander algorithm. (A and B) Population Phase: tRNA and tmRNA genes (tDNAs), tDNA fragments and integrase genes are placed on the chromosome, and each interval between a tDNA and its cognate fragments is considered a candidate island. (C) Filtering Phase: Candidates pass through several filters, including tests for an integrase gene, correct fragment/tDNA orientation and length. (D) Resolution Phase: Multiple candidates at the same tDNA are resolved, identifying tandem arrays when each island in the array has its own tDNA fragment and integrase gene. (Image (c) CM Hudson et. al, 2015)<sub>

### TIGER and TIGER2
TIGER utilizes the output of Islander for locational data regarding the position of integrase genes including serine integrases. It then uses the sequences 15kb to the left and right of the integrase gene midpoint to query against a BLAST database of species-specific reference genomes to search for query sequences in close proximity which signify the absence of an integrated genetic island and an aproximate attB sequence or, the site of putative island integration. 
<p align="center">
  <img src=image-2.png />
</p>

<sub> Graphical description of TIGER's usage of ping-pong BLAST for integrated genetic element (IGE) discovery. The corresponding regions of an IGE-bearing and uninterrupted reference genome pair (A) produce a sequence alignment pattern (B). Strand crossover presumably occurs with the direct repeat block (yellow). In TIGER (C), the first BLAST simultaneously locates the int-proximal end of the IGE and the attB, and the second locates the distal end of the IGE. (Image (c) CM Mageeney et. al, 2020)<sub>

TIGER2 improves on TIGER by enabling integrated genomic element identification in cases where the whole putative island sequence is not present on a singular scaffold/contig. Ultimately, the TIGER2 update introduces two new “split” modes that yield split GIs, in addition to the intact GIs. “CircleOrigin” mode finds split GIs that wrap around the origin of a circular replicon. “Cross” mode detects split GIs with termini on separate scaffolds. 
<p align="center">
  <img src=image-3.png />
</p>

<sub> Graphical Demonstration of Fragmented Islands which Can be Identified with Tiger2. The same circular chromosome with 3 (colored) GIs is shown with a complete (A,B) or fragmented assembly (C). With complete assembly, if the origin of the linearized sequence of the circle is randomly chosen, it will occasionally fall within a GI, splitting the GI (B). Yields are shown for the various TIGER modes. The original mode can only find intact GIs on a single scaffold, while the new modes, CircleOrigin (applied to complete assemblies) and Cross (applied to fragmented assemblies), can additionally find the split islands. Because TIGER focuses on GI-flanking sequences, the Cross-mode call for a multiply split GI (red in panel C) will only include the terminal fragments and exclude middle GI fragments. (Image (c) CM Mageeney et. al, 2022)<sub>

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

    git clone https://github.com/sandialabs/TIGER.git

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

### Islander Requirements:
 The only requirment of Islander is an output directory containing a single genomic file in fasta format with a '.fa' subscript.

 ***
 *Islander Flags:*
 - '-outDir':    Output directory. Default: same directory as GENOME_FASTA_FILE.
 - '-tax':       Taxonomic info for query genome. Enter file in outDir containing 
               NCBI taxonomy string with a '.tax' subscript, or use B for Bacteria, 
               A for Archaea, M for Mycoplasmatales/Entomoplasmatales, G for 
               Gracilibacteria/candidate division SR1. Automatically sets -gencode. 
               Default: B.
 - '-gencode':   Genetic code table (see NCBI). Default: 11
 - '-nickname':  Brief name for genome (as might be used to start a locus_tag).
 - '-criterion': Basis for overlap resolution, 3 options: random, score (7-test
               false positive formula), deltaGC. Default = score.
 - '-virus':     Comma-separated list of entries assigned as viruses.
 - '-complete':  Consider genome complete and categorize replicons. Default:
                consider genome incomplete and call all entries contigs.
 - '-force':     Overwrite current output files. Default: leave existing files.
 - '-cpu':       Number of cpus to use. Default: $cpu.
 - '-tateronly': Toggle to exit after running tater.pl annotator. Default: off.
 - Additional options: '-help', '-version', '-verbose', '-authors', '-license'
***

Executing Islander Example:

From your Conda Environment:
```bash
conda activate Tiger
```

```bash
cd <path to out directory containing singular .fa file> ; islander.pl -verbose <fasta file name>
```
<sub> Note: if you made the "path to TIGER"/bin a system-wide executable by enabling execute with chmod and adding the path to your bash profile, you will not need to specify perl (perl "path to TIGER"/bin/islander.pl) or use the path to the script in your launch command. If not, you will need to call the program as: cd "path to out directory containing singular .fa file" ; perl "path to TIGER"/bin/islander.pl -verbose <fasta file name><sub>

To rerun this program: delete genome.stats

### TIGER Requirements:
TIGER requires a singular genomic or metagenomic file in fasta format ('.fa' subscript) in an output directory and a BLAST database of species specific reference genomes. See https://github.com/sandialabs/SmartDBs for instructions on generating this database. Alternatively, if you only need a BLAST database for a single/few species you can reach out to Kelly Williams (kpwilli at sandia dot gov) or Katie Mageeney (cmmagee at sandia dot gov) for assisance. Please provide the GTDB taxonomic name for your species as well as the link to either a DropBox or GoogleDrive Folder for the database to be deposited in.

***
*TIGER Flags:*
 - '-fasta':    Genomic fasta DNA sequence file.
 - '-db':   Blast database of reference genomes, absolute path.
 - '-search':   Search type. Specify island or IS. Default: island.
 - '-tax':  Taxonomic info for query genome. Enter name of a file containing 
    NCBI taxonomy string, or use B for Bacteria, A for Archaea, M for 
    Mycoplasmatales/Entomoplasmatales, G for Gracilibacteria/candidate
    division SR1. Automatically sets -gencode. Default: B.
 - '-gencode':  Genetic code table to use (see NCBI). Default: 11.
 - '-nickname': Brief name for genome (as might be used to start a locus_tag).
 - '-circle':   Specify C if all genomic DNA sequences are circular, L if all DNAs 
    are linear, or a filename for a tab-delimited file of query and 
    circularity (eg. acc.vers[tab]circular/linear). Default: C.
 - '-cross':    Three options: intact, cross, or circleOrigin. Default: intact. 
 - '-complete': Consider genome complete and categorize replicons. Default:
             consider genome incomplete and call all entries contigs.
 - '-qlen':     Query length for islands. Default: 15000 (3000 is always used
                in the second pass test for IS's that rule out island artifacts).
 - '-force':    Overwrite current output files. Default: leave existing files.
 - '-outDir':   Output directory. Default: same directory as GENOME_FASTA_FILE.
 - '-cpu':      Number of cpus to use. Default: 1.
 - Additional options: '-help', '-version', '-verbose', '-authors', '-license'
***

Executing TIGER Example:

From your Conda Environment:
```bash
conda activate Tiger
```

```bash
cd <path to out directory containing singular .fa file> ; tiger.pl -verbose -db <path to reference genome database and database prefix> -fasta <fasta file name>
```
<sub> Note: if you made the "path to TIGER"/bin a system-wide executable by enabling execute with chmod and adding the path to your bash profile, you will not need to specify perl (perl "path to TIGER"/bin/tiger.pl) or use the path to the script in your launch command. If not, you will need to call the program as: cd "path to out directory containing singular .fa file" ; perl "path to TIGER"/bin/tiger.pl -verbose -db "path to reference genome database and database prefix" -fasta "fasta file name".

To rerun this program, delete 'genome.island.nonoverlap.gff'

## Contact
Please leave an issue on this GitHub repo or reach out to Kelly Williams (kpwilli at sandia dot gov) or Katie Mageeney (cmmagee at sandia dot gov) for questions or assistance.

