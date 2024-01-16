#!/usr/bin/perl
use strict; use warnings;
use File::Spec;

# OPTIONS
my $dir = File::Spec->rel2abs($0); $dir =~ s/\/([^\/]+)$//; my $scriptname = $1;
my ($verbose, $nickname, $force, $tax, $gencode, $extraDoms) = ('', '', '', '', '', '');
my $circ = 'C';
my $cpu = 1;
my $hmm = "$dir/db/Pfam-A.hmm";
Options(); # see bottom of script; help and other messages, Getopt
my ($dna, $dnapath) = (File::Spec->rel2abs($ARGV[0]), '.'); $dnapath = $1 if $dna =~ s/(.*\/)//; $dna =~ s/\.fa//;
my @doms; for (split ',', $extraDoms) {push @doms, File::Spec->rel2abs($_)}; $extraDoms = join(',', @doms);
$extraDoms = "-extraDoms $extraDoms" if $extraDoms;

# PROGRAM
chdir $dnapath;
die "Can't find $dnapath/$dna.fa file\n" unless -e "$dna.fa";
Tax();
FindGenes(); # find tRNAs, integrases
CollectGff();

sub Tax {
 if ($gencode and ($gencode == 4 or $gencode == 25)) {$tax = 'M' if $gencode == 4; $tax ='G' if $gencode == 25; return}
 ($gencode, $tax) = (11, 'B'); # Default
}

sub FindGenes {
 if (-f $circ) {open IN, $circ; $circ = 'C'; while (<IN>) {$circ = 'L', last if /^\S+\tlin/} close IN}
 mkdir 'trna'; chdir 'trna';
 RunCommand("perl $dir/bin/tfind.pl $tax .. &> tfind.log", "ttm.fa"); 
 mkdir '../protein'; chdir '../protein';
 my $kingdom = 'Bacteria'; $kingdom = 'Archaea' if $tax =~ /^A/;
 $nickname = '--locustag ' . $nickname if $nickname;
 RunCommand("prokka --rfam --prefix protein --locustag $dna --gcode $gencode --kingdom $kingdom --cpus $cpu --rnammer --notrna --outdir ./ --force --quiet $nickname ../$dna.fa", "protein.gff");
 unlink qw/protein.err protein.ffn protein.fna protein.fsa protein.gbk protein.sqn protein.tbl protein.txt/;
 RunCommand("$dir/bin/hmmsearch --domtblout protein.domtbl --cpu $cpu --cut_tc $hmm protein.faa &> /dev/null", "protein.domtbl");
 RunCommand("perl $dir/bin/integrase_finder.pl -faa protein.faa -gff protein.gff -integron $dir/db/intI_Cterm.hmm -pfam $hmm -xer $dir/db/xers.prf -intDoms $dir/db/integrase_domain_pfams.txt -cpu $cpu $extraDoms $force $verbose", "protein.phage");
 RunCommand("$dir/bin/hmmsearch --tblout tnp.tbl --noali --cut_tc --cpu $cpu $dir/db/TnpPred_HMM_Profiles.hmm protein.faa &> /dev/null", "tnp.tbl");
 RunCommand("$dir/bin/hmmsearch --domtbl is607.domtbl --noali --cut_ga --cpu $cpu $dir/db/is607.hmm protein.faa &> /dev/null", "is607.domtbl");
}

sub RunCommand {
 my ($command, $checkfile) = @_;
 if (not $force and ($checkfile and -e $checkfile)) {print "Skipping command: $command\n" if $verbose; return}
 print "Running command: $command\n" if $verbose;
 my $out = system $command;
 if ($out) {die "Command '$command' failed with error message $out\n"}
 else {print "Command '$command' succeeded\n"}
}

sub ReadFile {
 my @ret;
 open IN, $_[0] or die "Can't open $_[0]\n"; while (<IN>) {last if /^##FASTA/; next if /^#/; push @ret, $_} close IN;
 chomp @ret; return \@ret;
}

sub CollectGff { # Still in protein subdirectory
 my (%serials, %gff, %otherMob, %tnpCutoffs, %tnpPass, %tnpNicks, %is607);
 for (@{ReadFile("$dir/db/TnpPredCutoffs.txt")}) {$tnpCutoffs{$1} = $2 if /^(\S+)\t(\S+)/}
 for (@{ReadFile("$dir/db/TnpRename.txt")}) {$tnpNicks{$1} = $2 if /^(\S+)\t(\S+)/}
 for (@{ReadFile("protein.gff")}) {
  my @f = split "\t"; next if $f[2] =~ /^tm*RNA$/;
  s/inference=[^;]+;//; s/eC_number=[^;]+;//;
  my $id = '';
  if (s/locus_tag=([^;]+);//) {$id = $1}
  if ($id) {s/ID=[^;]+/ID=$id/} else {if (/ID=([^;]+)/) {$id = $1} else {$serials{NoID} ++; $id = "NoID.$serials{NoID}"}}
  %{$gff{$id}} = (dna => $f[0], L => $f[3], R => $f[4], line => $_, pfam => '', annot => '');
 }
 for (@{ReadFile("protein.phage")}) {next unless /ID=([^;]+)/; my $label = 'Y-Int'; $serials{$label} ++; $gff{$1}{annot} = "$label.$serials{$label}"}
 for (@{ReadFile("protein.xer")}) {next unless /ID=([^;]+)/; my $label = 'Xer'; $serials{$label} ++; $gff{$1}{annot} = "$label.$serials{$label}"}
 for (@{ReadFile("protein.integrons")}) {next unless /ID=([^;]+)/; my $label = 'Integron-Int'; $serials{$label} ++; $gff{$1}{annot} = "$label.$serials{$label}"}
 for (@{ReadFile("tnp.tbl")}) { # Hmmer table
  my @f = split /\s+/;
  next unless $f[2] and $tnpCutoffs{$f[2]} and $f[4] <= $tnpCutoffs{$f[2]};
  push @{$tnpPass{$f[0]}}, {tnp => $f[2], evalue => $f[4], score => $f[5]};
 }
 for (@{ReadFile("is607.domtbl")}) { # Hmmer domain table
  next if /^#/;
  my @f = split /\s+/;
  $is607{$f[0]} ++ if $f[13] > 64;  # Elsewhere cutoff given as 100
 }
 for my $prot (keys %tnpPass) {
  my $win = (sort {$$a{evalue} <=> $$b{evalue}} @{$tnpPass{$prot}})[0];
  $serials{$$win{tnp}} ++; 
  $gff{$prot}{annot} = "Tnp_$tnpNicks{$$win{tnp}}.$serials{$$win{tnp}}";
  my @g = split "\t", $gff{$prot}{line}; $g[5] = $$win{evalue}; $gff{$prot}{line} = join("\t", @g); # Include score
 }
 for (@{ReadFile("$dir/db/other_mobility.txt")}) {$otherMob{$1} = $2 if /(\S+)\t(\S+)/}
 for (@{ReadFile("protein.pfam.gff")}) {
  my @f = split "\t";
  next unless /ID=([^;]+);(pfam1=([^;]+).*)/;
  my ($id, $pfam, $label) = ($1, $2, $3);
  $gff{$id}{pfam} = $pfam;
  next if $gff{$id}{annot} or !$otherMob{$label};
  $label = 'Tnp_'.$label if $otherMob{$label} eq 'Tnp';
  if ($label eq 'Resolvase') {if (/pfam\d+=Recombinase/) {$label = 'S-Int'} else {$label = 'S-Core'}}
  $label = 'S-Int' if $label eq 'Recombinase'; # and /pfam\d+=Resolvase/;
  $label = 'S-Core_IS607' if $is607{$id};
  $serials{$label} ++; $gff{$id}{annot} = "$label.$serials{$label}";
 }
 chdir '..';
 for (@{ReadFile("trna/ttm.gff")}) {my @f = split "\t"; /ID=([^;]+)/; %{$gff{$1}} = (dna => $f[0], L => $f[3], R => $f[4], line => $_, pfam => '', annot => '')} 
 open OUT, ">$dna.gff";
 for my $gene (sort {$gff{$a}{dna} cmp $gff{$b}{dna} || $gff{$a}{L} <=> $gff{$b}{L} || $gff{$a}{R} <=> $gff{$b}{R}} keys %gff) {
  $gff{$gene}{annot} = "annot=$gff{$gene}{annot};" if $gff{$gene}{annot};
  $gff{$gene}{line} =~ s/;$//;
  print OUT "$gff{$gene}{line};$gff{$gene}{annot}$gff{$gene}{pfam}\n";
 }
 close OUT;
 print scalar(keys %gff), " gff lines\n";
}

sub Options {
my $version = '1.0 (Nov 2016)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] <GENOME.fa file>
  -tax:      Taxonomic info for query genome. Enter name of a file containing 
              NCBI taxonomy string, or use B for Bacteria, A for Archaea, M for 
              Mycoplasmatales/Entomoplasmatales, G for Gracilibacteria/candidate
              division SR1. Default: B.
  -gencode:  Genetic code table to use (see NCBI). Default: 11.
  -nickname: Brief name for genome (as might be used to start a locus_tag).
  -circle:   'C' if DNA is circular; 'L' if linear. Default: C.
  -force:    Overwrite current output files. Default: leave existing files.
  -outDir:   Output directory. Default: same directory as GENOME.fa.
  -extraDoms:Additional .domtbl files, comma-separated.
  -cpu:      Number of cpus to use. Default: 1.
  -hmm:      Main hmm database file. Default: $hmm
  Additional options: -help, -version, -verbose, -authors, -license

Example: perl $scriptname -verbose -gencode 11 -nickname Eco Eco.fa

END
my $authors = "AUTHORS: Kelly Williams (kpwilli\@sandia.gov), Corey Hudson, Britney Lau\n";
my $license = "LICENSE AND COPYRIGHT: Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.\n";

die $help if @ARGV == 0;
use Getopt::Long;
my $options_okay = GetOptions(
 'help' => sub {print $help; exit},
 'version' => sub {print "$scriptname version $version\n"; exit},
 'authors' => sub {print $authors; exit},
 'license' => sub {print $license; exit},
 'verbose' => sub {$verbose = "--verbose"},
 'force'   => sub {$force = "-force"},
 'gencode=i' => \$gencode,
 'tax=s' => \$tax,
 'nickname=s' => \$nickname,
 'extraDoms=s' => \$extraDoms,
 'hmm=s' => \$hmm,
 'circle=s' => \$circ,
 'cpu=i' => \$cpu,
);
die $help if !$options_okay;
}
