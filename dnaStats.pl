#!/usr/bin/perl
use strict; use warnings;

# OPTIONS
use Cwd 'abs_path';
my $dir = abs_path($0); $dir =~ s/\/([^\/]+)$//; my $scriptname = $1;
my ($verbose, $force, $complete, $virus, %virus) = (''. '', '', '');
my $circularity = 'C';
Options(); # see bottom of script; help and other messages, Getopt
if ($virus) {for (split ',', $virus) {$virus{$_} ++}}
my $dna = $ARGV[0]; $dna =~ s/\.fn*a$//;

# PROGRAM
my ($entry, $seq, $tax, %stats, %gff, %dom, %trna, %trna2dna); 
MoveDnas();
Tax();
for (keys %stats) {$stats{$_}{replicon} = Replicon($_)};
Stats();
Gffs();
Trna();
for (keys %stats) {PerDna($_) unless $_ eq 'all'}

# SUBROUTINES
sub MoveDnas {
 for (@{ReadFile("$dna.fa")}) {
  if (/^>(\S+)/) {
   WriteDna($entry, $seq) if $entry; $seq = ''; $entry = $1;
   $trna{$entry} = ''; $gff{$entry} = ''; $dom{$entry} = '';
  } else {$seq .= $_} 
 }
 WriteDna($entry, $seq) if $entry;
}

sub ReadFile {
 my @ret;
 print "Reading file $_[0]\n" if $verbose;
 open IN, $_[0] or die "Can't open $_[0]\n";
 while (<IN>) {last if /^##FASTA/; next if /^#/; push @ret, $_}
 close IN; chomp @ret; return \@ret;
}

sub WriteDna {
 my ($entry, $seq) = @_;
 return unless $seq;
 $stats{$entry}{len} = length $seq;
 mkdir $entry;
 open OUT, ">$entry/in.fa"; print OUT ">$entry\n$seq\n"; close OUT;
}

sub Tax {
 print "Reading taxonomy file\n" if $verbose;
 $tax = 'O';
 if (-f "$dna.tax") {
  my $text = do { local( @ARGV, $/ ) = "$dna.tax" ; <> } ;
  if ($text =~ /Viruses;| viruses;|Bacteriophages;/) {$tax = 'V'}
  elsif ($text =~ /Archaea;/) {$tax = 'A'}
  elsif ($text =~ /Bacteria;/) {$tax = 'B'}
 }
 if ($tax =~ s/^O/B/) {print "Taxonomy $tax not resolveable; using default, Bacteria\n" if $verbose}
}

sub Replicon {
 return 'Vir' if $tax eq 'V' or $virus{$_[0]};
 return 'Scf' unless $complete;
 return 'Chr' if keys %stats == 1;
 return 'Plm' if $stats{$_[0]}{len} <= 278000; 
 return 'Chr' if $stats{$_[0]}{len} >= 2580000; 
 return 'Int';
}

sub Stats {
 if (not $force and -f "$dna.stats") {for (@{ReadFile("$dna.stats")}) {$stats{$1}{line} = $_ if /^(\S+)/} return}
 my (%hskpEnrich, %fornEnrich);
 RunCommand("perl $dir/bin/relAbun.pl $dna.fa > $dna.dnaRA", "$dna.dnaRA");
 for (@{ReadFile("$dna.dnaRA")}) {next unless s/(\S+)\t//; $stats{$1}{ra} = $_}
 Circularity();
 for (@{ReadFile("$dir/db/housekeep_enrich.txt")}) {$hskpEnrich{$1} = $2 if /(\S+)\t(\S+)/}
 for (@{ReadFile("$dir/db/foreign_enrich.txt")})   {$fornEnrich{$1} = $2 if /(\S+)\t(\S+)/}
 for (@{ReadFile("$dna.gff")}) {
  next unless /^(\S+)/; $entry = $1;
  $gff{$entry} .= $_ . "\n";
  $trna2dna{$1} = $entry if /\ttm*RNA.*ID=([^;]+)/;
  next unless /\tCDS\t/; $stats{$entry}{cds} ++; 
  next unless /pfam1=([^;]+)/;   $stats{$entry}{pfam} ++;
  $stats{$entry}{hskp} += $hskpEnrich{$1} if $hskpEnrich{$1};
  $stats{$entry}{forn} += $fornEnrich{$1} if $fornEnrich{$1};
 }
 for my $entry (keys %stats) { # Includes 'all'
  for (qw/len pfam hypoth hskp forn cds/) {
   $stats{$entry}{$_} = 0 unless $stats{$entry}{$_};
   $stats{all}{$_} += $stats{$entry}{$_} unless $entry eq 'all';
  }
 }
 open OUT, ">$dna.stats";
 for my $entry (sort {$stats{$b}{len} <=> $stats{$a}{len}} keys %stats) {
  $stats{$entry}{hypoth} = ($stats{$entry}{cds}-$stats{$entry}{pfam})/$stats{$entry}{cds} if $stats{$entry}{cds};
  for (qw/hskp forn/) {$stats{$entry}{$_} /= $stats{$entry}{pfam} if $stats{$entry}{pfam}}
  $stats{$entry}{replicon} = 'NA' if $entry eq 'all';;
  $stats{$entry}{line} = $entry; for (qw/len circ replicon cds pfam hypoth hskp forn ra/) {$stats{$entry}{line} .= "\t$stats{$entry}{$_}"}
  print OUT "$stats{$entry}{line}\n";
 }
 close OUT;
}

sub RunCommand {
 my ($command, $checkfile) = @_;
 if (not $force and ($checkfile and -e $checkfile)) {print "Skipping command: $command\n" if $verbose; return}
 print "Running command: $command\n" if $verbose;
 my $out = system($command);
 die "Command '$command' failed with error message $out\n" if $out;
}

sub Circularity {
 if ($circularity =~ /^[lL]$/) {for (keys %stats) {$stats{$_}{circ} = 'Lin'} return}
 elsif (not -f $circularity) {for (keys %stats) {$stats{$_}{circ} = 'Cir'} return}
 for (@{ReadFile($circularity)}) {
  next unless /^(\S+)\t([CL])/i and $stats{$1};
  $stats{$1}{circ} = uc $2;
 }
 for (keys %stats) {$stats{$_}{circ} = 'Cir' unless $stats{$_}{circ}}
}

sub Gffs {
 my $gfftest = 0; for (keys %gff) {$gfftest ++ if $gff{$_}}
 unless ($gfftest) { # May have already collected during sub Stats
  for (@{ReadFile("$dna.gff")}) {$gff{$1} .= $_ . "\n" if /^(\S+)/; $trna2dna{$2} = $1 if /^(\S+).*\ttm*RNA.*ID=([^;]+)/}
 }
 for (@{ReadFile("protein/protein.domains.gff")}) {$dom{$1} .= $_ . "\n" if /^(\S+)/}
}

sub Trna {
 for (@{ReadFile("trna/ttm.fa")}) {
  $entry = $1 if /^>(\S+)/; 
  die scalar(keys %trna2dna), " trna2dna keys but none for $entry\n" unless $trna2dna{$entry};
  $trna{$trna2dna{$entry}} .= $_ . "\n";
 }
}

sub PerDna {
 $entry = $_[0];
 open OUT, ">$entry/in.stats"; print OUT $stats{$entry}{line}; close OUT;
 open OUT, ">$entry/in.gff";   print OUT $gff{$entry};         close OUT;
 open OUT, ">$entry/in.domains.gff"; print OUT $dom{$entry};   close OUT;
 open OUT, ">$entry/in.trna.fa";     print OUT $trna{$entry};  close OUT;
}

sub Options {
my $version = '1.0 (Nov 2016)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] GENOME_FASTA_FILE
  -circle:   Specify C if all genomic DNA sequences are circular, L if all DNAs 
              are linear, or a filename for a tab-delimited file of query and 
              circularity (eg. acc.vers[tab]circular/linear). Default: C.
  -virus:    Comma-separated list of entries assigned as viruses.
  -complete: Consider genome complete and categorize replicons. Default: 
              consider genome incomplete and call all entries contigs.
  -force:    Overwrite current output files. Default: leave existing files.
  Additional options: -help, -version, -verbose, -authors, -license

Example: perl $scriptname -circle L -verbose GENOME_FASTA_FILE
Note: run after tater.pl, from same folder as GENOME_FASTA_FILE 

END
my $authors = "AUTHORS: Kelly Williams (kpwilli\@sandia.gov), Corey Hudson, Britney Lau, Owen Solberg\n";
my $license = "LICENSE AND COPYRIGHT: Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.\n";

die $help if @ARGV == 0;
use Getopt::Long;
my $options_okay = GetOptions(
 'help'     => sub {print $help; exit},
 'version'  => sub {print "$scriptname version $version\n"; exit},
 'authors'  => sub {print $authors; exit},
 'license'  => sub {print $license; exit},
 'verbose'  => sub {$verbose = "--verbose"},
 'virus=s'  => \$virus,
 'force'    => \$force,
 'circle=s' => \$circularity,
 'complete' => \$complete,
);
die $help if !$options_okay;
}
